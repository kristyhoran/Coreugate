#!/usr/bin/env nextflow
nextflow.enable.dsl=2   

params.schema_path = file('/home/khhor/cgMLST/dev/vre_prepped')
params.ptf = file('/home/khhor/cgMLST/dev/vre_prepped/Enterococcus_faecium.trn')
params.threads = 16
params.force = true
params.profile_pass = 0.95
params.cluster = true
params.thresholds = [20,50]

contigs = Channel.fromPath('CONTIGS/*/*.fa')
            .map { def file -> tuple(file.baseName, file)}
// println contigs.view()

process PROFILES {
    publishDir "$task.process", mode: 'copy'
    maxForks 1
    cache 'lenient'
    cpus params.threads
    input:
        tuple val(sample_id), path(contig)
    output:
        tuple val(sample_id), path("*_results_alleles.tsv"), emit: alleles
        tuple val(sample_id), path("*_results_statistics.tsv"), emit: statistics

    script:
        """
        [ !  -f ptf.trn ] && ln -s ${params.ptf} ptf.trn
        chewie AlleleCall -i . -g $params.schema_path --ptf ptf.trn -o ${sample_id}_profile --cpu $task.cpus --fr
        cp ${sample_id}_profile/results_*/results_alleles.tsv ${sample_id}_results_alleles.tsv
        cp ${sample_id}_profile/results_*/results_statistics.tsv ${sample_id}_results_statistics.tsv
        rm -r ${sample_id}_profile
        """

}


process COLLATE_STATS {
    publishDir "$task.process", mode: 'copy'

    input:
        val(statistics)
    output:
        path('overall_statistics.txt'), emit: overall_statistics
        path('passed.txt'), emit:passed_profile
    script:
        """
        combine_statistics.py $task.process $params.profile_pass $statistics 
        """
}


process COLLATE_ALLELES {
    publishDir "$task.process", mode: 'copy'

    input:
        val(alleles)
        path passed_profile
    output:
        path('overall_alleles.txt'), emit: overall_alleles
    script:
        """
        combine_alleles.py $task.process $passed_profile $alleles
        """
}

process PAIRWISE_DISTANCE {
    publishDir "$task.process", mode: 'copy'

    input:
        path overall_alleles
    output:
        path('pad.txt'), emit: dists
    script:
        """
        cgmlst-dists $overall_alleles > pad.txt
        """

}

process CLUSTER {
    publishDir "$task.process", mode: 'copy'

    input:
        path pad
    output:
        path('clusters.txt'), emit: clusters
    script:
        """
        cluster.py $pad $params.thresholds
        """
}

workflow {
    PROFILES( contigs )
    COLLATE_STATS ( PROFILES.out
                                .statistics
                                .map { sample_id, file -> file}
                                .collect() )
    // collate alleles into a table.. only including ones that have more than the threshold set ie 0.95
    COLLATE_ALLELES ( PROFILES.out
                                .alleles
                                .map { sample_id, file -> file}
                                .collect(),
                      COLLATE_STATS.out.passed_profile)
    PAIRWISE_DISTANCE ( COLLATE_ALLELES.out.overall_alleles )
    
    if( params.cluster ) {
        CLUSTER ( PAIRWISE_DISTANCE.out.dists)
    }
}