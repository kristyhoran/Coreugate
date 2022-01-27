#!/usr/bin/env nextflow
nextflow.enable.dsl=2   

// params.schema_path = file(params.schema_path)
// params.ptf = file(params.ptf)


contigs = Channel.fromPath(params.input)
// println contigs.view()

process PROFILES {
    publishDir "${task.process}", mode: 'copy'
    maxForks 1
    cache 'lenient'
    cpus params.threads
    input:
        path(contigs)
    output:
        path("results_alleles.tsv"), emit: alleles
        path("results_statistics.tsv"), emit: statistics

    script:
        """
        [ !  -f ptf.trn ] && ln -s ${params.ptf} ptf.trn
        chewie AlleleCall -i $contigs -g $params.schema_path --ptf ${params.ptf} -o chewie_profile --cpu $task.cpus --fr
        cp chewie_profile/results_*/results_alleles.tsv results_alleles.tsv
        cp chewie_profile/results_*/results_statistics.tsv results_statistics.tsv
        rm -r chewie_profile
        """

}


process COLLATE_STATS {
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
        val(statistics)
    output:
        path('overall_statistics.txt'), emit: overall_statistics
        path('passed.txt'), emit:passed_profile
    script:
        
        """
        combine_statistics.py ${launchDir}/overall_statistics.txt $params.profile_pass $statistics 
        """
}


process COLLATE_ALLELES {
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
        val(alleles)
        path passed_profile
    output:
        path('overall_alleles.txt'), emit: overall_alleles
    script:
        """
        combine_alleles.py "${params.publish_dir}"/overall_alleles.txt $passed_profile $alleles
        """
}

process PAIRWISE_DISTANCE {
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
        path overall_alleles
    output:
        path('pad.txt'), emit: dists
    script:
        // def lines = file(overall_alleles).countLines() - 1
        """
        x=\$(csvtk nrows $overall_alleles)
        cgmlst-dists $overall_alleles | sed "s/cgmlst-dists/\$x/g" > pad.txt
        """

}

process CLUSTER {
    publishDir "${launchDir}", mode: 'copy'

    input:
        path pad
    output:
        path('clusters.txt'), emit: clusters
    script:
        """
        cluster.py $pad $params.thresholds
        """
}

process RAPIDNJ {
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
        path pad
    output:
        path('pad.nwk'), emit: tree
    script:
        """
        rapidnj -i pd -o t -x pad.nwk $pad
        """
}


workflow {
    PROFILES( contigs )
    COLLATE_STATS ( PROFILES.out
                                .statistics
                                .collect() )
    // collate alleles into a table.. only including ones that have more than the threshold set ie 0.95
    COLLATE_ALLELES ( PROFILES.out
                                .alleles
                                .collect(),
                      COLLATE_STATS.out.passed_profile)
    PAIRWISE_DISTANCE ( COLLATE_ALLELES.out.overall_alleles )
    RAPIDNJ ( PAIRWISE_DISTANCE.out.dists )
    if( params.cluster ) {
        CLUSTER ( PAIRWISE_DISTANCE.out.dists)
    }
}