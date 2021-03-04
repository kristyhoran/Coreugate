#!/usr/bin/env nextflow
nextflow.enable.dsl=2   
params.cluster = false
params.schema_path = file('/home/khhor/cgMLST/dev/vre_prepped')
params.ptf = file('/home/khhor/cgMLST/dev/vre_prepped/Enterococcus_faecium.trn')
params.threads = 36
contigs = Channel.fromPath('CONTIGS/*/*.fa')
            .map { file -> tuple(file.baseName, file)}
// println contigs.view()

process PROFILES {
    publishDir "$task.process", mode: 'copy'
    maxForks 1
    cpus params.threads
    input:
        tuple val(sample_id), path(contig)
    output:
        path("${sample_id}_results_alleles.tsv"), emit: alleles
        path("${sample_id}_results_statistics.tsv"), emit: statistics

    script:
        """
        [ !  -f ptf.trn ] && ln -s ${params.ptf} ptf.trn
        chewie AlleleCall -i . -g $params.schema_path --ptf ptf.trn -o ${sample_id}_profile --cpu $task.cpus --fr
        cp ${sample_id}_profile/results_*/results_alleles.tsv ${sample_id}_results_alleles.tsv
        cp ${sample_id}_profile/results_*/results_statistics.tsv ${sample_id}_results_statistics.tsv
        """

}

process COLLATE_STATS {
    publishDir "$task.process", mode: 'copy'

    input:
        statistics
    output:
        path('overall_statistics.txt'), emit: overall_statistics
    script:
        """

        """
}

workflow {
    PROFILES( contigs )
    PROFILES.out
        .alleles.collect()
        .view()
}