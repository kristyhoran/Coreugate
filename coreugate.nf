// set header for running COREugate"
coreugate_version = "1.0"
log.info "".center(60, "=")
log.info "COREugate".center(60)
log.info "Version ${coreugate_version}".center(60)
log.info "".center(60, "=")


println params.assemble


if (params.fasta != 'none') {
	
	qc_input = Channel.fromPath( 'assembly/*' )
	log.info "Set input channel".center(60)
			}
if (params.assemble == 'enable') {

	Channel
	        .fromFilePairs( params.fastq )
        	.ifEmpty { error "Can't find any reads matching : ${ params.fastq }" }
	        .set { fastq_input }
	log.info "Set input channel for assembly".center(60)
}




if (params.assemble == 'enable') {

process assembly {

	maxForks 2

	publishDir 'results/assemblies/', pattern: "*.fa"

	//when:
        //params.assemble == 'enable'	

	input:
	set sample_id, file( fastq_pair ) from fastq_input

	output:
	set sample_id, "${ sample_id }.fa" into qc_input
			    
	script:
	if( params.assembler == 'skesa' )           
        """
        skesa --fastq ${fastq_pair[0]},${fastq_pair[1]} --gz --use_paired_ends --cores 4 > ${sample_id}.fa
        """
	else if( params.assembler == 'shovill' )
	"""
        shovill --outdir results/ --R1 ${fastq_pair[0]} --R2 ${fastq_pair[1]} --force --minlen 200
        cp results/contigs.fa ${sample_id}.fa
        """
	else if( params.assembler == 'spades' )
	
	"""
	spades.py -o results/ -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} 
	cp results/contigs.fasta ${sample_id}.fa	

	"""
		
	
}
}
// add quality control step only use assemblies with less than a desired threshold default 85
min_contigs = params.contigs
contig_size = params.contig_size
log.info "QC of assemblies excluding isolates with greater than $min_contigs"
process assemblyqc {
	publishDir 'results/assemblies_qc/', pattern: "*.fa"
	
	
	input:
	set sample_id, file(assembly) from qc_input

	output:
	set sample_id, "${sample_id}.fa" into qc_output


	script:
	template "min_contig_size.py"
	
}
 
process chewie {
        
        maxForks 2

        publishDir 'results/chewie/'

        input:
        set sample_id, file(assembly) from qc_output

        
        output:
        set sample_id, "${sample_id}_*_results_alleles.tsv" into chewie_alleles_out
        set sample_id, "${sample_id}_*_results_statistics.tsv" into chewie_statistics_out

	workflow.onComplete {
                         println "chewBBACA completed at: $workflow.complete"
                         
                            }	

        """
        if [ ! -d genomes ]; then
                mkdir genomes
        fi
        
        cp -r ${workflow.projectDir}/${params.schemaPath} .

        mv ${sample_id}.fa genomes/
        chewBBACA.py AlleleCall -i ./genomes/ -g ${params.schemaPath} -o chewie_results/ --cpu 16 --fc
        cp chewie_results/*/results_alleles.tsv ${sample_id}_${params.assembler}_results_alleles.tsv
        cp chewie_results/*/results_statistics.tsv ${sample_id}_${params.assembler}_results_statistics.tsv
        
        """ 
}

process combine_stats {

        publishDir "results/combine/", overwrite: true

        input:
        file data from chewie_statistics_out.collect()


        output:
        file "overall_statistics.tsv"

        """
        
        for i in \$(ls *_results_*.tsv);
        do
                echo \$(sed -n '1p' \${i}) > overall_statistics.tsv
        done 

        for i in \$(ls *_results_*.tsv);
        do
                echo \$(sed -n '2p' \${i})
                
        done >> overall_statistics.tsv

        """     

        
}      

process combine_alleles {

        publishDir "results/combine/", overwrite: true

        input:
        file data from chewie_alleles_out.collect()


        output:
        file "overall_alleles.tsv"

        """

        for i in \$(ls *_results_*.tsv);
        do
                echo \$(sed -n '1p' \${i}) > overall_alleles.tsv
        done    
        
        for i in \$(ls *_results_*.tsv);
        do      
                echo \$(sed -n '2p' \${i})
                
        done >> overall_alleles.tsv

        """
}

