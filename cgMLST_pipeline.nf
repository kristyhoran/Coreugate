Channel
	.fromFilePairs( params.fastq )
	.ifEmpty { error "Can't find any reads matching : ${ params.fastq }" }
	.set { fastq_input }

s = params.fastq
fq = s.split("/")[-2]

process assembly {

	maxForks 2

	publishDir 'results/assemblies/', pattern: "*_contigs.fa"

	input:
	set fastq_id, file( fastq_pair ) from fastq_input

	output:
	set fastq_id, "${ fastq_id }_*_contigs.fa" into fastq_out
	
		workflow.onComplete {
 	  		 println "Assembly of $fq using $params.assembler completed at: $workflow.complete"
	 				 
			    }


	script:
	if( params.assembler == 'skesa' )           
        """
        skesa --fastq ${fastq_pair[0]},${fastq_pair[1]} --gz --use_paired_ends --cores 4 > ${fastq_id}_skesa_contigs.fa
        """
	else if( params.assembler == 'shovill' )
	"""
        shovill --outdir results/ --R1 ${fastq_pair[0]} --R2 ${fastq_pair[1]} --force
        cp results/contigs.fa ${fastq_id}_shovill_contigs.fa
        """
		
	
}
process chewie {
        
        maxForks 2

        publishDir 'results/chewie/'

        input:
        set fastq_id, file(assembly) from fastq_out

        
        output:
        set fastq_id, "${fastq_id}_*_results_alleles.tsv" into chewie_alleles_out
        set fastq_id, "${fastq_id}_*_results_statistics.tsv" into chewie_statistics_out

	workflow.onComplete {
                         println "chewBBACA completed at: $workflow.complete"
                         
                            }	

        """
        if [ ! -d genomes ]; then
                mkdir genomes
        fi
        
        cp -r ${workflow.projectDir}/${params.schemaPath} .

        mv ${fastq_id}_*_contigs.fa genomes/
        chewBBACA.py AlleleCall -i ./genomes/ -g ${params.schemaPath} -o chewie_results/ --cpu 16 --fc
        cp chewie_results/*/results_alleles.tsv ${fastq_id}_${params.assembler}_results_alleles.tsv
        cp chewie_results/*/results_statistics.tsv ${fastq_id}_${params.assembler}_results_statistics.tsv
        
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

