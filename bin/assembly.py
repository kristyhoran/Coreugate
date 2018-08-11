import subprocess


if snakemake.config['assembler'] == 'shovill':
    subprocess.run(['shovill', '--outdir ','temp','--R1', snakemake.input.r1, '--R2',snakemake.input.r2,'--force','--minlen', '200'])
    subprocess.run(['cp','temp/contigs.fa' ,snakemake.output])
    subprocess.run(['rm','-r' ,'temp'])
elif config['assembler'] == 'skesa':
    subprocess.run(['skesa', '--fastq', snakemake.input.r1,snakemake.input.r2, '--gz' ,'--use_paired_ends' ,'--cores' ,'4' ,'>',snakemake.output])
elif config['assembler'] == 'spades':
    subprocess.run(['spades.py','-o', 'temp/','-1', snakemake.input.r1,'-2',snakemake.input.r2])
    subprocess.run(['cp', 'temp/contigs.fasta', snakemake.output])
    subprocess.run(['rm','-r' ,'temp'])
