
import os, subprocess

VERSION = 2.0
HEADER = 'COREugate v ' + str(VERSION)
print(HEADER.center(80, '='))
print('A snakemake workflow for cgMLST'.center(80))
print('type `snakemake help` or https://github.com/kristyhoran/Coreugate'.center(80))
print('All issues can be filed at https://github.com/kristyhoran/Coreugate')
print('Good luck! May the force be with you.')

# set path to config file
configfile:'config.yaml'
# configs variables
# path to schema
SCHEMA_PATH = config['schemaPath']

# a tab separated file with col 1 = isolate id, col 2 = path to read 1, and col 3 = path to read 2
INPUT = config['isolate_file']
# ASSEMBLER
ASSEMBLER = config['assembler']
# CPU
prepCPU = config['prepCPU']
chewieCPU = config['chewieCPU']
# get current directory
cwd = os.getcwd()

# a list to place samples for wildcards
SAMPLES = []

# generate the output folder and symlinks to reads
# TODO try to make a rule for this
# # TODO: add logic preventing overwriting existing folder

for i in open(INPUT, 'r'):
    s = i.split()[0]
    r1 = i.split()[1]
    r2 = i.split()[2]
    SAMPLES.append(s)
    if not os.path.exists(cwd + '/coreugate/READS/' + s):
        print('Creating READS directory')
        os.makedirs(cwd + '/coreugate/READS/' + s)
    if not os.path.exists(cwd + '/coreugate/READS/'+ s + '/R1.fq.gz'):
        print('Creating symlink to ' + s + ' R1')
        subprocess.run(['ln', '-f', '-s', r1, cwd + '/coreugate/READS/'+ s+ '/R1.fq.gz'])
    if not os.path.exists(cwd + '/coreugate/READS/'+ s + '/R2.fq.gz'):
        print('Creating symlink to ' + s + ' R2')
        subprocess.run(['ln', '-f', '-s', r2, cwd + '/coreugate/READS/'+ s + '/R2.fq.gz'])


rule all: #
    input:
        'coreugate/overall_statistics.tsv', 'coreugate/overall_alleles.tsv'

rule help:
    """
    Coreugate v 2.0
    Thank you for using Coreugate for your cgMLST needs!
    All issues can be filed at https://github.com/kristyhoran/Coreugate
    Parameters can be added after the --config flag
    example
    --config schemaPath=path/to/schema
    """
    run:
        for rule in workflow.rules:
            print(rule.name)
            print(rule.docstring)


rule assembly:
    """
    Assemble genome using the chosen assembler
    assembler=str options: skesa, shovill (default), spades (version 3.7)
    """
    input:
        r1 = "coreugate/READS/{sample}/R1.fq.gz",
        r2 = "coreugate/READS/{sample}/R2.fq.gz"
    output:
        "coreugate/assemblies/{sample}.fa"
    singularity:
        "shub://phgenomics-singularity/multi_assembler_singularity@latest"

    shell:
        """
        if [ "{ASSEMBLER}" == "shovill" ]; then

            shovill --outdir temp --R1 {input.r1} --R2 {input.r2} --force --minlen 200
            cp temp/contigs.fa {output}
            rm -r temp

        elif [ "{ASSEMBLER}" == "skesa" ]; then
            skesa --fastq {input.r1},{input.r2} --gz --use_paired_ends --cores 4 > {output}

        elif [ "{ASSEMBLER}" == "spades" ]; then
            spades.py -o temp/ -1 {input.r1} -2 {input.r2}
            cp temp/contigs.fasta {output}
            rm -r temp
        fi
        """

rule assembly_qc:
    """
    min_contigs=int options: the min number of contigs desired (default: off)
    min_contig_size=int options: remove small contigs (default: 500)
    """

    input:
        "coreugate/assemblies/{sample}.fa"
    output:
        "coreugate/assemblies/{sample}.qc.fa"
    message:
        'Performing quality control of assemblies.'
    script:
        'bin/qc.py'



rule chewBBACA:
    """
    Prepare external schema (if needed) and perform AlleleCall
    schemaPath=str COMPULSORY path the folder containing cgMLST schema
    prepCPU=int option: cpus used for PrepExternalSchema (default: 36)
    chewieCPU=int option: cpus used for AlleleCall (default: 36)
    """
    input:
        isolate ="coreugate/assemblies/{sample}.qc.fa",
        schema = SCHEMA_PATH
    output:
        directory('coreugate/chewBBACA/{sample}')
    resources:
        wget_limit=1
    singularity:
        "docker://mickaelsilva/chewbbaca_py3"
    shell:
        """
            if [ ! -d {SCHEMA_PATH}/short ]; then
                chewBBACA.py PrepExternalSchema -i {input.schema}/ -cpu {prepCPU} -v
            fi

            rm -rf {SCHEMA_PATH}/temp
            chewBBACA.py AlleleCall -i {input.isolate} -g {input.schema}/ -o chew_temp --cpu {chewieCPU} --fc
            mv chew_temp/*/ {output}
            rm -r chew_temp
        """


rule combinechewie:
    """
    No need anything here!
    """
    input:
        expand("coreugate/chewBBACA/{sample}", sample = SAMPLES)
    output:
        'coreugate/overall_statistics.tsv', 'coreugate/overall_alleles.tsv'
    message:
        'Combining chewBBACA results.'
    script:
        'bin/combine_results.py'
