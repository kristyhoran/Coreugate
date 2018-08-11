from Bio import SeqIO

ASSEMBLY = str(snakemake.input)

# add logic for no min contig number
MIN_SIZE = snakemake.config['min_contig_size']
MIN_CONTIGS = snakemake.config['min_contigs']

records = list(SeqIO.parse(ASSEMBLY, "fasta"))
no_contigs = len(records)
if MIN_CONTIGS == 0:
    MIN = no_contigs + 2
else:
    MIN = MIN_CONTIGS

if no_contigs < MIN:
    for record in records:

        if len(record.seq) > MIN_SIZE:
            with open(str(snakemake.output), 'a') as out:
                out.write('>'+str(record.description) + '\n')
                out.write(str(record.seq) + '\n')
