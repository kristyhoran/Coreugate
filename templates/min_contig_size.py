#!/usr/bin/env python3

from Bio import SeqIO

ASSEMBLY = "$assembly"
MIN_SIZE = $contig_size
MIN_CONTIGS = $min_contigs
#ASSEMBLY = "ERR1102467.fa"
#MIN_SIZE = 500
#MIN_CONTIGS = 0
NEW_NAME = ASSEMBLY[:-2] + "QCd.fa"
#print(NEW_NAME, ASSEMBLY)
records = list(SeqIO.parse(ASSEMBLY, "fasta"))
no_contigs = len(records)
if MIN_CONTIGS == 0:
    MIN = no_contigs + 2
else:
    MIN = MIN_CONTIGS
out_put_list = []
print(MIN, no_contigs)
if no_contigs < MIN:
    print(MIN)
    for record in records:
        if len(record.seq) > MIN_SIZE:
            out_put_list.append(record)
  
#print("Found %i sequences" % len(out_put_list))
if len(out_put_list) > 0:
    SeqIO.write(out_put_list, NEW_NAME, "fasta")
