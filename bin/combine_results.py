import os

for i in snakemake.input:
    for f in os.listdir(i + '/'):
        in_path = i + '/' + f
        if f == 'results_statistics.tsv':
            stat = open(in_path, 'r')

            if not os.path.exists(snakemake.output[0]):
                with open(snakemake.output[0], 'w') as out:
                    for l in stat.readlines():
                        out.write(l)
                    out.write('\n')
            else:
                with open(snakemake.output[0], 'a') as out:
                    out.write(stat.readlines()[1])
                    out.write('\n')
        if f == 'results_alleles.tsv':
            allele = open(in_path, 'r')

            if not os.path.exists(snakemake.output[1]):
                with open(snakemake.output[1], 'w') as out:
                    for l in allele.readlines():
                        out.write(l)
                    out.write('\n')
            else:
                with open(snakemake.output[1], 'a') as out:
                    out.write(allele.readlines()[1])
                    out.write('\n')
