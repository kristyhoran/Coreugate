import os
import pathlib
import sys

input_file = sys.argv[1]

out_alleles = pathlib.Path('overall_alleles.tsv')
out_stats = pathlib.Path('overall_statistics.tsv')
with open(input_file, 'r') as f:
    for i in f.readlines():
        l = i.strip()
        p = pathlib.Path(l)
        for f in p.iterdir():
            if f.name == 'results_statistics.tsv':
                stat = open(f, 'r')
                if not out_stats.exists():
                    with open(out_stats, 'w') as out:
                        for l in stat.readlines():
                            out.write(l)
                        out.write('\n')
                else:
                    with open(out_stats, 'a') as out:
                        out.write(stat.readlines()[1])
                        out.write('\n')
            if f.name == 'results_alleles.tsv':
                allele = open(f, 'r')
                if not os.path.exists(out_alleles):
                    with open(out_alleles, 'w') as out:
                        for l in allele.readlines():
                            out.write(l)
                        out.write('\n')
                else:
                    with open(out_alleles, 'a') as out:
                        out.write(allele.readlines()[1])
                        out.write('\n')
