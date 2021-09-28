#!/usr/bin/env python3
import pandas, sys, ast, pathlib

output_file = pathlib.Path(sys.argv[1])
# print(sys.argv[2:])
if pathlib.Path(sys.argv[2]).exists:
    with open(sys.argv[2], 'r') as ps:
        passed = ps.read().split('\n')
alleles = []
for i in sys.argv[3:]:
    l = i.strip('[,]')
    # print(l)
    alleles.append(l)


if output_file.exists():
    df = pandas.read_csv(f"{output_file}", sep = '\t', dtype = str)
else:
    df = pandas.DataFrame()
print(passed)
for a in alleles:
    if pathlib.Path(a).exists():
        print(a)
        tmp = pandas.read_csv(a, sep = '\t', dtype = str)
        tmp['FILE'] = tmp['FILE'].str.replace("\.fa([a-z])*", "", regex = True)
        print(tmp)
        tmp = tmp[tmp['FILE'].isin(passed)]
        # if tmp.at[0,'FILE'] in passed:
            
        if df.empty:
            df = tmp
        else:
            df = df.append(tmp).drop_duplicates()



df.to_csv('overall_alleles.txt', sep = '\t', index = False)