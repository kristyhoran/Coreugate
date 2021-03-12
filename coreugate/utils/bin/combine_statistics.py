#!/usr/bin/env python3
import pandas, sys, ast, pathlib

publish_dir = pathlib.Path(sys.argv[1])
threshold = sys.argv[2]
# print(sys.argv[2:])
statistics = []
for i in sys.argv[3:]:
    l = i.strip('[,]')
    # print(l)
    statistics.append(l)



if publish_dir.exists():
    df = pandas.read_csv(f"{publish_dir}", sep = '\t')
else:
    df = pandas.DataFrame()

for s in statistics:
    if pathlib.Path(s).exists():
        tmp = pandas.read_csv(s, sep = '\t')

        if df.empty:
            df = tmp
        else:
            df = df.append(tmp)
df['ID'] = df['Genome'].str.replace("\.fa([a-z])*", "", regex = True)
df['TOTAL'] = df[[ "EXC","INF","LNF","PLOT","NIPH","ALM","ASM"]].apply(lambda x: sum(x), axis = 1)
df['PASSED'] = df[["EXC","INF","TOTAL"]].apply(lambda x:sum(x[:2])/x[2], axis = 1)
df = df[['ID', 'Genome',"EXC","INF","LNF","PLOT","NIPH","ALM","ASM", "TOTAL", 'PASSED']]
passed = list(df[df['PASSED'] > float(threshold)]['ID'])

pathlib.Path('passed.txt').write_text('\n'.join(passed))
df.to_csv(f'overall_statistics.txt', sep = '\t', index = False)