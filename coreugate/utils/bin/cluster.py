#!/usr/bin/env python3
import numpy,sys, pathlib, pandas
from sklearn.cluster import AgglomerativeClustering



def extract_unclustered(df, col, col1):

    u = df.groupby(col)[col1].agg('count').reset_index()
    
    uc = list(u[u[col1] == 1][col])
    
    return uc


thresholds = sys.argv[2].split(',')
# get distances into lists
tabfile = sys.argv[1]
if pathlib.Path(tabfile).exists():
    isos = []
    mat = []
    if isinstance(tabfile, str):
        print(f"opening {tabfile}")
        with open(tabfile, 'r') as d:
            tab = d.read().split('\n')
            
            # get isolates
            for line in tab[1:len(tab)]:
                ln = line.split('\t')
                if len(ln) > 1:
                    # isos.append(ln[0])
                    mat.append(ln[1:len(ln)])
       
            isos = tab[0].split('\t')[1:]
    
    clusters = pandas.DataFrame()
    # convery the lists to a numpy array
    print("converting matrix to numpy array")
    X = numpy.array(mat)
    X = X.astype(numpy.float64)
    
    for level in thresholds:
        print(f"clustering at {level}")
        clustering = AgglomerativeClustering(n_clusters = None, affinity = 'precomputed',linkage = 'complete', distance_threshold =int(level)).fit(X)
        df = pandas.DataFrame(data = {'ID': isos, f"Tx:{level}": clustering.labels_})
        print(len(df[f"Tx:{level}"].unique()))
        df[f"Tx:{level}"] = df[f"Tx:{level}"] + 1
        df = df.fillna('')
        df[f"Tx:{level}"] = df[f"Tx:{level}"].apply(lambda x: f"{x}")
       
        
        uc = extract_unclustered(df = df, col = f"Tx:{level}", col1 = 'ID')
        df[f"Tx:{level}"] = numpy.where(df[f"Tx:{level}"].isin(uc), 'UC', df[f"Tx:{level}"])
        if clusters.empty:
            clusters = df
        else:
            clusters = pandas.merge(df, clusters)
       
    clusters = clusters.fillna('')

    clusters.to_csv('clusters.txt', sep = '\t', index = False)