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
    print(len(mat))
    clusters = pandas.DataFrame()
    # convery the lists to a numpy array
    print("converting matrix to numpy array")
    X = numpy.array(mat)
    X = X.astype(numpy.float64)
    
    for level in thresholds:
        print(f"Clustering at {level}")
        clustering = AgglomerativeClustering(n_clusters = None, affinity = 'precomputed',linkage = 'complete', distance_threshold =int(level)).fit(X)
        df = pandas.DataFrame(data = {'ID': isos, f"Tx:{level}": clustering.labels_})
        print(len(clustering.labels_))
        print(len(isos))
        print(df.shape)
        
        print(len(df[f"Tx:{level}"].unique()))
        df[f"Tx:{level}"] = df[f"Tx:{level}"] + 1
        df = df.fillna('')
        df[f"Tx:{level}"] = df[f"Tx:{level}"].apply(lambda x: f"{x}")
        print(df.shape)
        
        uc = extract_unclustered(df = df, col = f"Tx:{level}", col1 = 'ID')
        df[f"Tx:{level}"] = numpy.where(df[f"Tx:{level}"].isin(uc), 'UC', df[f"Tx:{level}"])
        df = df.drop_duplicates(subset=['ID'])
        if clusters.empty:
            clusters = df
        else:
            clusters = df.merge(clusters, right_on = 'ID',left_on ='ID', how = 'inner')
        clusters = clusters.drop_duplicates()
        print(clusters.shape)
    clusters = clusters.fillna('')

    clusters.to_csv('clusters.txt', sep = '\t', index = False)