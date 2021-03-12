#!/usr/bin/env python3
import numpy,sys, pathlib, pandas
from sklearn.cluster import AgglomerativeClustering



def extract_unclustered(df, col, col1):

    u = df.groupby(col)[col1].agg('count').reset_index()
    # print(u)
    uc = list(u[u[col1] == 1][col])
    # print(uc)
    return uc

# get thresholds
thresholds = [sys.argv[2].split(',')]
# for i in sys.argv[2:]:
#     l = i.strip('[,]')
#     # print(l)
#     thresholds.append(l)
# print(thresholds)


# get distances into lists
tabfile = sys.argv[1]
if pathlib.Path(tabfile).exists():
    isos = []
    mat = []
    if isinstance(tabfile, str):
        
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
    # lv = {
    #     f: 'F',r: 'R', b: 'B'
    # }
    # convery the lists to a numpy array
    X = numpy.array(mat)
    X = X.astype(numpy.float64)
    print(X)
    for level in thresholds:
        # print(int(level))
        clustering = AgglomerativeClustering(n_clusters = None, affinity = 'precomputed',linkage = 'single', distance_threshold =int(level)).fit(X)
        df = pandas.DataFrame(data = {'ID': isos, f"Tx:{level}": clustering.labels_})
        df[f"Tx:{level}"] = df[f"Tx:{level}"] + 1
        df = df.fillna('')
        df[f"Tx:{level}"] = df[f"Tx:{level}"].apply(lambda x: f"{x}")
        # df[lv[level]] = df[lv[level]].apply(lambda x: f'{lv[level]}_{x}')
        
        uc = extract_unclustered(df = df, col = f"Tx:{level}", col1 = 'ID')
        df[f"Tx:{level}"] = numpy.where(df[f"Tx:{level}"].isin(uc), 'UC', df[f"Tx:{level}"])
        if clusters.empty:
            clusters = df
        else:
            clusters = pandas.merge(df, clusters)
        # if 'Genomic_cluster' in list(clusters.columns):
        #     clusters['Genomic_cluster'] = clusters[['Genomic_cluster', lv[level]]].apply(lambda x: '_'.join(x), axis = 1)
    # clusters['Genomic_clusters_001'] = clusters[['Cluster_001', 'Cluster_0005']].apply(lambda x: '_'.join(x), axis = 1)
    clusters = clusters.fillna('')

    clusters.to_csv('clusters.txt', sep = '\t', index = False)