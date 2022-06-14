import time
start_time = time.time()

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import TruncatedSVD

path = '/mnt/cnag/github/'

res = 5000   #resolution

### mm10 genome
chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
chrom_strings = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
chr_sizes = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299]

ChIPs = ['Cbx3', 'CTCF', 'Ezh2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'Nkx2.2', 'Nkx6.1', 'Olig2', 'Pcgf2', 'PolII', 'Ring1b', 'SMC1', 'Sox2', 'Sox3', 'Suz12', 'Zrf1']
ChIPs_lab = ['CBX3', 'CTCF', 'EZH2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'NKX2.2', 'NKX6.1', 'OLIG2', 'PCGF2', 'POLII', 'RING1B', 'SMC1', 'SOX2', 'SOX3', 'SUZ12', 'ZRF1']

nC = len(ChIPs)
n_components = nC-1   # the maximum nu

counter = 0
for chrom in chroms:
    chrom_name = chrom_strings[chrom]
    L = chr_sizes[chrom]
    nb_tot = L/res + 1
    df = np.load(path + 'data/NPC/CHROMATIC/3Dtypes/3Dint-factor_lsa_' + chrom_name + '.npy', allow_pickle=True)
    df = pd.DataFrame(df)
    ldf = len(df)
    if counter == 0:
        df_tot = df
        counter+=1
    else:
        df = df.set_index(np.arange(len(df_tot), len(df_tot)+ldf, step=1))
        df_tot = pd.concat([df_tot, df])

df_tot = df_tot.squeeze()
# document-term matrix
vectorizer = TfidfVectorizer(stop_words='english',
max_df = 1.0,
smooth_idf=True)
X = vectorizer.fit_transform(df_tot)

# topic modeling
svd_model = TruncatedSVD(n_components=n_components, algorithm='randomized', n_iter=100, random_state=122)
lsa = svd_model.fit_transform(X)

sparse.save_npz(path + 'data/NPC/CHROMATIC/3Dtypes/3Dint-3Dtype_allchrs.npz', sparse.csr_matrix(lsa))
np.save(path + 'data/NPC/CHROMATIC/3Dtypes/3Dtype-factor_allchrs.npy', svd_model.components_)

print("--- %s seconds ---" % (time.time() - start_time))

