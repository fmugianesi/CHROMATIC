import numpy as np
from scipy import sparse

path = '/mnt/cnag/github/'

res = 5000   #resolution

### mm10 genome
chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
chrom_strings = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
chr_sizes = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299]

ChIPs = ['Cbx3', 'CTCF', 'Ezh2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'Nkx2.2', 'Nkx6.1', 'Olig2', 'Pcgf2', 'PolII', 'Ring1b', 'SMC1', 'Sox2', 'Sox3', 'Suz12', 'Zrf1']
ChIPs_lab = ['CBX3', 'CTCF', 'EZH2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'NKX2.2', 'NKX6.1', 'OLIG2', 'PCGF2', 'POLII', 'RING1B', 'SMC1', 'SOX2', 'SOX3', 'SUZ12', 'ZRF1']

nC = len(ChIPs)
n_components = nC-1

int_type = sparse.load_npz(path + 'data/NPC/CHROMATIC/3Dtypes/3Dint-3Dtype_allchrs.npz')
ni = int_type.shape[0]  # total number of 3D interactions genome-wide
print 'Total pixels genome-wide: ' + str(ni) 

topic_counter = np.zeros(n_components)   # number of 3D interactions for each 3D type
for j in np.arange(ni):
    topic = np.where(int_type[j,].toarray()==np.amax(int_type[j,].toarray()))[1]  # the pixel belongs to the 3D type for which it has the maximum value
    topic_counter[topic] += 1
np.save(path + 'data/NPC/CHROMATIC/3Dtypes/3Dtype_counter_allchrs.npy', topic_counter)


