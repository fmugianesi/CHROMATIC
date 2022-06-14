import numpy as np
from scipy import sparse

path = '/mnt/cnag/github/'

res = 5000   #resolution

### mm10 genome
chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
chrom_strings = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
chr_sizes = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299]

ChIPs = ['Cbx3', 'CTCF', 'Ezh2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'Nkx2.2', 'Nkx6.1', 'Olig2', 'Pcgf2', 'PolII', 'Ring1b', 'SMC1', 'Sox2', 'Sox3', 'Suz12', 'Zrf1']

nC = len(ChIPs)

for chrom in chroms:
        chrom_name = chrom_strings[chrom]
        L = chr_sizes[chrom]
        nb_tot = L/res + 1
        P_sum =  sparse.csr_matrix((nb_tot, nb_tot)) # collects peaks of all proteins
        for c in np.arange(nC):
                ChIP = ChIPs[c]
                P = sparse.load_npz(path + 'data/NPC/CHROMATIC/detections/' + ChIP + '_' + chrom_reg + '_CHROMATIC_peaks.npz')   # load the map of CHROMATIC detections
                P_sum += P   # sum it in P_sum, so that P_sum[i,j]>=1 if in [i,j] there is a CHROMATIC detection for at least one factor
        sparse.save_npz(path + 'data/NPC/CHROMATIC/3Dtypes/sum_CHROMATIC_peaks_allfactors_' + chrom_name + '.npz', P_sum)
        