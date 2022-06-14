import numpy as np
import pandas as pd
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

chroms = [0]   # with this script we are processing only chr1
# if you want to use this same script for all chrs, uncomment the following line (it will take more time of course):
#chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

for chrom in chroms:
    chrom_name = chrom_strings[chrom]
    L = chr_sizes[chrom]
    nb_tot = L/res + 1
    P_sum = sparse.load_npz(path + 'data/NPC/CHROMATIC/3Dtypes/sum_CHROMATIC_peaks_allfactors_' + chrom_name + '.npz')
    P_tri = np.triu(P_sum.toarray())
    A,B = np.where(P_tri>0)
    nab = len(A)
    Obs_docu = {}
    P_checkfirst = sparse.lil_matrix((nb_tot, nb_tot))
    for c in np.arange(nC):
        ChIP = ChIPs[c]
        ChIP_lab = ChIPs_lab[c]
        P = sparse.load_npz(path + 'data/NPC/CHROMATIC/detections/' + ChIP + '_' + chrom_reg + '_CHROMATIC_peaks.npz')
        for i in np.arange(nab):
            if P[A[i],B[i]]>0 :
                if P_checkfirst[A[i],B[i]]==0 :
                    Obs_docu[A[i]*nb_tot + B[i]] = 'Apple' + ' ' + ChIP_lab    
                    P_checkfirst[A[i],B[i]] = 1
                elif P_checkfirst[A[i],B[i]]>0 :
                    Obs_docu[A[i]*nb_tot + B[i]] = Obs_docu[A[i]*nb_tot + B[i]] + ' ' + ChIP_lab

    df = pd.Series(Obs_docu)
    df.sort_index(inplace=True)
    np.save(path + 'data/NPC/CHROMATIC/3Dtypes/3Dint-factor_lsa_' + chrom_name + '.npy', df)


