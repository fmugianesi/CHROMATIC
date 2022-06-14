import numpy as np
import scipy
from scipy import sparse
from scipy import signal

path = '/mnt/cnag/github/'

res = 5000   #resolution

### mm10 genome
# chrs: names of chromosomes
chrom_strings = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
# chr sizes in base pairs
chr_sizes = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299]

path_hic = '/scratch/production/fmugian/NPs_all/05_sub-matrices/'
path = '/fastscratch/sgtteam/fmugian/Chromatic_NPs/'

chroms = [0]   # with this script we are processing only chr1
# if you want to use this same script for all chrs, uncomment the following line (it will take more time of course):
#chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

for chrom in chroms:
        chrom_reg = chrom_strings[chrom]
        L = chr_sizes[chrom]
        nb_tot = L/res + 1
        H = np.loadtxt(path + 'data/NPC/HiC/abc/HiC_' + chrom_reg + '_5kb.abc')  # sparse Hi-C matrix of chr1
		Hsp = sparse.lil_matrix((nb_tot, nb_tot))
		rows = np.arange(H.shape[0])
		for r in rows:
        		i = int(H[r,0])
        		j = int(H[r,1])
        		Hsp[i,j] = H[r,2]
		W = Hsp.tocsr()
		sparse.save_npz(path + 'data/NPC/HiC/npz/HiC_' + chrom_reg + '_5kb.npz', W)  # save it in csr format (to enhance speed). If you want, you can delete it.
        H = signal.medfilt(H.toarray(), kernel_size = 5)  # this median filter reduces noise
        H = sparse.csr_matrix(H)
        sparse.save_npz(path + 'data/NPC/HiC/npz/HiC_' + chrom_reg + '_median.npz', H)   # this is the input for CHROMATIC, keep it :)
