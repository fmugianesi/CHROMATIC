import numpy as np
import scipy
from scipy import signal
from scipy import ndimage

path = '/mnt/cnag/github/'

res = 5000   #resolution

### mm10 genome
# chrs: names of chromosomes
chrom_strings = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
# chr sizes in base pairs
chr_sizes = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299]

chroms = [0]   # with this script we are processing only chr1
# if you want to use this same script for all chrs, uncomment the following line (it will take more time of course):
#chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

ChIPs = ['Cbx3_Huang17']   # with this script we are processing only CBX3
# if you want to use this same script for all the factors, uncomment the following line (not recommended for the sake of time):
# ChIPs = ['Cbx3_Huang17', 'CTCF_Bonev17', 'Ezh2_Kloet16', 'H3K27ac_Bonev17', 'H3K27me3_Bonev17', 'H3K4me3_Bonev17', 'H3K9me3_Bonev17, 'Nkx2.2_Nishi15', 'Nkx6.1_Nishi15', 'Olig2_Nishi15', 'Pcgf2_Kloet16', 'PolII_Huang17', 'Ring1b_Bonev17', 'SMC1_Cremins13', 'Sox2_Lodato13', 'Sox3_McAninch14', 'Suz12_Kloet16', 'Zrf1_Aloia14']

########################### CODE

# parameters of CHROMATIC that are used later
a = 0.5
b = 0.7937
th = 0.5

for chrom in chroms:
        chrom_reg = chrom_strings[chrom]
        L = chr_sizes[chrom]
        nb_tot = L/res + 1  # total number of bins of the chromosome

        # load Hi-C
        W = scipy.sparse.load_npz(path + 'data/NPC/HiC/npz/HiC_' + chrom_reg + '_median.npz')

        for ChIP in ChIPs:

                w = 10000000/res # window of 10Mb (in our case res=5000, thus w=2000. if you significantly increase your resolution (e.g. res=1000 bp), you are also increasing w and the computational burden)
                t = np.arange(w)
                nb = w
                step = 2000000/res # step of 2Mb with which I move the window

                # s,e: indices of the bins of start and end of the sliding window in the reference system of the whole chromosome
                S = np.arange(0, (((nb_tot-w)/step)*step+step+1), step = step) # all the starting points of the sliding widows in a chromosome
                s_last = S[-1]

                Cw = np.load(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_reg + '_ChIPseq_track_norm_01.npy') # ChIPseq track of the full chromosome

                ChIP = ChIP.split('_')[0]   # important, otherwise error because the next file names are without the name-of-the-study. Comment this line if you don't have this problem.
                Pt = np.load(path + 'data/NPC/ChIPseq/ChIPseq_peaks/npy/' + ChIP + '_peaks_' + chrom_reg + '.npy') # ChIP-seq peaks of the full chromosome for CBX3

                M_all = sparse.csr_matrix((nb_tot, nb_tot)) # collects the CHROMATIC interaction map in sparse format, for the whole chromosome, for CBX3 
                Macs_all = sparse.csr_matrix((nb_tot, nb_tot)) # collects ChIP-seq peaks matrix in sparse format, for the whole chromosome, for CBX3

                for s in S:
                        if s == s_last:
                                e = nb_tot
                                nb = e-s
                                t = np.arange(nb)
                        else:
                                e = s+w

                        H = W[s:e,s:e]
                        C = Cw[s:e]

                        M = np.zeros([nb,nb])
                        for i in t:
                                for j in t:
                                        if np.sum(H[i,])>0 and np.sum(H[j,])>0:  # only if both bins are mappable 
                                                
                                                # the following transform applied to ChIP-seq values "stretches" their distribution. If you ever dealt with images, it is like a contrast enhancement. See the paper where there is a plot of this used transform, where it is visible that lower ChIP-seq values are lowered further, and high ChIP-seq values are increased further.
                                                if C[i]<th :
                                                        C1 = -(abs(C[i]-a))**(1./3)+b   # lowered
                                                elif C[i]>=th :
                                                        C1 = 3**C[i]                    # increased
                                                if C[j]<th :
                                                        C2 = -(abs(C[j]-a))**(1./3)+b   # lowered
                                                elif C[i]>=th :
                                                        C2 = 3**C[j]                    # increased
                                               
                                                M[i,j] = C1*C2*H[i,j]   # this is the new coefficient in the CHROMATIC interaction map

                        median = signal.medfilt(M, kernel_size = 5)
                        M = median

                        # here I create a sparse matrix Peaks_m that has value Peaks_m[i,j]=1 if the bins i-1,i,i+1 OR the bins j-1,j,j+1 contain at least one ChIP-seq peak, and Peaks_m[i,j]=0 otherwise. It will be used later.
                        Peaks = Pt[s:e]   # ChIP-seq peaks in the considered window
                        Peaks_m = np.zeros([nb,nb])
                        for i in t:
                                for j in t:
                                        if Peaks[i]==1 or Peaks[j]==1:
                                                si = i-1
                                                if si<0:
                                                        si = 0
                                                ei = i+2
                                                if ei>nb:
                                                        ei = nb
                                                sj = j-1
                                                if sj<0:
                                                        sj = 0
                                                ej = j+2
                                                if ej>nb:
                                                        ej = nb
                                                Peaks_m[si:ei, sj:ej] = np.ones((ei-si, ej-sj))

                        M = (M + np.transpose(M))/2.0   # symmetrize the matrix
                        mm = np.amax(M)
                        M = M/mm    
                        where_are_NaNs = np.isnan(M)
                        M[where_are_NaNs] = 0    # put NaNs at 0

                        M_all[s:e,s:e] = sparse.lil_matrix(M)  
                        Macs_all[s:e,s:e] = sparse.lil_matrix(Peaks_m)

                scipy.sparse.save_npz(path + 'data/NPC/CHROMATIC/int_maps/' + ChIP + '_' + chrom_reg + '_CHROMATIC_int_map.npz', M_all.tocsr())
                scipy.sparse.save_npz(path + 'data/NPC/ChIPseq/ChIPseq_peaks/npz/' + ChIP + '_peaks_' + chrom_reg + '.npz', Macs_all.tocsr())

