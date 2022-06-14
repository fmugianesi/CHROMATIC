import numpy as np
from scipy import ndimage
from scipy import sparse

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

ChIPs = ['Cbx3']   # with this script we are processing only CBX3
# if you want to process all factors, uncomment the following:
#ChIPs = ['Cbx3', 'CTCF', 'Ezh2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'Nkx2.2', 'Nkx6.1', 'Olig2', 'Pcgf2', 'PolII', 'Ring1b', 'SMC1', 'Sox2', 'Sox3', 'Suz12', 'Zrf1']


thrs = [0.2]   # threshold for peak calling. If you want to try more values, you can add them here: e.g., thrs = [0.2, 0.3]

# structuring elements that are used in the morphological image processing, in order to detect loops and hubs associated with CBX3
e1 = ndimage.generate_binary_structure(2, 1)
e1 = e1.astype(np.int)
e4 = np.array(np.ones([4,4]))
e5 = np.array(np.ones([5,5]))

for chrom in chroms:
        chrom_reg = chrom_strings[chrom]
        L = chr_sizes[chrom]
        nb_tot = L/res + 1  # total number of bins of the chromosome

        for ChIP in ChIPs:
                M_all = sparse.load_npz(path + 'data/NPC/CHROMATIC/int_maps/' + ChIP + '_' + chrom_reg + '_CHROMATIC_int_map.npz')   # CHROMATIC interaction map
                Macs_all = sparse.load_npz(path + 'data/NPC/ChIPseq/ChIPseq_peaks/npz/' + ChIP + '_peaks_' + chrom_reg + '.npz')   # matrix of the ChIP-seq peaks

                for thr in thrs:
                        w = 10000000/res # window of 10Mb (in our case res=5000, thus w=2000. if you significantly increase your resolution (e.g. res=1000 bp), you are also increasing w and the computational burden)
                        t = np.arange(w)
                        nb = w
                        step = 2000000/res # step of 2Mb with which I move the window


                        # s,e: indices of the bins of start and end of the sliding window in the reference system of the whole chromosome
                        S = np.arange(0, (((nb_tot-w)/step)*step+step+1), step = step) # all the starting points of the sliding widows in a chromosome
                        s_last = S[-1]

                        P_all = sparse.csr_matrix((nb_tot, nb_tot)) # collects the detected CHROMATIC peaks (loops/hubs) in sparse format, for the whole chromosome

                        for s in S:
                                if s == s_last:
                                        e = nb_tot
                                        nb = e-s
                                        t = np.arange(nb)
                                else:
                                        e = s+w

                                M = M_all[s:e,s:e].toarray()
                                Peaks_m = Macs_all[s:e,s:e].toarray()

                                M_th = M > thr
                                M_th = np.multiply(M_th, Peaks_m)   # this operation ensures that for each detected interaction (loop/hub) there is at least 1 ChIP-seq peak at one of the two interacting loci
                                P_fin = np.zeros([nb,nb])
                                # morphological image processing:
                                P = ndimage.binary_opening(M_th, structure = e4).astype(M.dtype)
                                P = ndimage.binary_closing(P, structure = e1).astype(M.dtype)
                                P = ndimage.binary_dilation(P, structure = e5).astype(M.dtype)
                                P = ndimage.binary_closing(P, structure = e5).astype(M.dtype)
                                P_fin[np.where(P>0)] = 1

                                P_all[s:e,s:e] = sparse.csr_matrix(P_fin)

                        sparse.save_npz(path + 'data/NPC/CHROMATIC/detections/' + ChIP + '_' + chrom_reg + '_CHROMATIC_peaks.npz', P_all)



