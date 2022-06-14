# compute and save Size of detections (in terms of number of pixels)

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

chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]   # you can process all chromosomes in this script
# if you want to process only one chr, such as chr 1, uncomment the following line:
#chroms = [0]

ChIPs = ['Cbx3']   # with this script we are processing only CBX3, but you can add more factors (it is pretty fast)

for chrom in chroms:
    chrom_reg = chrom_strings[chrom]
    L = chr_sizes[chrom]
    nb_tot = L/res + 1

    for ChIP in ChIPs:

        P = sparse.load_npz(path + 'data/NPC/CHROMATIC/detections/' + ChIP + '_' + chrom_reg + '_CHROMATIC_peaks.npz')   # load the map of CHROMATIC detections
        TP = P.toarray()
        TPt = np.triu(TP)
        label_TP, nTP = ndimage.label(TPt)
        GD_TP = np.zeros(nTP)
        i = 0
        for p in np.unique(label_TP)[1:] :  # for each label = for each detected patch
            x,y = np.where(label_TP==p)   # bins corresponding to the pixels of the detection
            gd = x.shape[0]   # number of pixels in the patch/detection
            GD_TP[i] = gd
            i += 1
        np.save(path + 'data/NPC/CHROMATIC/detections_analysis/size_CHROMATIC_peaks_' + ChIP + '_' + chrom_reg + '.npy', GD_TP)



