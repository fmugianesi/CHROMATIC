########
# create .npy array from the .bed file of the ChIP-seq peaks 
# (example of one line of the input .bed file: chr1 3977462 3977612)
########

path = '/mnt/cnag/github/'

import numpy as np 

res = 5000   #resolution

### mm10 genome
# chrs: from 0 to total number-1
chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
# chrs: names of chromosomes
chrom_strings = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
# chr sizes in base pairs
chr_sizes = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299]

# list of all ChIPseq tracks, controls excluded
# name-of-the-factor
ChIPs = ['Cbx3', 'CTCF', 'Ezh2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'Nkx2.2', 'Nkx6.1', 'Olig2', 'Pcgf2', 'PolII', 'Ring1b', 'SMC1', 'Sox2', 'Sox3', 'Suz12', 'Zrf1' ]

for ChIP in ChIPs:
    for chrom in chroms:
        L = chr_sizes[chrom]
        nb = L/res + 1
        C = np.zeros(nb)
        with open(path + '/data/NPC/ChIPseq/ChIPseq_peaks/' + ChIP + '_peaks_ucsc.bed') as tsvfile:   # .bed file with ChIP-seq peaks (for Cbx3, the file name is 'Cbx3_peaks_ucsc.bed')
            for line in tsvfile:
                chromosome, s1, s2 = line.split()
                if chromosome == chrom_strings[chrom]:
                    s = np.arange(int(s1), int(s2)+1) 
                    for s_i in s:
                        C[s_i/res] = 1

            np.save(path + '/data/NPC/ChIPseq/ChIPseq_peaks/npy/' + ChIP + '_peaks_' + chrom_strings[chrom] + '.npy', C)

