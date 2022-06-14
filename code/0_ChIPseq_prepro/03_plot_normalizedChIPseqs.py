########
# optional: visualize normalized ChIP-seq tracks on the HoxA locus
########

import numpy as np
import matplotlib.pyplot as plt

path = '/mnt/cnag/github/'

res = 5000 
chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
chrom_strings = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
chr_sizes = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299]

region = ['chr6', 48050001, 52750006] # HoxA
regions = [region]   
# if you have more regions of interest you can do: 
# region2 = ['chr8', 12050001, 22750006]
# regions = [region1, region2]

ChIPs = ['Nkx2.2_Nishi15', 'Nkx6.1_Nishi15', 'Olig2_Nishi15',
         'Ezh2_Kloet16', 'Pcgf2_Kloet16', 'Ring1b_Bonev17', 'Suz12_Kloet16', 
         'H3K27ac_Bonev17', 'H3K4me3_Bonev17', 'PolII_Huang17', 'Zrf1_Aloia14',
        'CTCF_Bonev17', 'SMC1_Cremins13', 
         'Sox2_Lodato13', 'Sox3_McAninch14',
         'Cbx3_Huang17', 'H3K27me3_Bonev17', 'H3K9me3_Bonev17']
# labels appearing in the plot
ChIPs_lab = ['NKX2.2', 'NKX6.1', 'OLIG2', 
             'EZH2', 'PCGF2', 'RING1B', 'SUZ12', 
             'H3K27ac', 'H3K4me3',  'POLII', 'ZRF1',
             'CTCF', 'SMC1',
             'SOX2', 'SOX3',
             'CBX3', 'H3K27me3', 'H3K9me3']
nc = len(ChIPs)

# I color ChIP-seqs according to their functional role
color_pc = 'Blues'   # polycomb
color_ptfs = 'Oranges'   # pTFs
color_ntfs = 'YlGnBu'   # nTFs
color_arch = 'Greens'   # architectural proteins    
color_active = 'Reds'   # factors related to transcriptional activity
color_repr = 'Purples'   # factors related to transcriptional repression
mypalette = {'CBX3': color_repr, 'CTCF': color_arch, 'EZH2': color_pc, 'H3K27ac': color_active, 
             'H3K27me3': color_repr, 'H3K4me3': color_active, 'H3K9me3': color_repr, 'NKX2.2': color_ntfs, 
             'NKX6.1': color_ntfs, 'OLIG2': color_ntfs, 'PCGF2': color_pc, 'POLII': color_active, 
             'RING1B': color_pc, 'SMC1': color_arch, 'SOX2': color_ptfs, 'SOX3': color_ptfs,
              'SUZ12': color_pc, 'ZRF1': color_active }

for region in regions:

    chrom = region[0] 
    start = region[1]
    end = region[2]
    start_bin = start/res
    end_bin = end/res

    fig = plt.figure(1, figsize=(30,20))
    title_position = (0.5, 2.05)

    fig.subplots_adjust(top=2)
    fsz = 54

    l = len(ChIPs)
    i=0

    while i < l:
        C = np.load(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_norm_01.npy')
        C = C[start_bin:end_bin]
        nb = C.shape[0]
        height = int(round(0.03*nb))
        h = np.arange(0, nb, nb/20)
        C_H = np.tile(C,(height,1))
        ax = fig.add_subplot(l,1,i+1)

        plt.xticks([])
        plt.yticks([])
    
        plt.ylabel(ChIPs_lab[i], fontsize = fsz, rotation=0)
        ax.yaxis.set_label_coords(-0.1, -0.05)
        im = ax.imshow(C_H, cmap = mypalette[ChIPs_lab[i]])
    
        i += 1
        
    plt.tight_layout()
    #plt.savefig(path + 'NPC/data/ChIPseq/ChIPseq_tracks/ChIPseq_panel_NPC_HoxA.pdf', bbox_inches='tight')
    plt.show()
    plt.close()

