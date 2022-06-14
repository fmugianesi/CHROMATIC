########
# create .npy array (normalized in 0-1 range of values, which later will be combined with Hi-C) from the .bedgraph files of the ChIP-seq tracks
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

# list of all ChIPseq tracks, controls included, used in our application in neural progenitor cells (NPCs)
# the format of the name is (name-of-the-factor)_(name-of-the-study)
# ChIP-seq tracks from the same study will be normalized with the same control
# for ex: CTCF_Bonev17 and H3K27ac_Bonev17 will be both normalized with Input_Bonev17 (see below)
ChIPs_all = ['Cbx3_Huang17', 'CTCF_Bonev17', 'Ezh2_Kloet16', 'H3K27ac_Bonev17', 'H3K27me3_Bonev17', 
             'H3K4me3_Bonev17', 'H3K9me3_Bonev17', 'IgG_Aloia14', 'Input_Bonev17',
            'Input_Cbx3_Huang17', 'Input_Cremins13', 'Input_Kloet16', 'Input_McAninch14', 'Input_Nkx2.2_Nishi15',
            'Input_Nkx6.1_Nishi15', 'Input_Olig2_Nishi15', 'Input_PolII_Huang17', 'Nkx2.2_Nishi15', 
             'Nkx6.1_Nishi15', 'Olig2_Nishi15', 'Pcgf2_Kloet16', 'PolII_Huang17', 'Ring1b_Bonev17',
             'SMC1_Cremins13', 'Sox2_Lodato13', 'Sox3_McAninch14', 'Suz12_Kloet16', 'WCE_Lodato13', 
             'Zrf1_Aloia14' ]

# list of all ChIPseq tracks, controls excluded
# name-of-the-factor
ChIPs = ['Cbx3', 'CTCF', 'Ezh2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'Nkx2.2', 'Nkx6.1', 'Olig2', 'Pcgf2', 'PolII', 'Ring1b', 'SMC1', 'Sox2', 'Sox3', 'Suz12', 'Zrf1']

res = 5000



##### from .bedgraph to .npy (raw)

for ChIP in ChIPs_all:
    for chrom in chroms:
        L = chr_sizes[chrom]
        nb = L/res + 1
        C = np.zeros(nb)
        # IT MUST BE mm10!!!!!

        with open(path + 'data/NPC/ChIPseq/ChIPseq_tracks/bedgraph/' + ChIP + '.bedgraph') as tsvfile:  # the .bedgraph file must be (name-of-the-factor).bedgraph
            for line in tsvfile:
                chromosome, s1, s2, count = line.split()
                if chromosome == chrom_strings[chrom]:
                    s = int(np.mean([int(s1),int(s2)]))
                    C[s/res] = C[s/res] + float(count)
                    
        np.save(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_raw.npy', C)  # raw ChIP-seq tracks (when you finish, you can delete them)



##### remove hot spots from raw data (hot spots are regions where ChIP-seq reads accumulate anomalously)

for ChIP in ChIPs_all:
    for chrom in chroms:
        hot = np.load(path + 'data/NPC/ChIPseq/hot_spots/hot_spots_' + chrom_strings[chrom] + '.npy')
        C = np.load(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_raw.npy')
        L = chr_sizes[chrom]
        nb = L/res + 1
        t = np.arange(nb)
        for i in t:
            if hot[i]==1:
                C[i]=0
        np.save(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_hotspotsremoved.npy', C)  # raw ChIP-seq tracks with hot spots removed (when you finish, you can delete them)



##### reduce the influence of outliers before the linear transformation, by holding out values that are more than 5 standard deviations higher than the avarage and setting these values to the relative-maximum (it will be =1 after IgG/WCE/control and linear 0-1 normalizations)

from scipy import stats

for ChIP in ChIPs_all:
    N_tot = []
    for chrom in chroms:
        N = np.load(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_hotspotsremoved.npy')
        N_tot = np.append(N_tot, N)

    Nsort =  sorted(N_tot)
    fit = stats.norm.pdf(Nsort, np.mean(Nsort), np.std(Nsort))
    th_max = np.mean(Nsort) + 5*np.std(Nsort) # 3x10^(-7) is the probabilty of values higher than 5 std dvs, 
        # if they are normally distributed (99.99994% is inside +-5 std dvs)   (from particle physics)
  


##### reduce the influence of outliers by holding out the top 0.1% of values before the linear transformation and setting these values to the relative-maximum (it will be =1 after IgG and linear 0-1 normalizations)
    for chrom in chroms:
        L = chr_sizes[chrom]
        nb = L/res + 1
        t = np.arange(nb)
        N = np.load(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_hotspotsremoved.npy')
        N_noout = N
        for i in t:
            if N[i] > th_max:
                N_noout[i] = th_max      
        np.save(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_nooutliers.npy', N_noout)  # raw ChIP-seq tracks with no hotspots and no outliers (when you finish, you can delete them)



##### normalize all ChIP-seq tracks by their corresponding control (IgG/WCE/GFP)

# these arrays correspond to ChIp-seq tracks that are from the same study, so have the same control
# place the control at the end of each array, as reported here:
Bonev17 = ['CTCF_Bonev17', 'H3K27ac_Bonev17', 'H3K27me3_Bonev17', 'H3K4me3_Bonev17', 'H3K9me3_Bonev17',
           'Ring1b_Bonev17', 'Input_Bonev17']
Huang17_Cbx3 = ['Cbx3_Huang17', 'Input_Cbx3_Huang17']
Huang17_PolII = ['PolII_Huang17', 'Input_PolII_Huang17']
Kloet16 = ['Ezh2_Kloet16', 'Pcgf2_Kloet16', 'Suz12_Kloet16', 'Input_Kloet16']
Aloia14 = ['Zrf1_Aloia14', 'IgG_Aloia14']
Cremins13 = ['SMC1_Cremins13', 'Input_Cremins13']
McAninch14 = ['Sox3_McAninch14', 'Input_McAninch14']
Nishi15_Nkx22 = ['Nkx2.2_Nishi15', 'Input_Nkx2.2_Nishi15']
Nishi15_Nkx61 = ['Nkx6.1_Nishi15', 'Input_Nkx6.1_Nishi15']
Nishi15_Olig2 = ['Olig2_Nishi15', 'Input_Olig2_Nishi15']
Lodato13 = ['Sox2_Lodato13', 'WCE_Lodato13']
# put all together
studies = [Bonev17, Huang17_Cbx3, Huang17_PolII, Kloet16, Aloia14, Cremins13, McAninch14, Nishi15_Nkx22, 
           Nishi15_Nkx61, Nishi15_Olig2, Lodato13]
    
for study in studies:
    s_t = len(study)
    for chrom in chroms:
        # load the control
        I = np.load(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + study[s_t-1] + '_' + chrom_strings[chrom] + '_ChIPseq_track_nooutliers.npy')

        # load the tracks to normalize
        for s in np.arange(s_t-1):
            L = chr_sizes[chrom]
            nb = L/res + 1
            t = np.arange(nb)
            N = np.zeros(nb)
            R = np.load(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + study[s] + '_' + chrom_strings[chrom] + '_ChIPseq_track_nooutliers.npy')
            for i in t:
                if I[i]!=0:
                    N[i] = float(R[i])/I[i]   # normalization
            np.save(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + study[s] + '_' + chrom_strings[chrom] + '_ChIPseq_track_norm.npy', N)  # normalized ChIP-seq tracks (when you finish, you can delete them)
            


##### linear transformation of the normalized tracks in the range 0-1 (in order to later be able to compare different outputs that we obtain for different ChIP-seq tracks)

for ChIP in ChIPs: 
    N_tot = []
    for chrom in chroms:
        N = np.load(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_norm.npy')
        N_tot = np.append(N_tot, N)
    
    Nsort =  sorted(N_tot)
    fit = stats.norm.pdf(Nsort, np.mean(Nsort), np.std(Nsort))
    th_max = np.mean(Nsort) + 10*np.std(Nsort) 
    mm = np.amax(N_tot)
    for chrom in chroms:
        L = chr_sizes[chrom]
        nb = L/res + 1
        t = np.arange(nb)
        N = np.load(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_norm.npy')
        N01 = np.zeros(nb)
        for i in t:
            if N[i]>th_max:
                N[i] = th_max
            N01[i] = float(N[i])/th_max
        np.save(path + 'data/NPC/ChIPseq/ChIPseq_tracks/' + ChIP + '_' + chrom_strings[chrom] + '_ChIPseq_track_norm_01.npy', N01)   # normalized ChIP-seq tracks in the range 0-1. This data will be combined with Hi-C, keep it!


