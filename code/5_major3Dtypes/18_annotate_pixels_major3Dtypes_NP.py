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

shift = 0

n_clu = 4   # number of clusters we previously identified. change this parameter according to the number of clusters you identified
counter = np.zeros(n_clu)

for chrom in chroms:
    chrom_name = chrom_strings[chrom]
    L = chr_sizes[chrom]
    nb_tot = L/res + 1
    Painted = sparse.lil_matrix((nb_tot, nb_tot))
    P_sum = sparse.load_npz(path + 'data/NPC/CHROMATIC/3Dtypes/sum_CHROMATIC_peaks_allfactors_' + chrom_name + '.npz')
    P_tri = np.triu(P_sum.toarray())
    A,B = np.where(P_tri>0)
    nab = len(A)  # number of detected pixels in the chromosome
    for i in np.arange(nab): 
        topic = np.where(int_type[shift+i,].toarray()==np.amax(int_type[shift+i,].toarray()))[1] # the 3D type associated with the pixel is the one for which the row has the maximum


# ! since we have 4 major clusters, we generate a map that can assume value 0, 1, 2, 3, 4. each value corresponds to a different classification: 0:unclassified, 1:nTFs, 2:Inactive, 3:Active, 4:PcG-Bivalent

# ! change the parameters relative to which 3D types compose every major 3D type according to your findings
# ! if for example you found 5 major 3D types, add a new "if" where you assign clu=5

        # in our application case in NPC, 3D types 0, 1, 5, 9 constituted the major 3D type of Neuronal TFs
        if topic==0 or topic==1 or topic==5 or topic==9 :
            clu = 1   # nTFs
            counter[clu-1] += 1   # counts interactions for this major type

        # in our application case in NPC, 3D types 2, 3, 12 constituted the major 3D type of Inactive
        if topic==2 or topic==3 or topic==12 :
            clu = 2   # Inactive
            counter[clu-1] += 1   # counts interactions for this major type
        
        # in our application case in NPC, 3D types 4, 6, 6=7, 8, 10, 11, 13, 15, 16 constituted the major 3D type of Active
        if topic==4 or topic==6 or topic==7 or topic==8 or topic==10 or topic==11 or topic==13 or topic==15 or topic==16:
            clu = 3   # Active
            counter[clu-1] += 1   # counts interactions for this major type
  
        # in our application case in NPC, 3D type 4 constituted the major 3D type of PcG-Bivalent
        if topic==14 :
            clu = 4   # PcG-Bivalent
            counter[clu-1] += 1   # counts interactions for this major type
        Painted[A[i],B[i]] = clu

    sparse.save_npz(path + 'data/NPC/CHROMATIC/major3Dtypes/painted_hicmap_major3Dtypes_' + chrom_name + '.npz', Painted.tocsr())
    shift += nab

np.save(path + 'data/NPC/CHROMATIC/major3Dtypes/counter_pixels_major3Dtypes_NP.npy', counter)

print 'Total number of classified pixels: ' + str(shift) 

