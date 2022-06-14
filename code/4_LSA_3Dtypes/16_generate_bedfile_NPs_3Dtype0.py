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

Bed = []

topics_all = [0]  # only do 3Dtype 0
# uncomment the following line to perform the script for all 3D types
# topics_all = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

for top in topics_all:
        shift = 0

        for chrom in chroms:
                chrom_name = chrom_strings[chrom]
                L = chr_sizes[chrom]
                nb_tot = L/res + 1
                P_sum = sparse.load_npz(path + 'data/NPC/CHROMATIC/3Dtypes/sum_CHROMATIC_peaks_allfactors_' + chrom_name + '.npz')   # this matrix has the locations of the detections made for all factors, in each chromosome
                P_tri = np.triu(P_sum.toarray())
                A,B = np.where(P_tri>0)
                nab = len(A)  # number of 3D interactions in the chromosome
                for i in np.arange(nab): 
                        topic = np.where(int_type[shift+i,].toarray()==np.amax(int_type[shift+i,].toarray()))[1]  # find the 3D type of the pixel
                        if topic == top:
                                # save A[i],B[i] locations:
                                if Bed==[]:
                                        Bed = [chrom_name, str(A[i]*res), str((A[i]+1)*res)]
                                        Bed = np.vstack((Bed, [chrom_name, str(B[i]*res), str((B[i]+1)*res)] ))
                                else:
                                        control = 0
                                        for line in Bed:
                                                if line[0]==chrom_name and line[1]==str(A[i]*res) and line[2]==str((A[i]+1)*res):
                                                        control += 1
                                        if control == 0 :
                                                Bed = np.vstack((Bed, [chrom_name, str(A[i]*res), str((A[i]+1)*res)] ))

                                        control = 0
                                        for line in Bed:
                                                if line[0]==chrom_name and line[1]==str(B[i]*res) and line[2]==str((B[i]+1)*res):
                                                        control += 1
                                        if control == 0 :
                                                Bed = np.vstack((Bed, [chrom_name, str(B[i]*res), str((B[i]+1)*res)] ))
                            
                shift += nab

        with open(path + 'data/NPC/CHROMATIC/3Dtypes/3Dtype_' + str(top) + '_allchrs.bed', 'w') as txt:
                for line in Bed:
                        txt.write('\t'.join(line) + '\n')

                        