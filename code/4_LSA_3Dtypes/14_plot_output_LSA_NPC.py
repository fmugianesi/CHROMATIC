import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

path = '/mnt/cnag/github/'

ChIPs = ['Cbx3', 'Ezh2', 'Pcgf2', 'Ring1b', 'Suz12', 
         'H3K27me3', 'H3K9me3', 'Nkx2.2', 'Nkx6.1', 'Olig2', 
         'CTCF', 'SMC1', 
         'H3K27ac', 'H3K4me3', 'PolII', 'Zrf1', 
         'Sox2', 'Sox3']
ChIPs_lab = ['CBX3', 'EZH2', 'PCGF2', 'RING1B', 'SUZ12', 'H3K27me3', 'H3K9me3', 'NKX2.2', 'NKX6.1', 'OLIG2', 
             'CTCF', 'SMC1', 'H3K27ac', 'H3K4me3', 'POLII', 'ZRF1', 'SOX2', 'SOX3']
nc = len(ChIPs)

# I color factors according to their functional role:
color_pc = '#8da0cb'  # polycomb
color_ptfs = '#fc8d62'  # pluripotency TFs
color_ntfs = '#02818a'  # neuronal TFs  
color_arch = '#66c2a5'  # architectural proteins
color_active = '#d7301f'    # factors associated with transcriptional activity
color_repr = '#984ea3'  # factors associated with transcriptional silencing
# array with the color asssociated to each factor, following the order in which factor appears in ChIPs/ChIPs_lab:
types_colors = [color_repr, color_pc, color_pc, color_pc, color_pc,
                color_repr, color_repr, color_ntfs, color_ntfs, color_ntfs, 
                color_arch, color_arch, 
                color_active, color_active, color_active, color_active,
                color_ptfs, color_ptfs ]
                     
data = np.load(path + 'data/NPC/CHROMATIC/3Dtypes/3Dtype-factor_allchrs.npy')
data = data[:,1:]
df = pd.DataFrame(data, columns=ChIPs_lab)

sns.set(font_scale=2.1)
fsz = 34
vm = 2  # adjust this parameter if needed for visualization purposes
C = sns.clustermap(df, metric='canberra', cmap="vlag", 
                     col_colors=types_colors, 
                     z_score=0,
                     row_cluster=False, 
                     vmin=-vm,vmax=vm, 
                     method = 'complete',
                     figsize = (12,12), annot_kws={"size": fsz})
                      

plt.savefig(path + 'data/NPC/CHROMATIC/3Dtypes/LSA_clusterheatmap_NPC.pdf', bbox_inches='tight')
plt.show()
plt.close()    


