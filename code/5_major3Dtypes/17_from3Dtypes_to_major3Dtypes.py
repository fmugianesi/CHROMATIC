import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from __future__ import division

path = '/mnt/cnag/github/'

res = 5000   #resolution

### mm10 genome
chroms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
chrom_strings = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
chr_sizes = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299]

ChIPs = ['Cbx3', 'CTCF', 'Ezh2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'Nkx2.2', 'Nkx6.1', 'Olig2', 'Pcgf2', 'PolII', 'Ring1b', 'SMC1', 'Sox2', 'Sox3', 'Suz12', 'Zrf1']
ChIPs_lab = ['CBX3', 'CTCF', 'EZH2', 'H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'NKX2.2', 'NKX6.1', 'OLIG2', 'PCGF2', 'POLII', 'RING1B', 'SMC1', 'SOX2', 'SOX3', 'SUZ12', 'ZRF1']
nc = len(ChIPs)
n_components = nc-1

labels = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16']
features = ['Active enhancers','Active promoters','Poised enhancers','Bivalent promoters','Superenhancers','Constitutive LADs']
n_features = len(features)

topic_counter = np.load(path + 'data/NPC/CHROMATIC/3Dtypes/3Dtype_counter_allchrs.npy')  # no. of pixels per 3D type 
int_tot = np.sum(topic_counter)  # total no. of pixels assigned to 3D types

AE = [4778, 5551, 1163, 996, 2177, 2651, 2454, 1349, 1907, 1059, 575, 1921, 270, 664, 506, 405, 305]    # no. of lines in .bed files of the overlap of 3D types with AE
AP = [2831, 3560, 1725, 752, 4382, 1282, 2068, 1616, 1909, 670, 2374, 1621, 132, 649, 402, 404, 319]    # same with AP 
PE = [577, 175, 205, 174, 240, 555, 297, 312, 102, 467, 208, 383, 268, 579, 1263, 31, 215]  # same with PE
BP = [1, 1, 0, 3, 6, 1, 3, 1, 3, 0, 0, 0, 0, 2, 4, 2, 0]    # same with BP
SE = [524, 732, 93, 78, 615, 227, 450, 249, 405, 115, 215, 348, 4, 95, 35, 70, 53]      # same with SE
cL = [25430, 15270, 33006, 39557, 3495, 16417, 6929, 4401, 5220, 9717, 6662, 3434, 9769, 1357, 6392, 459, 501]    # same with CL



##########
# Compute the Log(OddsRatio) of the overlaps    (innput for PCA)
##########

#### Active Enhancers

# only consider subset genome (interactions that I pick)

F_T = np.zeros(n_components)
F_notT = np.zeros(n_components)
notF_T = np.zeros(n_components)
notF_notT = np.zeros(n_components)
OR = np.zeros(n_components)
OR_log_AE = np.zeros(n_components)

for i in np.arange(n_components):
    F_T[i] = AE[i]   # AE in T0
    F_notT[i] = np.sum(AE)-AE[i]   # AE that are not in T0
    notF_T[i] = topic_counter[i]-AE[i]   # T0 that are not AE
    notF_notT[i] = np.sum(topic_counter)-topic_counter[i] - np.sum(AE)-AE[i]   # not in T0, not AE

    OR[i] = (F_T[i]/F_notT[i])/(notF_T[i]/notF_notT[i])
    OR_log_AE[i] = np.log10(OR[i])



#### Active Promoters

# only consider subset genome (interactions that I pick)

F_T = np.zeros(n_components)
F_notT = np.zeros(n_components)
notF_T = np.zeros(n_components)
notF_notT = np.zeros(n_components)
OR = np.zeros(n_components)
OR_log_AP = np.zeros(n_components)

for i in np.arange(n_components):
    F_T[i] = AP[i]   # AP in T0
    F_notT[i] = np.sum(AP)-AP[i]   # AP that are not in T0
    notF_T[i] = topic_counter[i]-AP[i]   # T0 that are not AP
    notF_notT[i] = np.sum(topic_counter)-topic_counter[i] - np.sum(AP)-AP[i]   # not in T0, not AP

    OR[i] = (F_T[i]/F_notT[i])/(notF_T[i]/notF_notT[i])
    OR_log_AP[i] = np.log10(OR[i])



#### Poised Enhancers

# only consider subset genome (interactions that I pick)

F_T = np.zeros(n_components)
F_notT = np.zeros(n_components)
notF_T = np.zeros(n_components)
notF_notT = np.zeros(n_components)
OR = np.zeros(n_components)
OR_log_PE = np.zeros(n_components)

for i in np.arange(n_components):
    F_T[i] = PE[i]   # PE in T0
    F_notT[i] = np.sum(PE)-PE[i]   # PE that are not in T0
    notF_T[i] = topic_counter[i]-PE[i]   # T0 that are not PE
    notF_notT[i] = np.sum(topic_counter)-topic_counter[i] - np.sum(PE)-PE[i]   # not in T0, not PE

    OR[i] = (F_T[i]/F_notT[i])/(notF_T[i]/notF_notT[i])
    OR_log_PE[i] = np.log10(OR[i])



#### Bivalent Promoters

F_T = np.zeros(n_components)
F_notT = np.zeros(n_components)
notF_T = np.zeros(n_components)
notF_notT = np.zeros(n_components)
OR = np.zeros(n_components)
OR_log_BP = np.zeros(n_components)

for i in np.arange(n_components):
    F_T[i] = BP[i]   # BP in T0
    F_notT[i] = np.sum(BP)-BP[i]   # BP that are not in T0
    notF_T[i] = topic_counter[i]-BP[i]   # T0 that are not BP
    notF_notT[i] = np.sum(topic_counter)-topic_counter[i] - np.sum(BP)-BP[i]   # not in T0, not BP

    OR[i] = (F_T[i]/F_notT[i])/(notF_T[i]/notF_notT[i])
    OR_log_BP[i] = np.log10(OR[i])
    OR_log_BP[np.isinf(OR_log_BP)]=0


#### Super-enhancers

F_T = np.zeros(n_components)
F_notT = np.zeros(n_components)
notF_T = np.zeros(n_components)
notF_notT = np.zeros(n_components)
OR = np.zeros(n_components)
OR_log_SE = np.zeros(n_components)

for i in np.arange(n_components):
    F_T[i] = SE[i]   # SE in T0
    F_notT[i] = np.sum(SE)-SE[i]   # SE that are not in T0
    notF_T[i] = topic_counter[i]-SE[i]   # T0 that are not SE
    notF_notT[i] = np.sum(topic_counter)-topic_counter[i] - np.sum(SE)-SE[i]   # not in T0, not SE

    OR[i] = (F_T[i]/F_notT[i])/(notF_T[i]/notF_notT[i])
    OR_log_SE[i] = np.log10(OR[i])



#### Constituve LADs

F_T = np.zeros(n_components)
F_notT = np.zeros(n_components)
notF_T = np.zeros(n_components)
notF_notT = np.zeros(n_components)
OR = np.zeros(n_components)
OR_log_cL = np.zeros(n_components)

for i in np.arange(n_components):
    F_T[i] = cL[i]   # cL in T0
    F_notT[i] = np.sum(cL)-cL[i]   # cL that are not in T0
    notF_T[i] = topic_counter[i]-cL[i]   # T0 that are not cL
    notF_notT[i] = np.sum(topic_counter)-topic_counter[i] - np.sum(cL)-cL[i]   # not in T0, not cL

    OR[i] = (F_T[i]/F_notT[i])/(notF_T[i]/notF_notT[i])
    OR_log_cL[i] = np.log10(OR[i])


data = np.column_stack((OR_log_AE, OR_log_AP, OR_log_SE, OR_log_PE, OR_log_BP, OR_log_cL))
lab_type = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16']
lab_feat = ['AE', 'AP', 'SE', 'PE', 'BP', 'cL']
df = pd.DataFrame(data, columns=lab_feat, index=lab_type)



##########
# PCA
##########

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

scaler = StandardScaler()
segmentation_std = scaler.fit_transform(df)

pca = PCA()
pca.fit(segmentation_std)

pca.explained_variance_ratio_
pca.explained_variance_ratio_.cumsum()  # choose number of PCs that explain at least the 80% of the total variance

pca = PCA(n_components=3)  # in our case, in NPCs, the first 3 PCs explained the 90.88% of the total variance
pca.fit(segmentation_std)
pca.transform(segmentation_std)
scores_pca = pca.transform(segmentation_std)



##########
# K-means clustering
##########

wcss = []
for i in range(1,18):
    kmeans_pca = KMeans(n_clusters=i, init='k-means++', random_state=42)
    kmeans_pca.fit(scores_pca)
    wcss.append(kmeans_pca.inertia_)


# choose the number of clusters for which you have the kink in the following plot:

fig, ax = plt.subplots(1,1, figsize=(10,10))
plt.plot(range(1,18), wcss, marker='o', linestyle='-', color = '#66c2a4')
fsz = 30
plt.xlabel('No. of clusters', fontsize = fsz)
plt.xticks(range(1,18), fontsize = fsz-10)
plt.yticks(np.arange(0, 101, step=20), fontsize = fsz-10)
plt.ylabel('WCSS', fontsize = fsz)
for _,s in ax.spines.items():
    s.set_linewidth(1.3)
ax.tick_params(axis='both', length = 17, width = 1.6, pad = 18) 
plt.savefig(path + 'data/NPC/CHROMATIC/major3Dtypes/kink_NPC.pdf', bbox_inches='tight')
plt.show()
plt.close()

# we got 4 optimal clusters (4 major 3D-types)
kmeans_pca = KMeans(n_clusters=4, init='k-means++', random_state=42)
kmeans_pca.fit(scores_pca)

# see how 3D-types got clustered into 4 major 3D-types
df_segm_pca_kmeans = pd.concat([df.reset_index(drop=True), pd.DataFrame(scores_pca)], axis=1)
df_segm_pca_kmeans.columns.values[-3:] = ['PC1', 'PC2', 'PC3']
df_segm_pca_kmeans['Segment K-means PCA'] = kmeans_pca.labels_
df_segm_pca_kmeans



##########
# Plot the clusters on PC1 and PC2
##########

df_segm_pca_kmeans['Major types'] = df_segm_pca_kmeans['Segment K-means PCA'].map({0:'PcG-Bivalent', 
                                                                              1:'Inactive',
                                                                              2:'Active',
                                                                              3:'nTFs'})
fsz = 30
fig, ax = plt.subplots(1,1, figsize=(13,13))
x_axis = df_segm_pca_kmeans['PC2']
y_axis = df_segm_pca_kmeans['PC1']
palette = ['#08818A', '#7fcdbb', '#D6311F', '#8E9FCB']
sns.scatterplot(x_axis, y_axis, hue = df_segm_pca_kmeans['Major types'], palette = palette, s = 160, ax = ax)
for i in np.arange(n_components):
    if i == 11:
        plt.annotate(str(i), (x_axis[i]+0.08, y_axis[i]-0.15), fontsize = 30, rotation = 0)
    else:
        plt.annotate(str(i), (x_axis[i]+0.08, y_axis[i]), fontsize = 30, rotation = 0)
ax.set_xticks(np.arange(-2, 4, step=1))
ax.set_xticklabels(np.arange(-2, 4, step=1), fontsize = fsz-5)
ax.set_yticks(np.arange(-2, 5.1, step=1))
ax.set_yticklabels(np.arange(-2, 5.1, step=1).astype(int), fontsize = fsz-5)
ax.set_xlim([-2.5, 3.5])
ax.set_ylim([-2.3, 4.5])
for _,s in ax.spines.items():
    s.set_linewidth(1.6)
    ax.tick_params(axis='both', length = 17, width = 1.6, pad = 18) 
plt.xlabel('PC2', fontsize = fsz)
plt.ylabel('PC1', fontsize = fsz)
plt.savefig(path + 'data/NPC/CHROMATIC/major3Dtypes/PCA_NPC.pdf', bbox_inches='tight')
plt.show()
plt.close()



