import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from matplotlib import colors
import copy

path = '/mnt/cnag/github/'

res = 5000   #resolution

region1 = ['chr18', 53800000, 56600000, 'Zfp608']    # region chr18:53800000-56600000 containing the gene Zfp608
regions = [region1]   # if you have more regions, do:   region2 = [...]   regions = [region1, region2]

# modify the color map according to your number of major 3D types
# colors that we associate to the 4 major 3D types:
Inactive = '#7fcdbb'
Active = '#D6311F'
Bivalent = '#8E9FCB'
nTFs = '#08818A'
mycmap_NPC = colors.ListedColormap(['white', nTFs, Inactive, Active, Bivalent])

for region in regions:
    chrom_name = region[0]
    start_bp = region[1]
    start_bin = start_bp/res
    end_bp = region[2]
    end_bin = end_bp/res
    nb_tot = end_bp-start_bp
    
    Painted = sparse.load_npz(path + 'data/NPC/CHROMATIC/major3Dtypes/painted_hicmap_major3Dtypes_' + chrom_name + '.npz')
    P_NPC = Painted[start_bin:end_bin, start_bin:end_bin].toarray()
    M_NPC = copy.deepcopy(P_NPC)
    f, ax = plt.subplots(1,1,figsize=(9,9))
    ax.imshow(M_NPC, cmap = mycmap_NPC, vmin = 0, vmax = 4.3)  # since we have 4 major clusters, M_NPC can assume value 0, 1, 2, 3, 4. we set vmax=4.3 so that each value has a different color: 0:white, 1:nTFs, 2:Inactive, 3:Active, 4:PcG-Bivalent
    f.suptitle(region[3] + ' - ' +chrom_name + ':' + str(start_bp) + '-' + str(end_bp), fontsize=14, weight='heavy', x=0.53, y=0.95)
    plt.savefig(path + 'data/NPC/CHROMATIC/major3Dtypes/painted_hicmap_Zfp608.png', transparent=True)
    plt.show()
    plt.close() 
    