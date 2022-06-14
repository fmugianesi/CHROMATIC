import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

main_clusters = 4  # number of clusters we previously identified. change this parameter according to the number of clusters you identified
counters = np.load(path + 'data/NPC/CHROMATIC/major3Dtypes/counter_pixels_major3Dtypes_NP.npy')  
total = np.sum(counters)

labels = ['nTFs', 'Inactive', 'Active', 'PcG-Bivalent']
colors_clu = ['#08818A', '#80CDBB', '#D6311F', '#8E9FCB']

fig, ax = plt.subplots(1,1, figsize=(10,10))
fsz = 20
df = pd.DataFrame(counters, index=labels, columns=[''])
df.plot(kind='pie', subplots=True, ax=ax, colors=colors_clu, fontsize=fsz, legend = False, autopct='%1.1f%%')

plt.savefig(path + 'data/NPC/CHROMATIC/major3Dtypes/piechart_major3Dtypes_NPC.pdf', bbox_inches='tight')
plt.show()
plt.close()
