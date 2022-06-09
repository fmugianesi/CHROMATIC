# CHROMATIC  
CHROMATIC is a computational tool that integrates Hi-C and ChIP-seq data to study chromatin three-dimensional (3D) interactions associated with any factor of interest.   
It is faster and less expensive than performing experiments that probe protein-directed genome architecture, such as HiChIP. 
Thanks to the deconvolution of the Hi-C data into factor-specific interactions, our strategy allows discerning the role of each studied factor in genome 3D structure in a cell-type-specific manner.  
Furthermore, the classification of 3D colocalization patterns of factors using CHROMATIC identifies types of functional 3D interactions, that we call ‘3D-types’. 3D-types may reflect already known interactions between different chromatin factors or may help discover new associations between molecules with specific functional roles.  
Overall, the CHROMATIC tool unifies factor occupancy and genome topology analyses, to shed light on their link with gene expression. 
 

# How to use

## ChIP-seq data pre-processing

### 01. Transform the ChIP-seq peaks .bed file into .npy array
Use the script '01_ChIpseq_peaks_fromBed_toNpy.py' contained in the subfolder './code/ChIPseq_prepro/' to generate .npy arrays from the .bed files of ChIP-seq peaks locations.  
The input is a .bed file, one for each factor you study. See the input I used in the application in NPCs in the subfolder './data/NPC/ChIPseq/ChIPseq_peaks/bed/', where the format of the name of files is '(name-of-the-factor)_peaks_ucsc.bed'.  
The ouput is one .npy array for every factor and every chromosome, saved in the subfolder './data/NPC/ChIPseq/ChIPseq_peaks/npy/'.  

### 02. From the .bedgraph files of ChIP-seq tracks obtain the .npy normalized array that later will be combined with Hi-C
Use the script '02_ChIpseq_tracks_fromBedgraph_toNormalized.py' contained in the subfolder './code/ChIPseq_prepro/'.  
The input is .bedgraph files, one for each factor you study. See the input I used in the application in NPCs in the subfolder './data/NPC/ChIPseq/ChIPseq_tracks/bedgraph/', where the format of the name of files is '(name-of-the-factor).bedgraph'.  
This scripts involves several steps, where intermediate outputs are stored in the subfolder '.data/NPC/ChIPseq/ChIPseq_tracks/', in case the run gets interrupted and you don't want to start from the beginning. You can later delete these intermediate .npy arrays. DO NOT delete the final .npy arrays, saved in './data/NPC/ChIPseq/ChIPseq_tracks/', whose name ends with '_norm_01.npy'. They are the ChIP-seq tracks normalized in the range 0-1, which will be next combined with Hi-C.  
In the subfolder './data/NPC/ChIPseq/ChIPseq_tracks/' you can find the final data I obtained in the application of CHROMATIC to 18 factors in NPCs.

### 03. Optional: Plot the normalized ChIP-seq tracks in a region of interest
Use the script '03_plot_normalizedChIPseqs.py' contained in the subfolder './code/ChIPseq_prepro/' to plot the normalized ChIP-seq tracks in a region of interest. The script is set to perform this plot on the HoxA locus of the mouse genome mm10.





