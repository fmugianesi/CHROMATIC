# CHROMATIC  
CHROMATIC is a computational tool that integrates Hi-C and ChIP-seq data to study chromatin three-dimensional (3D) interactions associated with any factor of interest.   
It is faster and less expensive than performing experiments that probe protein-directed genome architecture, such as HiChIP. 
Thanks to the deconvolution of the Hi-C data into factor-specific interactions, our strategy allows discerning the role of each studied factor in genome 3D structure in a cell-type-specific manner.  
Furthermore, the classification of 3D colocalization patterns of factors using CHROMATIC identifies types of functional 3D interactions, that we call ‘3D-types’. 3D-types may reflect already known interactions between different chromatin factors or may help discover new associations between molecules with specific functional roles.  
Overall, the CHROMATIC tool unifies factor occupancy and genome topology analyses, to shed light on their link with gene expression. 
 

# How to use

This tutorial explains the code, while referring to the data relative to the application of CHROMATIC in NPCs.  
You will find the ChIP-seq data for all the factors we studied at present.  
As regards Hi-C maps and downstream CHROMATIC analysis, you will find only the data for chromosome 1 and the factor CBX3.  
However, we uploaded also the data that define the 3D-types we identified genome-wide, as explained below.

## ChIP-seq data pre-processing

### 1. Transform the ChIP-seq peaks .bed file into .npy array
Use the script '01_ChIpseq_peaks_fromBed_toNpy.py' contained in the subfolder './code/ChIPseq_prepro/' to generate .npy arrays from the .bed files of ChIP-seq peaks locations.  
The input is a .bed file, one for each factor you study. See the input I used in the application in NPCs in the subfolder './data/NPC/ChIPseq/ChIPseq_peaks/bed/', where the format of the name of files is '(name-of-the-factor)_peaks_ucsc.bed'.  
The ouput is one .npy array for every factor and every chromosome, saved in the subfolder './data/NPC/ChIPseq/ChIPseq_peaks/npy/'. The generated file iss named e.g., 'Cbx3_peaks_chr1.npy'.

### 2. From the .bedgraph files of ChIP-seq tracks obtain the .npy normalized array that later will be combined with Hi-C
Use the script '02_ChIpseq_tracks_fromBedgraph_toNormalized.py' contained in the subfolder './code/ChIPseq_prepro/'.  
The input is .bedgraph files, one for each factor you study. See the input I used in the application in NPCs in the subfolder './data/NPC/ChIPseq/ChIPseq_tracks/bedgraph/', where the format of the name of files is '(name-of-the-factor).bedgraph'.  
This scripts involves several steps, where intermediate outputs are stored in the subfolder '.data/NPC/ChIPseq/ChIPseq_tracks/', in case the run gets interrupted and you don't want to start from the beginning. You can later delete these intermediate .npy arrays. DO NOT delete the final .npy arrays, saved in './data/NPC/ChIPseq/ChIPseq_tracks/', whose name ends with '_norm_01.npy'. They are the ChIP-seq tracks normalized in the range 0-1, which will be next combined with Hi-C.  
In the subfolder './data/NPC/ChIPseq/ChIPseq_tracks/' you can find the final data I obtained in the application of CHROMATIC to 18 factors in NPCs.

### 3. Optional: Plot the normalized ChIP-seq tracks in a region of interest
Use the script '03_plot_normalizedChIPseqs.py' contained in the subfolder './code/ChIPseq_prepro/' to plot the normalized ChIP-seq tracks in a region of interest. The script is set to perform this plot on the HoxA locus of the mouse genome mm10.

## Hi-C data pre-processing

### 4. Transform the intra-chromosomal Hi-C maps into .npz scipy sparse arrays
Use the script '04_HiC_NPs_from_decay_to_decaymedian_chr1.py' contained in the subfolder './code/HiC_prepro/'. You may create a script like this one for each chromosome (good for speed), or uncomment one line as indicated in the code to perform the same script for all chromosomes (bad for speed).
The input is the intra-chromosomal Hi-C map, in a sparse format. In our application, we work at 5kb resolution, the input map for chr1 is in the subfolder './data/NPC/HiC/abc/', and the name of the file is 'HiC_chr1_5kb.abc'. It is stored in a sparse, bed-like format, where each line describes the pair of interacting bins and their Hi-C value. Note that the bins are not indicated in base pairs, but the number of the bin is indicated. Thus, bin 0 stands for the region 0-4,999 bp, bin 1 represents the region 5,000-9,999 bp, and so on. A generic example of one line in the 'HiC_chr1_5kb.abc' file is: bin1 bin2 HiCvalue12.
The script generates a .npz compressed array with scipy.sparse package, that optimizes storage and speed. We apply a median filter to the Hi-C map, to reduce noise and enhance CHROMATIC performance. The obtained map is stored in the subfolder './data/NPC/HiC/npz/' with the name 'HiC_chr1_median.npz'. Keep this file, since it will be combined with ChIP-seq in the next steps.

## CHROMATIC

### 5. Combine ChIP-seq and Hi-C data to identify the CHROMATIC interaction maps 
Use the script '05_CHROMATIC_NPs_Cbx3_chr1.py' contained in the subfolder './code/CHROMATIC/' to integrate ChIP-seq and Hi-C data. You may create a script like this one for each chromosome (good for speed), or uncomment one line as indicated in the code to perform the same script for all chromosomes (bad for speed). Also, this scripts performs CHROMATIC only for CBX3, again for the purpose of performance. Thus, you may create a script like this one for each chromosome and for each studied factor.  
This scripts takes as input the normalized ChIP-seq tracks (produced by the script '02_ChIpseq_tracks_fromBedgraph_toNormalized.py') and the .npy ChIP-seq peaks (produced by the script '01_ChIpseq_peaks_fromBed_toNpy.py').  
A sliding window  of 10 Mb is used to improve performance. It is indicated in the line 'w = 10000000/res', you can change this parameter (carefully to not excessively increase the computation time).  
The distribution of the ChIP-seq values is "stretched" by a specific transform that is reported in the paper as a Supplementary Figure. The goal is that lower ChIP-seq values are lowered further, and high ChIP-seq values are increased further. In image processing, this operation corresponds to contrast enhancement.  
The script generates to compressed-sparse arrays:  
1. the CHROMATIC interaction map (one for each factor and for each chromosome, saved in './data/NPC/CHROMATIC/int_maps/', see e.g., 'Cbx3_chr1_CHROMATIC_int_map.npz')
2. the "matrix of ChIP-seq peaks", a matrix that has value=1 in the pixel corresponding to bins i and j, if the bins i-1,i,i+1 OR the bins j-1,j,j+1 contain at least one ChIP-seq peak, and has value=0 otherwise (one for each factor and for each chromosome, saved in './data/NPC/ChIPseq/ChIPseq_peaks/npz/', see e.g., 'Cbx3_peaks_chr1.npz'). This matrix will be useful in the next step of detecting loops/hubs.  

### 6. Detect loops/hubs in the CHROMATIC interaction maps 
Use the script '06_CHROMATIC_NPs_peaks_Cbx3_chr1.py' contained in the subfolder './code/CHROMATIC/' to detect loops and hubs associated with the factor of interest. You may create a script like this one for each chromosome (good for speed), or uncomment one line as indicated in the code to perform the same script for all chromosomes (bad for speed). Also, this scripts performs CHROMATIC-peak-detection only for CBX3, again for the purpose of performance. Thus, you may create a script like this one for each chromosome and for each studied factor.  
The input of this step is the output of the previous step ('05_CHROMATIC_NPs_Cbx3_chr1.py'):
1. the CHROMATIC interaction map (one for each factor and for each chromosome, saved in './data/NPC/CHROMATIC/int_maps/', see e.g., 'Cbx3_chr1_CHROMATIC_int_map.npz')
2. the matrix of ChIP-seq peaks (one for each factor and for each chromosome, saved in './data/NPC/ChIPseq/ChIPseq_peaks/npz/', see e.g., 'Cbx3_peaks_chr1.npz').  
Again, a sliding window is used. We used a window of 10 Mb at 5 kb resolution (of 2000 bins). The size of the window is indicated in the line 'w = 10000000/res'. You can change this parameter (carefully to not excessively increase the computation time). It must have the same value as used in '05_CHROMATIC_NPs_Cbx3_chr1.py'.  
In brief, the script thresholds the CHROMATIC interaction map, multiplies it by the matrix of ChIP-seq peaks (this operation ensures that for each detected interaction there is at least one ChIP-seq peak at one of the two interacting loci), and performs a series of operations from morphological image processing (opening, closing, dilation) to produce a binary matrix, that has value=1 in correspondance of detected loops/hubs, and 0 elsewhere.  
This matrix is saved in './data/NPC/CHROMATIC/detections/' and represents the locations of the detected interactions for the considered factors, within a certain chromosome (see e.g., 'Cbx3_chr1_CHROMATIC_peaks.npz').  






