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
However, we uploaded also the data that define the 3D-types that we identified genome-wide, as explained below.

## ChIP-seq data pre-processing

### 1. Transform the ChIP-seq peaks .bed file into .npy array
Use the script '01_ChIpseq_peaks_fromBed_toNpy.py' contained in the subfolder './code/0_ChIPseq_prepro/' to generate .npy arrays from the .bed files of ChIP-seq peaks locations.  
The input is a .bed file, one for each factor you study. See the input I used in the application in NPCs in the subfolder './data/NPC/ChIPseq/ChIPseq_peaks/bed/', where the format of the name of files is '(name-of-the-factor)_peaks_ucsc.bed'.  
The ouput is one .npy array for every factor and every chromosome, saved in the subfolder './data/NPC/ChIPseq/ChIPseq_peaks/npy/'. The generated file iss named e.g., 'Cbx3_peaks_chr1.npy'.

### 2. From the .bedgraph files of ChIP-seq tracks obtain the .npy normalized array that later will be combined with Hi-C
Use the script '02_ChIpseq_tracks_fromBedgraph_toNormalized.py' contained in the subfolder './code/0_ChIPseq_prepro/'.  
The input is .bedgraph files, one for each factor you study. See the input I used in the application in NPCs in the subfolder './data/NPC/ChIPseq/ChIPseq_tracks/bedgraph/', where the format of the name of files is '(name-of-the-factor).bedgraph'.  
This scripts involves several steps, where intermediate outputs are stored in the subfolder '.data/NPC/ChIPseq/ChIPseq_tracks/', in case the run gets interrupted and you don't want to start from the beginning. You can later delete these intermediate .npy arrays. DO NOT delete the final .npy arrays, saved in './data/NPC/ChIPseq/ChIPseq_tracks/', whose name ends with '_norm_01.npy'. They are the ChIP-seq tracks normalized in the range 0-1, which will be next combined with Hi-C.  
In the subfolder './data/NPC/ChIPseq/ChIPseq_tracks/' you can find the final data I obtained in the application of CHROMATIC to 18 factors in NPCs.

### 3. Optional: Plot the normalized ChIP-seq tracks in a region of interest
Use the script '03_plot_normalizedChIPseqs.py' contained in the subfolder './code/0_ChIPseq_prepro/' to plot the normalized ChIP-seq tracks in a region of interest. The script is set to perform this plot on the HoxA locus of the mouse genome mm10.

## Hi-C data pre-processing

### 4. Transform the intra-chromosomal Hi-C maps into .npz scipy sparse arrays
Use the script '04_HiC_NPs_from_decay_to_decaymedian_chr1.py' contained in the subfolder './code/1_HiC_prepro/'. 
You may create a script like this one for each chromosome (good for speed), or uncomment one line as indicated in the code to perform the same script for all chromosomes (bad for speed).  
The input is the intra-chromosomal Hi-C map, in a sparse format.  
In our application, we work at 5kb resolution, the input map for chr1 is in the subfolder './data/NPC/HiC/abc/', and the name of the file is 'HiC_chr1_5kb.abc'.  
It is stored in a sparse, bed-like format, where each line describes the pair of interacting bins and their Hi-C value. Note that the bins are not indicated in base pairs, but the number of the bin is indicated. Thus, bin 0 stands for the region 0-4,999 bp, bin 1 represents the region 5,000-9,999 bp, and so on. A generic example of one line in the 'HiC_chr1_5kb.abc' file is: bin1 bin2 HiCvalue12.  
The script generates a .npz compressed array with scipy.sparse package, that optimizes storage and speed.  
We apply a median filter to the Hi-C map, to reduce noise and enhance CHROMATIC performance.  
The obtained map is stored in the subfolder './data/NPC/HiC/npz/' with the name 'HiC_chr1_median.npz'. Keep this file, since it will be combined with ChIP-seq in the next steps.


## CHROMATIC

### 5. Combine ChIP-seq and Hi-C data to identify the CHROMATIC interaction maps 
Use the script '05_CHROMATIC_NPs_Cbx3_chr1.py' contained in the subfolder './code/2_CHROMATIC/' to integrate ChIP-seq and Hi-C data.  
You may create a script like this one for each chromosome (good for speed), or uncomment one line as indicated in the code to perform the same script for all chromosomes (bad for speed).  
Also, this scripts performs CHROMATIC only for CBX3, again for the purpose of performance. Thus, you may create a script like this one for each chromosome and for each studied factor.  
This scripts takes as input the normalized ChIP-seq tracks (produced by the script '02_ChIpseq_tracks_fromBedgraph_toNormalized.py') and the .npy ChIP-seq peaks (produced by the script '01_ChIpseq_peaks_fromBed_toNpy.py').  
A sliding window  of 10 Mb is used to improve performance. It is indicated in the line 'w = 10000000/res', you can change this parameter (carefully to not excessively increase the computation time).  
The distribution of the ChIP-seq values is "stretched" by a specific transform that is reported in the paper as a Supplementary Figure. The goal is that lower ChIP-seq values are lowered further, and high ChIP-seq values are increased further. In image processing, this operation corresponds to contrast enhancement.  
The script generates to compressed-sparse arrays:  
1. the CHROMATIC interaction map (one for each factor and for each chromosome, saved in './data/NPC/2_CHROMATIC/int_maps/', see e.g., 'Cbx3_chr1_CHROMATIC_int_map.npz')
2. the "matrix of ChIP-seq peaks", a matrix that has value=1 in the pixel corresponding to bins i and j, if the bins i-1,i,i+1 OR the bins j-1,j,j+1 contain at least one ChIP-seq peak, and has value=0 otherwise (one for each factor and for each chromosome, saved in './data/NPC/ChIPseq/ChIPseq_peaks/npz/', see e.g., 'Cbx3_peaks_chr1.npz'). This matrix will be useful in the next step of detecting loops/hubs.  

### 6. Detect loops/hubs in the CHROMATIC interaction maps 
Use the script '06_CHROMATIC_NPs_peaks_Cbx3_chr1.py' contained in the subfolder './code/2_CHROMATIC/' to detect loops and hubs associated with the factor of interest. 
You may create a script like this one for each chromosome (good for speed), or uncomment one line as indicated in the code to perform the same script for all chromosomes (bad for speed).  
Also, this scripts performs CHROMATIC-peak-detection only for CBX3, again for the purpose of performance. Thus, you may create a script like this one for each chromosome and for each studied factor.  
The input of this step is the output of the previous step ('05_CHROMATIC_NPs_Cbx3_chr1.py'):
1. the CHROMATIC interaction map (one for each factor and for each chromosome, saved in './data/NPC/CHROMATIC/int_maps/', see e.g., 'Cbx3_chr1_CHROMATIC_int_map.npz')
2. the matrix of ChIP-seq peaks (one for each factor and for each chromosome, saved in './data/NPC/ChIPseq/ChIPseq_peaks/npz/', see e.g., 'Cbx3_peaks_chr1.npz').  
Again, a sliding window is used. We used a window of 10 Mb at 5 kb resolution (of 2000 bins).  
The size of the window is indicated in the line 'w = 10000000/res'.  
You can change this parameter (carefully to not excessively increase the computation time). It must have the same value as used in '05_CHROMATIC_NPs_Cbx3_chr1.py'.  
In brief, the script thresholds the CHROMATIC interaction map, multiplies it by the matrix of ChIP-seq peaks (this operation ensures that for each detected interaction there is at least one ChIP-seq peak at one of the two interacting loci), and performs a series of operations from morphological image processing (opening, closing, dilation) to produce a binary matrix, that has value=1 in correspondance of detected loops/hubs, and 0 elsewhere.  
This matrix is saved in './data/NPC/CHROMATIC/detections/' and represents the locations of the detected interactions for the considered factors, within a certain chromosome (see e.g., 'Cbx3_chr1_CHROMATIC_peaks.npz').  


## Analysis of the CHROMATIC output

If you have detected peaks of interactions (loops/hubs) for at least one factor and at least on one chromosome, you may now want to analyse them in terms of different aspects.  
In the subfolder './code/3_analysis_CHROMATIC_output/' you find some scripts to measure: the genomic distance of the detected interactions, their size, the Hi-C intensity in correspondance of the detection (before their combination with ChIP-seq), and the CHROMATIC score contained within the CHROMATIC interaction map (after the combination of Hi-C and ChIP-seq data) in correspondance of the detections.

### 7. Genomic distance of the detected interactions
Use the script '07_genomic_distance_NPs_peaks_Cbx3_allchrs.py' contained in the subfolder './code/3_analysis_CHROMATIC_output/' to obtain a numpy array of the genomic distance of the detected interactions.  
This script processes all chromosomes. If you want to process only one chromosome, you can easily do it as indicated within the script.  
Also, the script processes only 1 factor (CBX3 in our example). If you want to process all factors, you can follow the indications within the script.  
For each chromosome, the matrix of the detected interactions is loaded (e.g., 'Cbx3_chr1_CHROMATIC_peaks.npz', which was generated by the script '06_CHROMATIC_NPs_peaks_Cbx3_chr1.py').  
As described above, this is a binary matrix. Each label found in the matrix corresponds to a patch of detected interactions (e.g., a hub). The script computes the average distance between the interacting loci and stores this value in a numpy array, which is saved in './data/NPC/CHROMATIC/detections_analysis/' (see 'genomic_distance_CHROMATIC_peaks_Cbx3_chr1.npy' as an example).

### 8. Size of the detected interactions
Use the script '08_size_NPs_peaks_Cbx3_allchrs.py' contained in the subfolder './code/3_analysis_CHROMATIC_output/' to obtain a numpy array of the size of the detected interactions.  
This script processes all chromosomes. If you want to process only one chromosome, you can easily do it as indicated within the script.  
Also, the script processes only 1 factor (CBX3 in our example). If you want to process all factors, you can follow the indications within the script.  
For each chromosome, the matrix of the detected interactions is loaded (e.g., 'Cbx3_chr1_CHROMATIC_peaks.npz', which was generated by the script '06_CHROMATIC_NPs_peaks_Cbx3_chr1.py').  
As described above, this is a binary matrix. Each label found in the matrix corresponds to a patch of detected interactions (e.g., a hub). The script computes the number of pixels corresponding to each label/detection and stores this value in a numpy array, which is saved in './data/NPC/CHROMATIC/detections_analysis/' (see 'size_CHROMATIC_peaks_Cbx3_chr1.npy' as an example).  
In our application case, the size is in terms of 5kbx5kb pixels.

### 9. Hi-C intensity in correspondance of the detected interactions
Use the script '09_hic_intensity_NPs_peaks_Cbx3_allchrs.py' contained in the subfolder './code/3_analysis_CHROMATIC_output/' to obtain a numpy array of the Hi-C intensity in correspondance of the detected interactions (before their combination with ChIP-seq).  
This script processes all chromosomes. If you want to process only one chromosome, you can easily do it as indicated within the script.  
Also, the script processes only 1 factor (CBX3 in our example). If you want to process all factors, you can follow the indications within the script.  
For each chromosome, the matrix of the detected interactions is loaded (e.g., 'Cbx3_chr1_CHROMATIC_peaks.npz', which was generated by the script '06_CHROMATIC_NPs_peaks_Cbx3_chr1.py').  
As described above, this is a binary matrix. Each label found in the matrix corresponds to a patch of detected interactions (e.g., a hub). The script computes the average original Hi-C value in the pixels of the detection and stores this value in a numpy array, which is saved in './data/NPC/CHROMATIC/detections_analysis/' (see 'hic_intensity_CHROMATIC_peaks_Cbx3_chr1.npy' as an example).  

### 10. CHROMATIC score in correspondance of the detected interactions
Use the script '10_ChromaticScore_NPs_peaks_Cbx3_allchrs' contained in the subfolder './code/3_analysis_CHROMATIC_output/' to obtain a numpy array of the CHROMATIC score contained within the CHROMATIC interaction map (after the combination of Hi-C and ChIP-seq data) in correspondance of the detections. 
This script processes all chromosomes. If you want to process only one chromosome, you can easily do it as indicated within the script.  
Also, the script processes only 1 factor (CBX3 in our example). If you want to process all factors, you can follow the indications within the script.  
For each chromosome, the matrix of the detected interactions is loaded (e.g., 'Cbx3_chr1_CHROMATIC_peaks.npz', which was generated by the script '06_CHROMATIC_NPs_peaks_Cbx3_chr1.py').  
As described above, this is a binary matrix. Each label found in the matrix corresponds to a patch of detected interactions (e.g., a hub). The script computes the average of the CHROMATIC score in the pixels of the detection and stores this value in a numpy array, which is saved in './data/NPC/CHROMATIC/detections_analysis/' (see 'chromatic_score_CHROMATIC_peaks_Cbx3_chr1.npy' as an example).  


## Latent Semantic Analysis to find types of 3D interactions

### 11. Sum of the CHROMATIC detections made for all factors 
Use the script '11_sum_Chromatic_NPs_peaks_allfactors.py' contained in the subfolder './code/4_LSA_3Dtypes/' to compute the matrix of the overlap of the CHROMATIC detections of all factors ("overlap-matrix").  
This step is important for the subsequent ones, which will perform Latent Semantic Analysis (LSA) to identify types of 3D interactions.  
The script processes all chromosome and all factors.  
For each chromosome, a new matrix P_sum is produced, which is:  
P_sum[i,j]=0, if the pixel [i,j] does not contain any detection made for any of the studied factors,  
P_sum[i,j]=1, if the pixel [i,j] contains a detection for 1 factor,  
P_sum[i,j]=2, if the pixel [i,j] contains a detection for 2 factors,  
and so on.  
This overlap-matrix is saved in './data/NPC/CHROMATIC/3Dtypes/' (with the name 'sum_CHROMATIC_peaks_allfactors_chr1.npz', for chromosome 1).

### 12. Prepare the input for LSA
Use the script '12_lsa_NPs_prepare_input_chr1.py' contained in the subfolder './code/4_LSA_3Dtypes/' to prepare the input that is needed for LSA.  
The script processes all factors in chromosome 1.  
You may create a script like this one for each chromosome (good for speed), or uncomment one line as indicated in the code to perform the same script for all chromosomes (bad for speed).    
The input of the script is the matrix of the overlap of the CHROMATIC detections of all factors (saved in './data/NPC/CHROMATIC/3Dtypes/', with the name 'sum_CHROMATIC_peaks_allfactors_chr1.npz', for chromosome 1) and the CHROMATIC detections made for each factor (saved in './data/NPC/CHROMATIC/detections/', with the name 'Cbx3_chr1_CHROMATIC_peaks.npz' for CBX3 in chromosome 1).  
The ouput is a pandas array (pd.Series), which has a row for each non-zero pixel of the overlap-matrix (corresponding to a detected interaction for one factor or more), where the row contains the name of the factors participating in that same interaction.  
For example, if P_sum[i,j]=1 because factor CTCF has a detected interaction in the pixel [i,j], the corresponding line in the pandas array is 'CTCF'.  
Otherwise, if P_sum[i,j]=2 because factors H3K27ac and SMC1 have a detected interaction in the pixel [i,j], the corresponding line in the pandas array is 'H3K27ac SMC1'.  
The generated array is saved in './data/NPC/CHROMATIC/3Dtypes/' (with the name '3Dint-factor_lsa_chr1.npy', for chromosome 1).  

### 13. Perform LSA to identify the types of 3D interactions
Use the script '13_lsa_NPs.py' contained in the subfolder './code/4_LSA_3Dtypes/' to perform LSA genome-wide for all the studied factors.  
The input of the script is the output of step no.12, performed for all chromosomes, saved in './data/NPC/CHROMATIC/3Dtypes/'.  
The output is saved in './data/NPC/CHROMATIC/3Dtypes/' and is:  
1. a compressed array describing to what 3D type each pixel with a detected CHROMATIC interaction belongs. Each pixel with a detected CHROMATIC interaction corresponds to a row, and the row contains its enrichment in the identified 3D type. Each pixel belongs to the 3D type for which corresponding to the maximum value of the row. It is saved as '3Dint-3Dtype_allchrs.npz'. 
2. an array that defines the identified types of 3D interactions (3D types), based on specific combinations of the factors in 3D. 3D types correspond to rows and factors correspond columns. Thus, each row describes the composition of the 3D type in terms of enrichment (positive score) or depletion (negative score) of the factors. It is saved as '3Dtype-factor_allchrs.npy'.  
The script also prints its computing time.  

### 14. Plot the 3D types defined by LSA
Use the script '14_plot_output_LSA_NPC.py' contained in the subfolder './code/4_LSA_3Dtypes/' to plot the heatmap that defines 3D types.  
The input is the output of step no.13 '3Dtype-factor_allchrs.npy', saved in './data/NPC/CHROMATIC/3Dtypes/'.  
The script generates a heatmap where each row corresponds to a 3D type, defined in terms of enrichment (in red) or depletion (in blue) of the factors.  
On top, colors express the functional role of factors. Adjust the palette according to the functional role of your factors.  
Also, a dendogram done by unsupervised hierarchical clustering describes the association between factors based on the 3D types.  
The plot is saved as .pdf in './data/NPC/CHROMATIC/3Dtypes/' as 'LSA_clusterheatmap_NPC.pdf'.


### 15. Count the number of pixels/3D interactions corresponding to each identified 3D type
Use the script '15_count_3Dinteractions_per_3Dtype_NPs.py' contained in the subfolder './code/4_LSA_3Dtypes/' to count how many pixels correspond to each 3D type.  
The input was obtained at step no.13.  
The output is saved in './data/NPC/CHROMATIC/3Dtypes/' as '3Dtype_counter_allchrs.npy'. It is an array whose length is equal to the number of 3D types identified (which is N-1 and is the maximum number of 3D types that is possible to obtain with LSA, having N factors).  
The script prints the total number of pixels (genome-wide) with a detected CHROMATIC interaction

### 16. Get the 1D loci (on the linear chromatin) participating in the identified 3D types
Use the script '16_generate_bedfile_NPs_3Dtype0.py' contained in the subfolder './code/4_LSA_3Dtypes/' to generate a .bed file of the 1D loci that interact in 3D through 3D type-0.  
You may create a script like this one for each 3D type (good for speed), or uncomment one line as indicated in the code to perform the same script for all 3D types (bad for speed).  
The input is:  
1. one of the two output of step no.13, named '3Dint-3Dtype_allchrs.npz' and saved in './data/NPC/CHROMATIC/3Dtypes/'
2. the overlap-matrix generated at step no.11 for each chromosome, saved in './data/NPC/CHROMATIC/3Dtypes/' (with the name 'sum_CHROMATIC_peaks_allfactors_chr1.npz', for chromosome 1).  
The output is saved in './data/NPC/CHROMATIC/3Dtypes/' with the name '3Dtype_0_allchrs.bed' for 3D type-0. 


##  Identify the major types of 3D interaction
In order to better capture the biological meaning of 3D types, we mapped 3D types into lists of 1D loci (as in step no.16) and measured their overlap with the functional genomic features of active enhancers (AE), active promoters (AP), poised enhancers (PE), bivalent promoters (BP), super-enhancers (SE), and constitutive LADs (CL). 
This overlap was measured with the command *intersect -v* from the BEDTools toolkit (Quinlan, A.R. & Hall, I.M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841-2 (2010)).  
As an example, you can find our data of the overlap between '3Dtype_0_allchrs.bed' and AE in the subfolder './data/NPC/CHROMATIC/major3Dtypes/' with the name '3Dtype0_AE_overlap.bed'. We counted the number of entries in this file and we repeated the same for all 3D types, and for AP, PE, BP, SE, and CL.  
Then, we used the following scripts to compute the Log of the Odds Ratio of each overlap, use this data as input for Principal Component Analysis and K-means unsupervised clustering, in order to find the major types of 3D interactions, given by clusters of the previously identified 3D types.  
You can follow this process or you can modify it to include other functional genomic features in your classification (beyond AE, AP, etc.).  

### 17. Clustering of 3D types into major 3D types  
Use the script '17_from3Dtypes_to_major3Dtypes.py' contained in the subfolder './code/5_major3Dtypes/' to obtain the clusters of 3D types corresponding to the major 3D types.  
Prepare the input of the overlap of functional genomic features with 3D types as described above.  
In the part of the script of PCA, follow the comments to choose the number of PCs that explain at least the 80% of the total variance. In our application case in NPC, we set *PCA(n_components=3)* because the first 3 PCs explained the 90.88% of the total variance. Then, PCA is performed with the chosen number of components and the obtained data is used for K-means clustering.  
To determine the number of clusters to compute, the K-means algorithm is run multiple times with a different number of clusters (from 1 to 17 in NPC, where 17 3D types were classified). For each solution, the Within Cluster Sum of Squares (WCSS) is computed. To determine the number of clusters to use, we use the approach known as the Elbow method, which consists of looking for a kink or elbow in the plot of the values of WCSS against the number of clusters (which is produced by the script in the part of K-means, and saved in './data/NPC/CHROMATIC/major3Dtypes/' as 'kink_NPC.pdf'). The elbow point is identified by the different exponential of the descent on the left and the right of the plot. In our application in NPC, the elbow appeared in correspondence with 4 clusters of major types of interactions, so we do *KMeans(n_clusters=4)*.  
Finally, the script plots the identified clusters of 3D types on the PC1 and PC2. Even though it is not included in this script, we could associate the 4 major 3D types with their functional outcome, by studying:  
1. their enrichment/depletion in AE, AP, PE, BP, SE, CL
2. looking at their corresponding genes, the proportion of genes that are highly-, lowly-expressed or silent
3. the proportion of their overlap with A/B compartments.  
Thanks to this, we called the obtained 4 major 3D types as Active, Neuronal TFs, PcG-Bivalent and Inactive, as it appears in the final plot generated by the script '17_from3Dtypes_to_major3Dtypes.py', saved in './data/NPC/CHROMATIC/major3Dtypes/' as 'PCA_NPC.pdf'. You can modify names, colors and other parameters according to your findings.  

### 18. Count pixels associated with the major 3D types and annotate the contact maps with the major 3D types
Use the script '18_annotate_pixels_major3Dtypes_NP.py' contained in the subfolder './code/5_major3Dtypes/' to count the number of pixels classified in each major 3D type and annotate the contact maps with the major 3D types.
The input is:  
1. one of the two output of step no.13, named '3Dint-3Dtype_allchrs.npz' and saved in './data/NPC/CHROMATIC/3Dtypes/'
2. the overlap-matrix generated at step no.11 for each chromosome, saved in './data/NPC/CHROMATIC/3Dtypes/' (with the name 'sum_CHROMATIC_peaks_allfactors_chr1.npz', for chromosome 1).  
In our application to NPC, we found that the optimal number of clusters (major 3D types) is 4, thus in the script we set *n_clu = 4*. Change this parameter according to the number of clusters you identified.  
Also, change the parameters relative to which 3D types compose every major 3D type according to your findings.  
The output is:  
1. a compressed array for each chromosome, representing intra-chromosomal contact maps annotating the identified major 3D types, saved in './data/NPC/CHROMATIC/major3Dtypes/' (with the name 'painted_hicmap_major3Dtypes_chr1.npz' for chromosome 1). In this folder you can find an example of it, for our application in NPC. Since in NPC we have 4 major clusters, the generated map can assume value 0, 1, 2, 3, 4, where each value corresponds to a different classification: 0:unclassified, 1:nTFs, 2:Inactive, 3:Active, 4:PcG-Bivalent. See step no.20 for more details about this map. 
2. a numpy array that describes the number of pixels associated to each 3D type, saved in './data/NPC/CHROMATIC/major3Dtypes/' as 'counter_pixels_major3Dtypes_NP.npy'.  
The script also prints the total number of classified pixels genome-wide. 

### 19. Optional: Pie-chart of the major 3D types
Use the script '19_plot_major3Dtypes_piechart.py' contained in the subfolder './code/5_major3Dtypes/' to plot with a pie-chart the number of pixels associated with the major 3D types.  
The input is the output of step no.18, saved in './data/NPC/CHROMATIC/major3Dtypes/' as 'counter_pixels_major3Dtypes_NP.npy'.  
The plot is saved in './data/NPC/CHROMATIC/major3Dtypes/' as 'piechart_major3Dtypes_NPC.pdf'.

### 20. Optional: Paint the contact maps with the major 3D types
To look at specific loci, you can use the script '20_paint_hicmap_major3Dtypes_NP.py' contained in the subfolder './code/5_major3Dtypes/' to "paint" your contact maps with the classified major 3D types, each one with a different color.  
The input is the intra-chromosomal contact maps annotating the identified major 3D types generated in step no.18, saved in './data/NPC/CHROMATIC/major3Dtypes/'.  
In the script we plot the major 3D types in the region chr18:53800000-56600000 containing the gene Zfp608. Thus, the input map is 'painted_hicmap_major3Dtypes_chr18.npz'. You run this script for several different regions as indicated in the comments.   
The script is currently designed to plot 4 major 3D types, so, if you have a different number of major 3D types, you have to modify some parameters, such as the colormap and *vmax* in the *imshow* command.  
The output is saved in './data/NPC/CHROMATIC/major3Dtypes/' as 'painted_hicmap_Zfp608.png', plotting the 3D interactions that were classified in different major 3D typess with different colors. The white background is set to transparent to allow you to overlap this map with the original Hi-C map, and discern which Hi-C interactions have finally been classified.  











