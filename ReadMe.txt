Instructions for the codes:
1. Run the shell script first to obtain TPM matrix of fastq data of SRP118628/GSE104140 
(raw data can be downloaded from SRA database . https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP118628&o=acc_s%3Aa, sample_list.txt contains file names of fastq file)

(Codes below were run in R software)

2. Merge TPM matrix of SRP118628/GSE104140 and expression matrix of GSE28829 via common gene symbol using "VLOOKUP" function in excel (or merge() function in R) 
(The probe matrix and platform data table of GSE28829 can be downloaded from GEO database. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28829)
3. Conduct batch normlization for the merged matrix to get a batch normalized version of the merged matrix. 
4. Conduct DEG analysis and draw a volcano plot and heatmap for DEGs. 
(Grouping information for advanced plaque and early plaque of samples can be copied from GEO database and SRA database)
5. Conduct WGCNA analysis. 
6. Conduct Functional enrichment analysis for gene module
7. Conduct lasso regression and LDA analysis for MCODE genes obained from Cytoscape software. 

Note: 
a. The input files of each step were described in detail in the annotation of R code (annotations just below codes reading input file into R). 
b. For validation set GSE43292, the data acquisition is the same as GSE28829 except for the difference in accession number.