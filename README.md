# Project Description

This project conducts and reproduces analyses presented in Marisa et al. Gene Expression
Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic
Value. PLoS Medicine, May 2013. Data processing and analysis was performed on the BU Shared Computing Cluster 
using R version 4.1.2.

# Contributors

Data Curator: Urvy Mudgal

Programmer: Jason Yeung

Analyst: Manas Dhanuka

# Repository Contents

**preprocessing.R**

Data preprocessing and quality control

Normalization

RLE/NUSE

Batch effect correction

PCA

**Noise_Filtering_and_Dimension_reduction.R**
R script to perfrom the following filters and return  a csv file for heirarchal clustering.
Expressed in at least 20% of samples (i.e. for each gene, at least 20% of the gene-expression values must be > log2(15)
Have a variance significantly different from the median variance of all probe sets using a threshold of p<0.01.
Have a coefficient of variation > 0.186


**Heirarchal_clustering_and_subtype_disc.R**
R script to perform hierarchical clustering on filtered data from **Noise_Filtering_and_Dimension_reduction.R**. 
Plot a heatmap of the gene-expression of each gene across all samples grouping data into C3 molecular cancer subtype or others. 
Identifying genes differentially expressed between the two clusters using a Welch t-test with a p adjusted p<0.05.
Select the most differentially expressed genes that you feel best define the clusters and explain your selection.

