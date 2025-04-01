# Crohn’s disease subtypes identified by single cell RNA-sequencing signatures across treatment naïve intestinal samples 

## This project examines the variability of treatment naive Crohn's disease across donors and the small and large intestines using single-cell RNA sequencing. A clustering stability assessment workflow was developed to examine repeatability and robustness of clustering and downstream results. 

## Methods: 

The general clustering stability assessment workflow can be found in the template directory. The stability assessment can be performed in 4 steps. 

1. Perform the clustering stability assessment by randomly subsetting the whole data and clustering each subset. 
- This can be repeated n times. 

2. Compare the subset clustering results (query) to the whole data clustering results (reference). 
- This is accomplished using Seurat's TransferData function.
- Repeat for each query.

3. Assess clustering similarity by re-clustering the stably assigned cells with different clustering parameters. Then use the Adjusted Rand Index (computed with the R package Dune) to quantify similarity. 

4. Evaluate clustering repeatability of clustering results and downstream results using the Jaccard Index and Spearman Correlation.
- This can be performed comparing the stably assigned cells to the reference dataset. 
- Jaccard Index is used to evaluate the similarity of cell type annotation across different clustering results.
- Spearman correlation is used to evaluate the similarity of downstream results across different clustering results. 

## Downstream Analysis

Ileum, colon, and rectum specific analyses can be found in the respective directories. This includes the clustering stability assessment as well as downstream analysis. Code for the integrative GWAS summary statistic and single cell RNA seq analysis was adapted from [Duncan et al.](https://github.com/Integrative-Mental-Health-Lab/linking_cell_types_to_brain_phenotypes.git).


  

 
