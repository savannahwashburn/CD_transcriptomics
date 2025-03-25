#!/usr/bin/env Rscript

#script to create top results for conditional analysis - rectum

#load libraries 
library(readr)
library(dplyr)


#---------------rectum results--------------#
#read in specificity matrix 
df <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/rectum_conti_specificity_matrix.txt")

#read in step 2 results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_rectum_conti_spe_2025.gsa.out", skip=4) %>%
  arrange(P) %>%
  filter(P < 0.05/n()) #keep significant clusters after bonferroni correction

sig_list <- df_marg$VARIABLE

df_sig <- df[c("GENE", sig_list)] #Monocyte

#save the output 
write_tsv(df_sig, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_rectum_conti_matrix_sig.txt")


#relax thresholds and keep cell types with p value less than 0.05

df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_rectum_conti_spe_2025.gsa.out", skip=4) %>%
  arrange(P) %>%
  filter(P < 0.05) #keep significant clusters 

sig_list <- df_marg$VARIABLE

df_sig <- df[c("GENE", sig_list)] #Monocyte, mast cell, CD4 t cell, ILC3, CD8 T cell, Ambig T cell, Macrophage

#save the output 
write_tsv(df_sig, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_rectum_conti_matrix_relax_sig.txt")
