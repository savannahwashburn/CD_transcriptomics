#!/usr/bin/env Rscript

#script to create top results for conditional analysis - colon

#load libraries 
library(readr)
library(dplyr)


#---------------colon results--------------#
#read in specificity matrix 
df <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/colon_conti_specificity_matrix.txt")

#read in step 2 results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_colon_conti_spe_2025.gsa.out", skip=4) %>%
  arrange(P) %>%
  filter(P < 0.05/n()) #keep significant clusters after bonferroni correction

sig_list <- df_marg$VARIABLE

df_sig <- df[c("GENE", sig_list)] #Monocyte, Treg, and Mast cell 

#save the output 
write_tsv(df_sig, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_colon_conti_matrix_sig.txt")

