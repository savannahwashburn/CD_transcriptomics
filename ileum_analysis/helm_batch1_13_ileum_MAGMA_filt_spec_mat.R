#!/usr/bin/env Rscript

#script to create top results for conditional analysis

#load libraries 
library(readr)
library(dplyr)


#---------------ileum results--------------#
#read in specificity matrix 
df <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/conti_specificity_matrix.txt")

#read in step 2 results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_conti_2025.gsa.out", skip=4) %>%
  arrange(P) %>%
  filter(P < 0.05/n()) #keep significant clusters after bonferroni correction

sig_list <- df_marg$VARIABLE

df_sig <- df[c("GENE", sig_list)] #filter to keep Pro-inflammatory monocytes, CD4 T cells and monocytes 

#save the output 
write_tsv(df_sig, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_conti_matrix_sig.txt")


#-----------high resolution stably assigned cells--------------#

#read in specificity matrix 
df <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/ileum_high_res_conti_specificity_matrix.txt")

#read in step 2 results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_high_res_conti_spe_2025.gsa.out", skip=4) %>%
  arrange(P) %>%
  filter(P < 0.05/n()) #keep significant clusters after bonferroni correction

sig_list <- df_marg$VARIABLE

df_sig <- df[c("GENE", sig_list)] #filter to keep Pro-inflammatory monocytes, Treg and CD8 T cell 

#save the output 
write_tsv(df_sig, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_high_res_conti_matrix_sig.txt")


#------------reference--------------#

#read in specificity matrix 
df <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/ileum_ref_conti_specificity_matrix.txt")

#read in step 2 results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_ref_conti_spe_2025.gsa.out", skip=4) %>%
  arrange(P) %>%
  filter(P < 0.05/n()) #keep significant clusters after bonferroni correction

sig_list <- df_marg$VARIABLE

df_sig <- df[c("GENE", sig_list)] #filter to keep Pro-inflammatory monocytes, and Treg

#save the output 
write_tsv(df_sig, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_ref_conti_matrix_sig.txt")

