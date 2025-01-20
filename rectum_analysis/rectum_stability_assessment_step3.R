#!/usr/bin/env Rscript

#script to perform ARI on low resolution rectum stably assigned cells (batch1-13)

#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Dune)

#read in dataframe 
df <- read.csv("/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/helm_batch1_13_rectum_low_res_cluster_param_11_4_24.csv")

#remove "X"
df <- df[, -1]

#get column names 
columns <- colnames(df)

merger <- Dune(clusMat = df %>% select(columns), parallel = TRUE, metric = "ARI")

#save the output
df_merge <- merger$currentMat

write.csv(df_merge, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/df_ri_helm_batch1_13_rectum_low_res_dune_merge_11_4_24.csv", row.names = T)
