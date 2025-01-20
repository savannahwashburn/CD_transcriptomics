#!/usr/bin/env Rscript

#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Dune)

#script to compute ARI with Dune

#---------low resolution stably assigned cells clustering parameter results-----------#

#read in dataframe 
df <- read.csv("/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/helm_batch1_13_ileum_low_res_cluster_param_10_24_24.csv")

#remove "X"
df <- df[, -1]

#get column names 
columns <- colnames(df)

merger <- Dune(clusMat = df %>% select(columns), parallel = TRUE, metric = "ARI")

#save the output
df_merge <- merger$currentMat

write.csv(df_merge, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/df_ri_helm_batch1_13_ileum_low_res_dune_merge_10_25_24.csv", row.names = T)

#----------reference clustering parameter results----------#

#read in dataframe 
df <- read.csv("/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/helm_batch1_13_ileum_ref_cluster_param_10_25_24.csv")

#remove "X"
df <- df[, -1]

#get column names 
columns <- colnames(df)

merger <- Dune(clusMat = df %>% select(columns), parallel = TRUE, metric = "ARI")

#save the output
df_merge <- merger$currentMat

write.csv(df_merge, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/df_ri_helm_batch1_13_ileum_ref_dune_merge_10_25_24.csv", row.names = T)

#---------high resolution stably assigned cells clustering parameter results-----------#

df <- read.csv("/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/helm_batch1_13_ileum_high_res_cluster_param_10_28_24.csv")

#remove "X"
df <- df[, -1]

#get column names 
columns <- colnames(df)

merger <- Dune(clusMat = df %>% select(columns), parallel = TRUE, metric = "ARI")

#save the output
df_merge <- merger$currentMat

write.csv(df_merge, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/df_ri_helm_batch1_13_ileum_high_res_dune_merge_10_29_24.csv", row.names = T)
