#!/usr/bin/env Rscript

#script to perform FindAllMarkers on the stable colon clusters Batch 1-13

#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)

#stable colon without doublet cells 
colon <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_30_2000_rpca_sub_cluster_11_4_24.rds")

#use the orginial data for DEG testing 
DefaultAssay(colon) <- "RNA"

#set idents to seurat clusters
Idents(colon) <- colon$seurat_clusters

#run FindAllMarkers to determine the DEGs for each cluster 
colon.markers <- FindAllMarkers(colon, min.pct = 0.25)
colon.markers %>%
  group_by(cluster)

#write the output to a file 
write.table(colon.markers, "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/helm_batch1_13_colon_low_res_stable_fm_11_5_24.txt", sep="\t", quote=F, row.names = T)
