#!/usr/bin/env Rscript

#script to perform FindAllMarkers on the stable rectum clusters Batch 1-13

#11/5/24

#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)

#stable rectum without doublet cells 
rectum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_75_15_2000_rpca_sub_cluster_11_6_24.rds")

#use the orginial data for DEG testing 
DefaultAssay(rectum) <- "RNA"

#set idents to seurat clusters
Idents(rectum) <- rectum$seurat_clusters

#run FindAllMarkers to determine the DEGs for each cluster 
rectum.markers <- FindAllMarkers(rectum, min.pct = 0.25)
rectum.markers %>%
  group_by(cluster)

#write the output to a file 
write.table(rectum.markers, "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/helm_batch1_13_rectum_low_res_stable_fm_11_6_24.txt", sep="\t", quote=F, row.names = T)
