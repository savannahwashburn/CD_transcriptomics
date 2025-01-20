#!/usr/bin/env Rscript

#script to perform FindAllMarkers on the stable ileum clusters Batch 1-13

#11/13/2024 - perform find markers on high res stable cluster parameter 1

#11/17/2024 - find markers on broad cell type annotation 

#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)


#DEG analysis on cluster

#ileum high res stable parameter 1
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_13_ileum_75_15_2000_rpca_11_11_24.rds")

#use the orginial data for DEG testing 
DefaultAssay(ileum) <- "RNA"

#set idents to seurat clusters
Idents(ileum) <- ileum$seurat_clusters


#run FindAllMarkers to determine the DEGs for each cluster 
ileum.markers <- FindAllMarkers(ileum, min.pct = 0.25)
ileum.markers %>%
  group_by(cluster)

#write the output to a file 
write.table(ileum.markers, "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/helm_batch1_13_ileum_high_res_param_75_15_2000_stable_fm_11_13_24.txt", sep="\t", quote=F, row.names = T)

#DEG analysis on cluster annotation

#ileum high res stable parameter 1
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_13_ileum_75_15_2000_rpca_11_11_24.rds")

#use the orginial data for DEG testing 
DefaultAssay(ileum) <- "RNA"

#set idents to cell type v2
ileum@meta.data$cell_typev2 <- ileum@meta.data$seurat_clusters
ileum$cell_typev2 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36'),
  to = c("Naive B Cell", "Memory B Cell", "Enterocyte", "Stem/Paneth Cell", "CD4 T Cell", "Enterocyte", "Stem/Paneth Cell", 'CD8 T Cell', "Plasma Cell", "Enterocyte", "Monocyte", "Plasma Cell", "Enterocyte", "Goblet Cell", "Stem Cell", "Regulatory T Cell", "Ribosomal Cluster", "Goblet Cell", "Inflammatory Macrophage", "Macrophage", "Ambiguous T Cell", "Inflammatory Monocyte", "TA Cell", "NK/T Cell", "Mesenchymal Cell", "Germinal Center B Cell", "Endothelial Cell", "Mesenchymal Cell", "Mast Cell", "Tuft Cell", "Doublet Cluster", "Plasma Cell", "Naive B Cell", "ILC3", "Enteroendocrine", "Plasmablast", "Cycling NK/T Cell")
)

Idents(ileum) <- ileum$cell_typev2

#run FindAllMarkers to determine the DEGs for each cluster 
ileum.markers <- FindAllMarkers(ileum, min.pct = 0.25)
ileum.markers %>%
  group_by(cluster)

#write the output to a file 
write.table(ileum.markers, "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/helm_batch1_13_ileum_high_res_param_75_15_2000_stable_fm_cell_typev2_11_17_24.txt", sep="\t", quote=F, row.names = T)

