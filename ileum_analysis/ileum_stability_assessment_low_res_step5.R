#!/usr/bin/env Rscript

#script to perform FindAllMarkers on the stable ileum clusters Batch 1-13

#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)

#perform DEG comparison per cluster

#low res 
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_50_15_2000_rpca_10_29_24.rds")

#use the orginial data for DEG testing 
DefaultAssay(ileum) <- "RNA"

#set idents to seurat clusters
Idents(ileum) <- ileum$seurat_clusters


#run FindAllMarkers to determine the DEGs for each cluster 
ileum.markers <- FindAllMarkers(ileum, min.pct = 0.25)
ileum.markers %>%
  group_by(cluster)

#write the output to a file 

write.table(ileum.markers, "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/helm_batch1_13_ileum_low_stable_50_15_2000_fm_11_14_24.txt", sep="\t", quote=F, row.names = T)


#perfom DEG comparison per cluster

#low res
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_50_15_2000_rpca_10_29_24.rds")

#use the orginial data for DEG testing 
DefaultAssay(ileum) <- "RNA"

#set idents to cell type v2
ileum@meta.data$cell_typev2 <- ileum@meta.data$seurat_clusters
ileum$cell_typev2 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24'),
  to = c("Naive B Cell", 'Stem/Paneth Cell', 'Enterocyte', 'CD4 T Cell', 'Enterocyte', 'Macrophage', "CD8 T Cell", "Plasma Cell", "Enterocyte", "Plasma Cell", "Monocyte", "Goblet Cell", "Goblet Cell", "Ambiguous T Cell", "Mesenchymal Cell", "NK/T Cell", "Inflammatory Monocyte", "TA Cell", "Memory B Cell", 'Endothelial Cell', "Mast Cell", "Tuft Cell", "Doublet Cluster", "ILC3", "Plasmablast")
)

#set idents
Idents(ileum) <- ileum$cell_typev2


#run FindAllMarkers to determine the DEGs for each cluster 
ileum.markers <- FindAllMarkers(ileum, min.pct = 0.25)
ileum.markers %>%
  group_by(cluster)

write.table(ileum.markers, "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/helm_batch1_13_ileum_low_stable_50_15_2000_fm_cell_typev2_11_17_24.txt", sep="\t", quote=F, row.names = T)
