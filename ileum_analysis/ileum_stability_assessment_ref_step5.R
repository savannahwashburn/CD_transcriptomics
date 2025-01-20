#!/usr/bin/env Rscript

#script to perform FindAllMarkers on the stable ileum clusters Batch 1-13

#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)


#DEG analysis on cluster

#ileum reference parameter 1
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/reference/seurat_object/helm_batch1_13_ileum_75_30_2000_ref_rpca_11_11_24.rds")

#use the original data for DEG testing 
DefaultAssay(ileum) <- "RNA"

#set idents to seurat clusters
Idents(ileum) <- ileum$seurat_clusters

#run FindAllMarkers to determine the DEGs for each cluster 
ileum.markers <- FindAllMarkers(ileum, min.pct = 0.25)
ileum.markers %>%
  group_by(cluster)

#write the output to a file 
write.table(ileum.markers, "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/reference/helm_batch1_13_ileum_ref_param_75_30_2000_fm_11_13_24.txt", sep="\t", quote=F, row.names = T)


#DEG analysis on annotated cell types

#use the original data for DEG testing 
DefaultAssay(ileum) <- "RNA"


#set idents to broad cell type 
ileum@meta.data$cell_typev2 <- ileum@meta.data$seurat_clusters
ileum$cell_typev2 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35'),
  to = c("Naive B Cell", "Enterocyte", 'Memory B Cell', 'Enterocyte', 'Stem/Paneth Cell', "Stem/Paneth Cell", 'CD4 T Cell', 'Macrophage', "Plasma Cell", 'Monocyte', 'Stem Cell', 'Plasma Cell', 'Goblet Cell', 'Enterocyte', 'Goblet Cell', 'Ribosomal Cluster', 'Regulatory T Cell', 'NK/CD8 T Cell', 'CD8 T Cell', 'Ambiguous T Cell', 'Inflammatory Monocyte', 'TA Cell', 'Cycling B Cell', 'CD4 T Cell', 'Mesenchymal Cell', 'CD8 T Cell', 'Mesenchymal Cell', 'Endothelial Cell', 'Tuft Cell', 'Mast Cell', 'Enterocyte', 'Naive B Cell', 'Doublet Cluster', 'Enteroendocrine', 'ILC3', 'Plasmablast')
)

Idents(ileum) <- ileum$cell_typev2

#run FindAllMarkers to determine the DEGs for each cluster 
ileum.markers <- FindAllMarkers(ileum, min.pct = 0.25)
ileum.markers %>%
  group_by(cluster)

#write the output to a file 
write.table(ileum.markers, "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/reference/helm_batch1_13_ileum_ref_param_75_30_2000_fm_celltypev2_11_17_24.txt", sep="\t", quote=F, row.names = T)
