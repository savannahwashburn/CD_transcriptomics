#!/usr/bin/env Rscript

#Batch 1 - 13 test different cluster parameters of the stably assigned rectum cells

#load in libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)


#------low resolution stable cells-------#
rectum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_low_res_stable_cell_11_4_24.rds")

resolution <- c(0.25, 0.5, 0.75, 1)

dimension <- c(10, 15, 30)

feature <- c(500, 2000, 5000)

#try saving as a list
clustering_results_list <- list()

for (res in resolution) {
  for (dim in dimension) {
    for (hvg in feature) {
      print(res)
      print(dim)
      print(hvg)
      column_name <- paste("cluster", res, dim, hvg, sep = "_")
      
      DefaultAssay(rectum) <- "RNA"
      # split the dataset into a list of two seurat objects (stim and CTRL)
      rectum.list <- SplitObject(rectum, split.by = "batch")
      
      # normalize and identify variable features for each dataset independently
      rectum.list <- lapply(X = rectum.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg)
      })
      
      use <- "this used "
      print(paste(use, hvg))
      
      # select features that are repeatedly variable across datasets for integration run PCA on each
      # dataset using these features
      features <- SelectIntegrationFeatures(object.list = rectum.list)
      rectum.list <- lapply(X = rectum.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      
      #perform integration
      rectum.anchors <- FindIntegrationAnchors(object.list = rectum.list, anchor.features = features, reduction = "rpca")
      
      # this command creates an 'integrated' data assay
      rectum.combined <- IntegrateData(anchorset = rectum.anchors)
      
      # specify that we will perform downstream analysis on the corrected data note that the
      # original unmodified data still resides in the 'RNA' assay
      DefaultAssay(rectum.combined) <- "integrated"
      
      # Run the standard workflow for visualization and clustering
      rectum.combined <- ScaleData(rectum.combined, verbose = FALSE)
      rectum.combined <- RunPCA(rectum.combined, npcs = dim, verbose = FALSE)
      rectum.combined <- RunUMAP(rectum.combined, reduction = "pca", dims = 1:dim)
      print(paste(use, dim))
      rectum.combined <- FindNeighbors(rectum.combined, reduction = "pca", dims = 1:dim)
      rectum.combined <- FindClusters(rectum.combined, resolution = res)
      print(paste(use, res))
      
      # Save clustering results to metadata matrix
      clustering_results <- rectum.combined$seurat_clusters
      
      clustering_results_list[[column_name]] <- clustering_results
    }
  }
}

cluster_results_df <- as.data.frame(clustering_results_list)
print('converted list to dataframe')

write.csv(cluster_results_df, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/helm_batch1_13_rectum_low_res_cluster_param_11_4_24.csv")
rectum.combined <- AddMetaData(object = rectum.combined, metadata = cluster_results_df)
saveRDS(rectum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_low_res_cluster_param_11_4_24.rds")

