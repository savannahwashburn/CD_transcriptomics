#!/usr/bin/env Rscript

#Batch 1 - 13 test different cluster parameters of the stably assigned colon cells

#load in libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)


#------low resolution stable cells-------#
colon <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_low_res_stable_cell_10_29_24.rds")

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
      
      DefaultAssay(colon) <- "RNA"
      # split the dataset into a list of two seurat objects (stim and CTRL)
      colon.list <- SplitObject(colon, split.by = "batch")
      
      # normalize and identify variable features for each dataset independently
      colon.list <- lapply(X = colon.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg)
      })
      
      use <- "this used "
      print(paste(use, hvg))
      
      # select features that are repeatedly variable across datasets for integration run PCA on each
      # dataset using these features
      features <- SelectIntegrationFeatures(object.list = colon.list)
      colon.list <- lapply(X = colon.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      
      #perform integration
      colon.anchors <- FindIntegrationAnchors(object.list = colon.list, anchor.features = features, reduction = "rpca")
      
      # this command creates an 'integrated' data assay
      colon.combined <- IntegrateData(anchorset = colon.anchors)
      
      # specify that we will perform downstream analysis on the corrected data note that the
      # original unmodified data still resides in the 'RNA' assay
      DefaultAssay(colon.combined) <- "integrated"
      
      # Run the standard workflow for visualization and clustering
      colon.combined <- ScaleData(colon.combined, verbose = FALSE)
      colon.combined <- RunPCA(colon.combined, npcs = dim, verbose = FALSE)
      colon.combined <- RunUMAP(colon.combined, reduction = "pca", dims = 1:dim)
      print(paste(use, dim))
      colon.combined <- FindNeighbors(colon.combined, reduction = "pca", dims = 1:dim)
      colon.combined <- FindClusters(colon.combined, resolution = res)
      print(paste(use, res))
      
      # Save clustering results to metadata matrix
      clustering_results <- colon.combined$seurat_clusters
      
      clustering_results_list[[column_name]] <- clustering_results
    }
  }
}

cluster_results_df <- as.data.frame(clustering_results_list)
print('converted list to dataframe')

write.csv(cluster_results_df, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/helm_batch1_13_colon_low_res_cluster_param_10_29_24.csv")
colon.combined <- AddMetaData(object = colon.combined, metadata = cluster_results_df)
saveRDS(colon.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_low_res_cluster_param_10_29_24.rds")

