#!/usr/bin/env Rscript

#Batch 1 - 13 test different cluster parameters of the "stable" ileum cells

#load in libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)

#low resolution stably assigned cells, reference, and high resolution stably assigned cells


#------low resolution stable cells-------#
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_low_res_stable_cell_10_14_24.rds")

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
      
      DefaultAssay(ileum) <- "RNA"
      # split the dataset into a list of two seurat objects (stim and CTRL)
      ileum.list <- SplitObject(ileum, split.by = "batch")
      
      # normalize and identify variable features for each dataset independently
      ileum.list <- lapply(X = ileum.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg)
      })
      
      use <- "this used "
      print(paste(use, hvg))
      
      # select features that are repeatedly variable across datasets for integration run PCA on each
      # dataset using these features
      features <- SelectIntegrationFeatures(object.list = ileum.list)
      ileum.list <- lapply(X = ileum.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      
      #perform integration
      ileum.anchors <- FindIntegrationAnchors(object.list = ileum.list, anchor.features = features, reduction = "rpca")
      
      # this command creates an 'integrated' data assay
      ileum.combined <- IntegrateData(anchorset = ileum.anchors)
      
      # specify that we will perform downstream analysis on the corrected data note that the
      # original unmodified data still resides in the 'RNA' assay
      DefaultAssay(ileum.combined) <- "integrated"
      
      # Run the standard workflow for visualization and clustering
      ileum.combined <- ScaleData(ileum.combined, verbose = FALSE)
      ileum.combined <- RunPCA(ileum.combined, npcs = dim, verbose = FALSE)
      ileum.combined <- RunUMAP(ileum.combined, reduction = "pca", dims = 1:dim)
      print(paste(use, dim))
      ileum.combined <- FindNeighbors(ileum.combined, reduction = "pca", dims = 1:dim)
      ileum.combined <- FindClusters(ileum.combined, resolution = res)
      print(paste(use, res))
      
      # Save clustering results to metadata matrix
      clustering_results <- ileum.combined$seurat_clusters
      
      clustering_results_list[[column_name]] <- clustering_results
    }
  }
}

cluster_results_df <- as.data.frame(clustering_results_list)
print('converted list to dataframe')

write.csv(cluster_results_df, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/helm_batch1_13_ileum_low_res_cluster_param_10_24_24.csv")
ileum.combined <- AddMetaData(object = ileum.combined, metadata = cluster_results_df)
saveRDS(ileum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_9_ileum_low_res_cluster_param_10_24_24.rds")


#-------reference ileum data------------#
#ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_low_res_ref_10_21_24.rds")

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
      
      DefaultAssay(ileum) <- "RNA"
      # split the dataset into a list of two seurat objects (stim and CTRL)
      ileum.list <- SplitObject(ileum, split.by = "batch")
      
      # normalize and identify variable features for each dataset independently
      ileum.list <- lapply(X = ileum.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg)
      })
      
      use <- "this used "
      print(paste(use, hvg))
      
      # select features that are repeatedly variable across datasets for integration run PCA on each
      # dataset using these features
      features <- SelectIntegrationFeatures(object.list = ileum.list)
      ileum.list <- lapply(X = ileum.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      
      #perform integration
      ileum.anchors <- FindIntegrationAnchors(object.list = ileum.list, anchor.features = features, reduction = "rpca")
      
      # this command creates an 'integrated' data assay
      ileum.combined <- IntegrateData(anchorset = ileum.anchors)
      
      # specify that we will perform downstream analysis on the corrected data note that the
      # original unmodified data still resides in the 'RNA' assay
      DefaultAssay(ileum.combined) <- "integrated"
      
      # Run the standard workflow for visualization and clustering
      ileum.combined <- ScaleData(ileum.combined, verbose = FALSE)
      ileum.combined <- RunPCA(ileum.combined, npcs = dim, verbose = FALSE)
      ileum.combined <- RunUMAP(ileum.combined, reduction = "pca", dims = 1:dim)
      print(paste(use, dim))
      ileum.combined <- FindNeighbors(ileum.combined, reduction = "pca", dims = 1:dim)
      ileum.combined <- FindClusters(ileum.combined, resolution = res)
      print(paste(use, res))
      
      # Save clustering results to metadata matrix
      clustering_results <- ileum.combined$seurat_clusters
      
      clustering_results_list[[column_name]] <- clustering_results
    }
  }
}

cluster_results_df <- as.data.frame(clustering_results_list)
print('converted list to dataframe')

write.csv(cluster_results_df, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/helm_batch1_13_ileum_ref_cluster_param_10_25_24.csv")
ileum.combined <- AddMetaData(object = ileum.combined, metadata = cluster_results_df)
saveRDS(ileum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/helm_batch1_9_ileum_ref_cluster_param_10_25_24.rds")

#--------high resolution stable cells-------#
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_13_ileum_high_res_stable_cell_10_28_24.rds")


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
      
      DefaultAssay(ileum) <- "RNA"
      # split the dataset into a list of two seurat objects (stim and CTRL)
      ileum.list <- SplitObject(ileum, split.by = "batch")
      
      # normalize and identify variable features for each dataset independently
      ileum.list <- lapply(X = ileum.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg)
      })
      
      use <- "this used "
      print(paste(use, hvg))
      
      # select features that are repeatedly variable across datasets for integration run PCA on each
      # dataset using these features
      features <- SelectIntegrationFeatures(object.list = ileum.list)
      ileum.list <- lapply(X = ileum.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      
      #perform integration
      ileum.anchors <- FindIntegrationAnchors(object.list = ileum.list, anchor.features = features, reduction = "rpca")
      
      # this command creates an 'integrated' data assay
      ileum.combined <- IntegrateData(anchorset = ileum.anchors)
      
      # specify that we will perform downstream analysis on the corrected data note that the
      # original unmodified data still resides in the 'RNA' assay
      DefaultAssay(ileum.combined) <- "integrated"
      
      # Run the standard workflow for visualization and clustering
      ileum.combined <- ScaleData(ileum.combined, verbose = FALSE)
      ileum.combined <- RunPCA(ileum.combined, npcs = dim, verbose = FALSE)
      ileum.combined <- RunUMAP(ileum.combined, reduction = "pca", dims = 1:dim)
      print(paste(use, dim))
      ileum.combined <- FindNeighbors(ileum.combined, reduction = "pca", dims = 1:dim)
      ileum.combined <- FindClusters(ileum.combined, resolution = res)
      print(paste(use, res))
      
      # Save clustering results to metadata matrix
      clustering_results <- ileum.combined$seurat_clusters
      
      clustering_results_list[[column_name]] <- clustering_results
    }
  }
}

cluster_results_df <- as.data.frame(clustering_results_list)
print('converted list to dataframe')

write.csv(cluster_results_df, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/helm_batch1_13_ileum_high_res_cluster_param_10_28_24.csv")
ileum.combined <- AddMetaData(object = ileum.combined, metadata = cluster_results_df)
saveRDS(ileum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_9_ileum_high_res_cluster_param_10_28_24.rds")
