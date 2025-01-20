#!/usr/bin/env Rscript

#preliminary clustering of ileum data after QC; test different clustering and batch effect correction methods

#libraries
#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)x


#load in seurat object
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_10_17_24.rds")

#without batch effect correction 
ileum <- NormalizeData(ileum)
ileum <- FindVariableFeatures(ileum, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ileum)
ileum <- ScaleData(ileum, features = all.genes)
ileum <- RunPCA(ileum, features = VariableFeatures(object = ileum))
ileum <- FindNeighbors(ileum, dims = 1:15)
ileum <- FindClusters(ileum, resolution = 0.2)
ileum <- RunUMAP(ileum, dims = 1:15)

saveRDS(ileum, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_clustered_10_17_24.rds")

#----------------CCA Batch Effect Correction--------------------#

#split the dataset into a list of 6 seurat objects into batch
ileum.list <- SplitObject(ileum, split.by = "batch")

# normalize and identify variable features for each dataset independently
ileum.list <- lapply(X = ileum.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ileum.list)

#perform integration
ileum.anchors <- FindIntegrationAnchors(object.list = ileum.list, anchor.features = features)

# this command creates an 'integrated' data assay
ileum.combined <- IntegrateData(anchorset = ileum.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(ileum.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
ileum.combined <- ScaleData(ileum.combined, verbose = FALSE)
ileum.combined <- RunPCA(ileum.combined, npcs = 15, verbose = FALSE)
ileum.combined <- RunUMAP(ileum.combined, reduction = "pca", dims = 1:15)
ileum.combined <- FindNeighbors(ileum.combined, reduction = "pca", dims = 1:15)
ileum.combined <- FindClusters(ileum.combined, resolution = 0.2)

#save CCA seurat object 
saveRDS(ileum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_cca_10_17_24.rds")

#----------------rPCA Batch Effect Correction-----------------------#

# split the dataset into a list of two seurat objects (stim and CTRL)
ileum.list <- SplitObject(ileum, split.by = "batch")

# normalize and identify variable features for each dataset independently
ileum.list <- lapply(X = ileum.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

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
ileum.combined <- RunPCA(ileum.combined, npcs = 15, verbose = FALSE)
ileum.combined <- RunUMAP(ileum.combined, reduction = "pca", dims = 1:15)
ileum.combined <- FindNeighbors(ileum.combined, reduction = "pca", dims = 1:15)
ileum.combined <- FindClusters(ileum.combined, resolution = 0.2)

#save rPCA seurat object 
#saveRDS(ileum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_rpca_10_17_24.rds")

#---------------Harmony Batch Effect Correction------------------#
#correct for batch 
ileum_h <- ileum %>%
  RunHarmony("batch")

harmony_embeddings <- Embeddings(ileum_h, 'harmony')

ileum_h <- ileum_h %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.2) %>%
  identity()

#save the Harmony seurat object 
saveRDS(ileum_h, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_harmony_10_17_24.rds")
