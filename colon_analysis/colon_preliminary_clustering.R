#!/usr/bin/env Rscript

#preliminary clustering for the colon data 

#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)


#load in seurat object
colon <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_10_17_24.rds")

colon <- NormalizeData(colon)
colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(colon)
colon <- ScaleData(colon, features = all.genes)
colon <- RunPCA(colon, features = VariableFeatures(object = colon))
colon <- FindNeighbors(colon, dims = 1:15)
colon <- FindClusters(colon, resolution = 0.2)
colon <- RunUMAP(colon, dims = 1:15)

#save the colon seurat object 
saveRDS(colon, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_clustered_10_17_24.rds")

#----------------CCA Batch Effect Correction--------------------#

# split the dataset into a list of 6 seurat objects into batch
colon.list <- SplitObject(colon, split.by = "batch")

# normalize and identify variable features for each dataset independently
colon.list <- lapply(X = colon.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = colon.list)

#perform integration
colon.anchors <- FindIntegrationAnchors(object.list = colon.list, anchor.features = features)

# this command creates an 'integrated' data assay
colon.combined <- IntegrateData(anchorset = colon.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(colon.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
colon.combined <- ScaleData(colon.combined, verbose = FALSE)
colon.combined <- RunPCA(colon.combined, npcs = 15, verbose = FALSE)
colon.combined <- RunUMAP(colon.combined, reduction = "pca", dims = 1:15)
colon.combined <- FindNeighbors(colon.combined, reduction = "pca", dims = 1:15)
colon.combined <- FindClusters(colon.combined, resolution = 0.2)

#save CCA seurat object 
saveRDS(colon.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_cca_10_17_24.rds")


#----------------rPCA Batch Effect Correction-----------------------#

# split the dataset into a list of two seurat objects (stim and CTRL)
colon.list <- SplitObject(colon, split.by = "batch")

# normalize and identify variable features for each dataset independently
colon.list <- lapply(X = colon.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

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
colon.combined <- RunPCA(colon.combined, npcs = 15, verbose = FALSE)
colon.combined <- RunUMAP(colon.combined, reduction = "pca", dims = 1:15)
colon.combined <- FindNeighbors(colon.combined, reduction = "pca", dims = 1:15)
colon.combined <- FindClusters(colon.combined, resolution = 0.2)

#save rPCA seurat object 
saveRDS(colon.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_rpca_10_17_24.rds")

#---------------Harmony Batch Effect Correction------------------#

#use clustered data 
colon <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_clustered_10_17_24.rds")

#correct for batch 
colon_h <- colon %>%
  RunHarmony("batch")

harmony_embeddings <- Embeddings(colon_h, 'harmony')

colon_h <- colon_h %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.2) %>%
  identity()

#save the Harmony seurat object 
saveRDS(colon_h, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_harmony_10_17_24.rds")
