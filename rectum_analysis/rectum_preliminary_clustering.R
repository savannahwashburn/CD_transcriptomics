#!/usr/bin/env Rscript

#preliminary clustering for rectum batch1-13

#libraries
#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)

#load in seurat object
rectum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_10_17_24.rds")

#without batch effect correction 
rectum <- NormalizeData(rectum)
rectum <- FindVariableFeatures(rectum, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rectum)
rectum <- ScaleData(rectum, features = all.genes)
rectum <- RunPCA(rectum, features = VariableFeatures(object = rectum))
rectum <- FindNeighbors(rectum, dims = 1:15)
rectum <- FindClusters(rectum, resolution = 0.2)
rectum <- RunUMAP(rectum, dims = 1:15)

#save the rectum seurat object 
saveRDS(rectum, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_clustered_10_17_24.rds")

#----------------CCA Batch Effect Correction--------------------#

# split the dataset into a list of 6 seurat objects into batch
rectum.list <- SplitObject(rectum, split.by = "batch")

# normalize and identify variable features for each dataset independently
rectum.list <- lapply(X = rectum.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = rectum.list)

#perform integration
rectum.anchors <- FindIntegrationAnchors(object.list = rectum.list, anchor.features = features)

# this command creates an 'integrated' data assay
rectum.combined <- IntegrateData(anchorset = rectum.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(rectum.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
rectum.combined <- ScaleData(rectum.combined, verbose = FALSE)
rectum.combined <- RunPCA(rectum.combined, npcs = 15, verbose = FALSE)
rectum.combined <- RunUMAP(rectum.combined, reduction = "pca", dims = 1:15)
rectum.combined <- FindNeighbors(rectum.combined, reduction = "pca", dims = 1:15)
rectum.combined <- FindClusters(rectum.combined, resolution = 0.2)

#save CCA seurat object 
saveRDS(rectum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_cca_10_17_24.rds")

#----------------rPCA Batch Effect Correction-----------------------#

# split the dataset into a list of two seurat objects (stim and CTRL)
rectum.list <- SplitObject(rectum, split.by = "batch")

# normalize and identify variable features for each dataset independently
rectum.list <- lapply(X = rectum.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

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
rectum.combined <- RunPCA(rectum.combined, npcs = 15, verbose = FALSE)
rectum.combined <- RunUMAP(rectum.combined, reduction = "pca", dims = 1:15)
rectum.combined <- FindNeighbors(rectum.combined, reduction = "pca", dims = 1:15)
rectum.combined <- FindClusters(rectum.combined, resolution = 0.2)

#save rPCA seurat object 
saveRDS(rectum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_rpca_10_17_24.rds")


#---------------Harmony Batch Effect Correction------------------#

#use clustered results
rectum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_clustered_10_17_24.rds")

#correct for batch 
rectum_h <- rectum %>%
  RunHarmony("batch")

harmony_embeddings <- Embeddings(rectum_h, 'harmony')

rectum_h <- rectum_h %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.2) %>%
  identity()

#save the Harmony seurat object 
saveRDS(rectum_h, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_harmony_10_17_24.rds")
