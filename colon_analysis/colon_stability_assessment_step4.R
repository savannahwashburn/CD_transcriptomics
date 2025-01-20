#!/usr/bin/env Rscript

#libraries
#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)

#-------load in colon low resolution stably assigned cells--------_#
colon <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_30_2000_sub_11_4_24.rds")

#-------change these parameters as needed----------#

#top parameters assessed: 
  #resolution: 0.25, PCs: 1-30, HVG: 2000
  #Resolution: 0.50, PCs: 1-30, HVGs: 2000
  #Resolution: 0.25, PCs: 1-15, HVGs: 5000

#change these parameters as needed 

#resolution 
res <- 0.25

#PCs
dim <- 30

#HVG
hvg <- 2000

DefaultAssay(colon) <- "RNA"

# split the dataset into a list of two seurat objects (stim and CTRL)
colon.list <- SplitObject(colon, split.by = "batch")

# normalize and identify variable features for each dataset independently
colon.list <- lapply(X = colon.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg)
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
colon.combined <- RunPCA(colon.combined, npcs = dim, verbose = FALSE)
colon.combined <- RunUMAP(colon.combined, reduction = "pca", dims = 1:dim)
colon.combined <- FindNeighbors(colon.combined, reduction = "pca", dims = 1:dim)
colon.combined <- FindClusters(colon.combined, resolution = res)

#-----save the clustered seurat object---------#
saveRDS(colon.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_30_2000_rpca_11_2_24.rds")
