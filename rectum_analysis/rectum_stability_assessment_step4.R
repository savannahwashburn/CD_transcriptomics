#!/usr/bin/env Rscript

#11/5/2024 - cluster rectum low res top param

#libraries
#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)

#-------load in rectum low resolution stably assigned cells--------_#
rectum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_low_res_stable_cell_11_4_24.rds")

#-------change these parameters as needed----------#

#top parameters: 
  #Resolution: 0.75, PCs: 1-15, HVG: 2000
  #Resolution: 1, PCs: 1-15, HVG: 2000
  #Resolution: 0.50, PCs: 1-15, HVG: 2000

#resolution 
res <- 0.75

#PCs
dim <- 15

#HVG
hvg <- 2000

DefaultAssay(rectum) <- "RNA"

# split the dataset into a list of two seurat objects (stim and CTRL)
rectum.list <- SplitObject(rectum, split.by = "batch")

# normalize and identify variable features for each dataset independently
rectum.list <- lapply(X = rectum.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg)
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
rectum.combined <- RunPCA(rectum.combined, npcs = dim, verbose = FALSE)
rectum.combined <- RunUMAP(rectum.combined, reduction = "pca", dims = 1:dim)
rectum.combined <- FindNeighbors(rectum.combined, reduction = "pca", dims = 1:dim)
rectum.combined <- FindClusters(rectum.combined, resolution = res)

#-----save the clustered seurat object---------#
saveRDS(rectum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_75_15_2000_rpca_11_5_24.rds")
