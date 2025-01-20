#!/usr/bin/env Rscript

#11/11/2024 - cluster the top reference parameters 

#libraries
#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)

#-------load in ileum reference--------_#
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_low_res_ref_10_21_24.rds")

#-------change these parameters as needed----------#

#top clustering results from ARI score: 
  #resolution: 0.75, PCs: 1-30, HVG: 2000
  #resolution: 0.75, PCs: 1-15, HVG: 500
  #resolution: 1.0, PCs: 1-30, HVG: 5000

#resolution 
res <- 0.75

#PCs
dim <- 30

#HVG
hvg <- 2000

DefaultAssay(ileum) <- "RNA"

# split the dataset into a list of two seurat objects (stim and CTRL)
ileum.list <- SplitObject(ileum, split.by = "batch")

# normalize and identify variable features for each dataset independently
ileum.list <- lapply(X = ileum.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg)
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
ileum.combined <- RunPCA(ileum.combined, npcs = dim, verbose = FALSE)
ileum.combined <- RunUMAP(ileum.combined, reduction = "pca", dims = 1:dim)
ileum.combined <- FindNeighbors(ileum.combined, reduction = "pca", dims = 1:dim)
ileum.combined <- FindClusters(ileum.combined, resolution = res)

#-----save the clustered seurat object---------#
saveRDS(ileum.combined, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/reference/seurat_object/helm_batch1_13_ileum_75_30_2000_ref_rpca_11_11_24.rds")
