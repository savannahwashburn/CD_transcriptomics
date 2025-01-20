#!/usr/bin/env Rscript

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(tools)

#perform clustering stability assessment on colon samples batch 1-13 

#Low resolution stability assessment 

random_split <- function(){
  for (i in 1:5) {
    print(i)
    #read in low res reference SO 
    colon <- readRDS(file = '/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_low_res_ref_10_17_24.rds')
    Idents(colon) <- colon$orig.ident
    cells.to.sample <- round(ncol(colon) * 0.50)
    set.seed(Sys.time())
    sub_sample <- colon[, sample(Cells(colon), size = cells.to.sample), seed = NULL]
    sub_other <- colon[, !(colnames(colon) %in% colnames(sub_sample))]
    file_name <- paste0("colon_sub_res20_", i, ".rds")
    saveRDS(sub_sample, file = file_name)
    file_name2 <- paste0('colon_sub2_res20_', i, ".rds")
    saveRDS(sub_other, file = file_name2)
  }
}

cluster_query <- function(colon_files){
  for (i in colon_files) {
    #read in the SO
    colon <- readRDS(file = i)
    #just get the name of the file wihtout ".rds" extension 
    name = file_path_sans_ext(basename(i))
    #print the name of the file
    print(name)
    #print info about seurat object
    print(colon)
    
    #low resolution
    
    DefaultAssay(colon) <- "RNA"
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
    #resolution of 0.20
    colon.combined <- FindClusters(colon.combined, resolution = 0.2)
    
    #create unique file name for each SO
    file_name <- paste0(name, "_low_res_cluster.rds")
    #print the new seurat object file 
    print(file_name)
    #save the seurat object
    saveRDS(colon.combined, file = file_name)
  }
}

#function to map the reference cluster annotations to query clusters
#change this to go to the integrated assay instead of pca 
predict_labels <- function(clustered_files, colon){
  for (i in clustered_files) {
    cluster <- readRDS(file = i)
    name = file_path_sans_ext(basename(i))
    cluster.anchors <- FindTransferAnchors(reference = colon, query = cluster, dims = 1:15, reference.reduction = "pca")
    predictions <- TransferData(anchorset = cluster.anchors, refdata = colon$seurat_clusters, dims = 1:15)
    cluster <- AddMetaData(cluster, metadata = predictions)
    file_name <- paste0(name, "_pred.rds")
    saveRDS(cluster, file = file_name)
  }
}

#----------perform stability assessment---------#

#set the working directory - low resolution 
setwd("/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object")

#randomly split the seurat objects 
random_split()

#set directory path
directory = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object"

#list of randomly subsetted seurat objects used for clustering  
colon_files <- list.files(directory, pattern = ".rds", full.names = TRUE)

print(colon_files)

#cluster each random subset 
cluster_query(colon_files)

#list of clustered files to be used for clustering 
clustered_files <- list.files(directory, pattern = "cluster.rds", full.names = TRUE)

print(clustered_files)

#load in low resolution reference 
colon <- readRDS(file = '/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_low_res_ref_10_17_24.rds')

#predict labels 
predict_labels(clustered_files, colon)
