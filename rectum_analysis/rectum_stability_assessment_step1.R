#!/usr/bin/env Rscript

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(tools)

#perform clustering stability assessment on rectum samples batch 1-13 

#low resolution stability assessment 
random_split <- function(){
  for (i in 1:5) {
    print(i)
    #read in low res reference SO 
    rectum <- readRDS(file = '/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_low_res_ref_10_21_24.rds')
    Idents(rectum) <- rectum$orig.ident
    cells.to.sample <- round(ncol(rectum) * 0.50)
    set.seed(Sys.time())
    sub_sample <- rectum[, sample(Cells(rectum), size = cells.to.sample), seed = NULL]
    sub_other <- rectum[, !(colnames(rectum) %in% colnames(sub_sample))]
    file_name <- paste0("rectum_sub_res50_", i, ".rds")
    saveRDS(sub_sample, file = file_name)
    file_name2 <- paste0('rectum_sub2_res50_', i, ".rds")
    saveRDS(sub_other, file = file_name2)
  }
}

#function to cluster each randomly split seurat object
cluster_query <- function(rectum_files){
  for (i in rectum_files) {
    #read in the SO
    rectum <- readRDS(file = i)
    #just get the name of the file wihtout ".rds" extension 
    name = file_path_sans_ext(basename(i))
    #print the name of the file
    print(name)
    #print info about seurat object
    print(rectum)
    
    # low resolution 
    
    DefaultAssay(rectum) <- "RNA"
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
    rectum.combined <- FindClusters(rectum.combined, resolution = 0.50)
    
    #create unique file name for each SO
    #change file name to high res 
    file_name <- paste0(name, "_low_res_cluster.rds")
    #print the new seurat object file 
    print(file_name)
    #save the seurat object
    saveRDS(rectum.combined, file = file_name)
  }
}

#function to map the reference cluster annotations to query clusters
#change this to go to the integrated assay instead of pca 
predict_labels <- function(clustered_files, rectum){
  for (i in clustered_files) {
    cluster <- readRDS(file = i)
    name = file_path_sans_ext(basename(i))
    cluster.anchors <- FindTransferAnchors(reference = rectum, query = cluster, dims = 1:15, reference.reduction = "pca")
    predictions <- TransferData(anchorset = cluster.anchors, refdata = rectum$seurat_clusters, dims = 1:15)
    cluster <- AddMetaData(cluster, metadata = predictions)
    file_name <- paste0(name, "_pred.rds")
    saveRDS(cluster, file = file_name)
  }
}

#----------perform stability assessment---------#

#set the working directory
setwd("/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object")

#randomly split the seurat objects 
random_split()

#low resolution directory
directory = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object"

#list of randomly subsetted seurat objects used for clustering  
rectum_files <- list.files(directory, pattern = ".rds", full.names = TRUE)

print(rectum_files)

#cluster each random subset 
cluster_query(rectum_files)

#list of clustered files to be used for clustering 
clustered_files <- list.files(directory, pattern = "cluster.rds", full.names = TRUE)

print(clustered_files)

#load in reference 
rectum <- readRDS(file = '/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_low_res_ref_10_21_24.rds')

#predict labels 
predict_labels(clustered_files, rectum)
