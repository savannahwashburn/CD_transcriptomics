#!/usr/bin/env Rscript

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(tools)

#perform clustering stability assessment on ileum samples batch 1-13 


#-----------low resolution clustering stability assessment-------------#

#use rPCA batch effect correction 
random_split <- function(){
  for (i in 1:5) {
    print(i)
    #10/21/24 - read in low res reference SO 
    ileum <- readRDS(file = '/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_low_res_ref_10_21_24.rds')
    Idents(ileum) <- ileum$orig.ident
    cells.to.sample <- round(ncol(ileum) * 0.50)
    set.seed(Sys.time())
    sub_sample <- ileum[, sample(Cells(ileum), size = cells.to.sample), seed = NULL]
    sub_other <- ileum[, !(colnames(ileum) %in% colnames(sub_sample))]
    file_name <- paste0("ileum_sub_res20_", i, ".rds")
    saveRDS(sub_sample, file = file_name)
    file_name2 <- paste0('ileum_sub2_res20_', i, ".rds")
    saveRDS(sub_other, file = file_name2)
  }
}


#function to cluster each randomly split seurat object
cluster_query <- function(ileum_files){
  for (i in ileum_files) {
    #read in the SO
    ileum <- readRDS(file = i)
    #just get the name of the file wihtout ".rds" extension 
    name = file_path_sans_ext(basename(i))
    #print the name of the file
    print(name)
    #print info about seurat object
    print(ileum)
    
    
    DefaultAssay(ileum) <- "RNA"
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
    
    #create unique file name for each SO
    file_name <- paste0(name, "_low_res_cluster.rds")
    #print the new seurat object file 
    print(file_name)
    #save the seurat object
    saveRDS(ileum.combined, file = file_name)
  }
}

#function to map the reference cluster annotations to query clusters
predict_labels <- function(clustered_files, ileum){
  for (i in clustered_files) {
    cluster <- readRDS(file = i)
    name = file_path_sans_ext(basename(i))
    cluster.anchors <- FindTransferAnchors(reference = ileum, query = cluster, dims = 1:15, reference.reduction = "pca")
    predictions <- TransferData(anchorset = cluster.anchors, refdata = ileum$seurat_clusters, dims = 1:15)
    cluster <- AddMetaData(cluster, metadata = predictions)
    file_name <- paste0(name, "_pred.rds")
    saveRDS(cluster, file = file_name)
  }
}


#-----------Stability Assessment-------------#


#set the low resolution working directory 
setwd("/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/seurat_object")

#randomly split the seurat objects 
random_split()

directory = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/seurat_object"

#list of randomly subsetted seurat objects used for clustering  
ileum_files <- list.files(directory, pattern = ".rds", full.names = TRUE)

print(ileum_files)

#cluster each random subset 
cluster_query(ileum_files)

#list of clustered files to be used for clustering 
clustered_files <- list.files(directory, pattern = "cluster.rds", full.names = TRUE)

#print(clustered_files)

#load in reference 
ileum <- readRDS(file = '/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_low_res_ref_10_21_24.rds')

#predict labels 
predict_labels(clustered_files, ileum)

#----------high resolution clustering stability assessment--------------#


#function to cluster each randomly split seurat object
cluster_query <- function(ileum_files){
  for (i in ileum_files) {
    #read in the SO
    ileum <- readRDS(file = i)
    #just get the name of the file wihtout ".rds" extension 
    name = file_path_sans_ext(basename(i))
    #print the name of the file
    print(name)
    #print info about seurat object
    print(ileum)
    
    
    DefaultAssay(ileum) <- "RNA"
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
    ileum.combined <- FindClusters(ileum.combined, resolution = 1.5)
    
    #create unique file name for each SO
    file_name <- paste0(name, "_high_res_cluster.rds")
    #print the new seurat object file 
    print(file_name)
    #save the seurat object
    saveRDS(ileum.combined, file = file_name)
  }
}

#function to map the reference cluster annotations to query clusters
predict_labels <- function(clustered_files, ileum){
  for (i in clustered_files) {
    cluster <- readRDS(file = i)
    name = file_path_sans_ext(basename(i))
    cluster.anchors <- FindTransferAnchors(reference = ileum, query = cluster, dims = 1:15, reference.reduction = "pca")
    predictions <- TransferData(anchorset = cluster.anchors, refdata = ileum$seurat_clusters, dims = 1:15)
    cluster <- AddMetaData(cluster, metadata = predictions)
    file_name <- paste0(name, "_pred.rds")
    saveRDS(cluster, file = file_name)
  }
}


#-----------Stability Assessment-------------#

#set the working directory - high resolution 
setwd("/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/seurat_object")

#randomly split the seurat objects 
#random_split() - use the same random splits for high resolution stability assessment

#set the working directory
directory_high_res = "/storage/home/swashburn30/Helmsley/Batch1_9/stability_test/3_21_24/ileum/seurat_object/high_res"

#list of randomly subsetted seurat objects used for clustering  
ileum_files <- list.files(directory_high_res, pattern = ".rds", full.names = TRUE)

print(ileum_files)

#cluster each random subset 
cluster_query(ileum_files)


#list of clustered files to be used for clustering 
clustered_files <- list.files(directory_high_res, pattern = "cluster.rds", full.names = TRUE)

print(clustered_files)

#load in high resolution reference (rPCA clustering with resolution of 1.0)
ileum <- readRDS(file = '/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_rpca_res1_10_21_24.rds')

#predict labels 
predict_labels(clustered_files, ileum)

