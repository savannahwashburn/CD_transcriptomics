#Template for Clustering Stability Assessment

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(tools)
library(pheatmap)
library(dbplyr)
library(Dune)

########Step 1: Perform Clustering Stability Assessment############

#can adjust clustering parameters in function

#use rPCA batch effect correction 
random_split <- function(){
  for (i in 1:5) {
    print(i)
    #10/21/24 - read in low res reference SO 
    seurat_object <- readRDS(file = '/path/to/seurat_object.rds')
    Idents(seurat_object) <- seurat_object$orig.ident
    cells.to.sample <- round(ncol(seurat_object) * 0.50)
    set.seed(Sys.time())
    sub_sample <- seurat_object[, sample(Cells(seurat_object), size = cells.to.sample), seed = NULL]
    sub_other <- seurat_object[, !(colnames(seurat_object) %in% colnames(sub_sample))]
    file_name <- paste0("seurat_object_query_", i, ".rds")
    saveRDS(sub_sample, file = file_name)
    file_name2 <- paste0('seurat_object_query2', i, ".rds")
    saveRDS(sub_other, file = file_name2)
  }
}


#function to cluster each randomly split seurat object
cluster_query <- function(so_files){
  for (i in so_files) {
    #read in the SO
    seurat_object <- readRDS(file = i)
    #just get the name of the file wihtout ".rds" extension 
    name = file_path_sans_ext(basename(i))
    #print the name of the file
    print(name)
    #print info about seurat object
    print(seurat_object)
    
    #cluster with Seurat rPCA clustering pipeline - can substitute with other batch effect correction methods 
    DefaultAssay(seurat_object) <- "RNA"
    # split the dataset into a list of two seurat objects (stim and CTRL)
    seurat_object.list <- SplitObject(seurat_object, split.by = "batch")
    
    # normalize and identify variable features for each dataset independently
    seurat_object.list <- lapply(X = seurat_object.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    
    # select features that are repeatedly variable across datasets for integration run PCA on each
    # dataset using these features
    features <- SelectIntegrationFeatures(object.list = seurat_object.list)
    seurat_object.list <- lapply(X = seurat_object.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
    
    #perform integration
    seurat_object.anchors <- FindIntegrationAnchors(object.list = seurat_object.list, anchor.features = features, reduction = "rpca")
    
    # this command creates an 'integrated' data assay
    seurat_object.combined <- IntegrateData(anchorset = seurat_object.anchors)
    
    # specify that we will perform downstream analysis on the corrected data note that the
    # original unmodified data still resides in the 'RNA' assay
    DefaultAssay(seurat_object.combined) <- "integrated"
    
    # Run the standard workflow for visualization and clustering
    seurat_object.combined <- ScaleData(seurat_object.combined, verbose = FALSE)
    seurat_object.combined <- RunPCA(seurat_object.combined, npcs = 15, verbose = FALSE)
    seurat_object.combined <- RunUMAP(seurat_objectcombined, reduction = "pca", dims = 1:15)
    seurat_object.combined <- FindNeighbors(seurat_object.combined, reduction = "pca", dims = 1:15)
    seurat_object.combined <- FindClusters(seurat_object.combined, resolution = 0.2)
    
    #create unique file name for each SO
    
    file_name <- paste0(name, "_low_res_cluster.rds")
    #print the new seurat object file 
    print(file_name)
    #save the seurat object
    saveRDS(seurat_object.combined, file = file_name)
  }
}

#function to map the reference cluster annotations to query clusters
predict_labels <- function(clustered_files, seurat_object){
  for (i in clustered_files) {
    cluster <- readRDS(file = i)
    name = file_path_sans_ext(basename(i))
    cluster.anchors <- FindTransferAnchors(reference = seurat_object, query = cluster, dims = 1:15, reference.reduction = "pca")
    predictions <- TransferData(anchorset = cluster.anchors, refdata = seurat_object$seurat_clusters, dims = 1:15)
    cluster <- AddMetaData(cluster, metadata = predictions)
    file_name <- paste0(name, "_pred.rds")
    saveRDS(cluster, file = file_name)
  }
}


#-----------Stability Assessment-------------#


#set the working directory
setwd("/path/to/working_directory")

#randomly split the seurat objects 
random_split()

directory = "/path/to/directory/with/query/seurat_object"

#list of randomly subsetted seurat objects used for clustering  
so_files <- list.files(directory, pattern = ".rds", full.names = TRUE)

print(so_files)

#cluster each random subset 
cluster_query(so_files)

#list of clustered files to be used for clustering 
clustered_files <- list.files(directory, pattern = "cluster.rds", full.names = TRUE)

#print(clustered_files)

#load in reference 
seurat_object <- readRDS(file = '/path/to/reference/seurat_object.rds')

#predict labels 
predict_labels(clustered_files, seurat_object)


##########Step 2:Compare Query to Reference Cluster Annotations########

#example for one query; repeat for each query

cluster1 <- readRDS(file = "/path/to/query1_seurat_object.rds")

#make sure the predicted.id (from TransferData function) is in order - change based on number of clusters in reference
cluster1$predicted.id <- factor(cluster1$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster1, reduction = 'umap', label = TRUE) + DimPlot(cluster1, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table1 <- table(cluster1$predicted.id, cluster1$seurat_clusters)
table1 <- as.data.frame.matrix(table1)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table1_prop <- as.data.frame.matrix(apply(table1, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster1) <- cluster1$seurat_clusters
sub0 <- WhichCells(cluster1, idents = c('0'))

Idents(cluster1) <- cluster1$predicted.id
pred0 <- WhichCells(cluster1, idents = '0')

int0 <- intersect(sub0, pred0) 

#cluster 1
Idents(cluster1) <- cluster1$seurat_clusters
sub1 <- WhichCells(cluster1, idents = c('1'))

Idents(cluster1) <- cluster1$predicted.id
pred1 <- WhichCells(cluster1, idents = '1')

int1 <- intersect(sub1, pred1) 

#cluster 2
Idents(cluster1) <- cluster1$seurat_clusters
sub2 <- WhichCells(cluster1, idents = c('2', '9'))

Idents(cluster1) <- cluster1$predicted.id
pred2 <- WhichCells(cluster1, idents = '2')

int2 <- intersect(sub2, pred2) 

#cluster 3
Idents(cluster1) <- cluster1$seurat_clusters
sub3 <- WhichCells(cluster1, idents = c('3'))

Idents(cluster1) <- cluster1$predicted.id
pred3 <- WhichCells(cluster1, idents = '3')

int3 <- intersect(sub3, pred3) 

#cluster 4
Idents(cluster1) <- cluster1$seurat_clusters
sub4 <- WhichCells(cluster1, idents = c('4'))

Idents(cluster1) <- cluster1$predicted.id
pred4 <- WhichCells(cluster1, idents = '4')

int4 <- intersect(sub4, pred4) 

#cluster 5
Idents(cluster1) <- cluster1$seurat_clusters
sub5 <- WhichCells(cluster1, idents = c('5'))

Idents(cluster1) <- cluster1$predicted.id
pred5 <- WhichCells(cluster1, idents = '5')

int5 <- intersect(sub5, pred5) 

#cluster 6
Idents(cluster1) <- cluster1$seurat_clusters
sub6 <- WhichCells(cluster1, idents = c('6'))

Idents(cluster1) <- cluster1$predicted.id
pred6 <- WhichCells(cluster1, idents = '6')

int6 <- intersect(sub6, pred6) 

#cluster 7
Idents(cluster1) <- cluster1$seurat_clusters
sub7 <- WhichCells(cluster1, idents = c('1', '7'))

Idents(cluster1) <- cluster1$predicted.id
pred7<- WhichCells(cluster1, idents = '7')

int7 <- intersect(sub7, pred7) 

#cluster 8 
Idents(cluster1) <- cluster1$seurat_clusters
sub8 <- WhichCells(cluster1, idents = c('8'))

Idents(cluster1) <- cluster1$predicted.id
pred8 <- WhichCells(cluster1, idents = '8')

int8 <- intersect(sub8, pred8) 

#cluster 9
Idents(cluster1) <- cluster1$seurat_clusters
sub9 <- WhichCells(cluster1, idents = c('10'))

Idents(cluster1) <- cluster1$predicted.id
pred9 <- WhichCells(cluster1, idents = '9')

int9 <- intersect(sub9, pred9) 

#cluster 10
Idents(cluster1) <- cluster1$seurat_clusters
sub10 <- WhichCells(cluster1, idents = c('11'))

Idents(cluster1) <- cluster1$predicted.id
pred10 <- WhichCells(cluster1, idents = '10')

int10 <- intersect(sub10, pred10) 

#cluster 11
Idents(cluster1) <- cluster1$seurat_clusters
sub11 <- WhichCells(cluster1, idents = c('12'))

Idents(cluster1) <- cluster1$predicted.id
pred11 <- WhichCells(cluster1, idents = '11')

int11 <- intersect(sub11, pred11) 

#cluster 12
Idents(cluster1) <- cluster1$seurat_clusters
sub12 <- WhichCells(cluster1, idents = c('14'))

Idents(cluster1) <- cluster1$predicted.id
pred12 <- WhichCells(cluster1, idents = '12')

int12 <- intersect(sub12, pred12) 

#cluster 13
Idents(cluster1) <- cluster1$seurat_clusters
sub13 <- WhichCells(cluster1, idents = c('13'))

Idents(cluster1) <- cluster1$predicted.id
pred13 <- WhichCells(cluster1, idents = '13')

int13 <- intersect(sub13, pred13) 

#cluster 14
Idents(cluster1) <- cluster1$seurat_clusters
sub14 <- WhichCells(cluster1, idents = c('15'))

Idents(cluster1) <- cluster1$predicted.id
pred14 <- WhichCells(cluster1, idents = '14')

int14 <- intersect(sub14, pred14) 

#cluster 15
Idents(cluster1) <- cluster1$seurat_clusters
sub15 <- WhichCells(cluster1, idents = c('10'))

Idents(cluster1) <- cluster1$predicted.id
pred15 <- WhichCells(cluster1, idents = '15')

int15 <- intersect(sub15, pred15) 

#total 
#pre:43892
#42348
total1 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total1,file = '/path/to/save/barcodes.csv', row.names=F)

#repeat this process for each query 

#-----load in the stable barcodes----------#
#load in barcodes for each query 
total1 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)
total2 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)
total3 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)
total4 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)
total5 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)
total6 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)
total7 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)
total8 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)
total9 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)
total10 <- scan('/path/to/save/barcodes.csv', what = "", sep = ",", skip = 1)

#merge the barcodes into a list 
total_merge <- c(total1, total2, total3, total4, total5, total6, total7, total8, total9, total10)

#counts the number of times each barcode is present
value_counts <- table(total_merge)
table(value_counts)

#define a threshold - more times a barcode is present, more stable it is
threshold <- 4

#filter values that meet your threshold 
filtered_values <- names(value_counts[value_counts >= threshold])  

#1: 87438 - keeps ~99.9% of data
#2: 86601 - keeps ~98.65% of data 
#3: 85493 - keeps ~97.38% of data
#4: 83702 - keeps ~95.34% of data ----select threshold of 4 (biggest jump after threshold of 4)
#5: 79775 - keeps ~90.8% of data

#save the filtered barcodes 
write.csv(filtered_values, file = '/path/to/stably_assigned_cells.csv', row.names = F)

###########Step 3: assess clustering similarity ############

#use Dune to calculate adjusted rand index (ARI)

#first cluster with variety of parameters

#load in stably assigned cell (or reference) seurat object 
seurat_object <- readRDS(file = "/path/to/stably_assigned_cell_seurat_object.rds")


resolution <- c(0.25, 0.5, 0.75, 1)

dimension <- c(10, 15, 30)

feature <- c(500, 2000, 5000)

#save clustering results as a list
clustering_results_list <- list()

for (res in resolution) {
  for (dim in dimension) {
    for (hvg in feature) {
      print(res)
      print(dim)
      print(hvg)
      column_name <- paste("cluster", res, dim, hvg, sep = "_")
      
      #perform typical seurat rPCA batch effect clustering pipeline
      
      DefaultAssay(seurat_object) <- "RNA"
      # split the dataset into a list of two seurat objects (stim and CTRL)
      seurat_object.list <- SplitObject(seurat_object, split.by = "batch")
      
      # normalize and identify variable features for each dataset independently
      seurat_object.list <- lapply(X = seurat_object.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg)
      })
      
      use <- "this used "
      print(paste(use, hvg))
      
      # select features that are repeatedly variable across datasets for integration run PCA on each
      # dataset using these features
      features <- SelectIntegrationFeatures(object.list = seurat_object.list)
      ileum.list <- lapply(X = seurat_object.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      
      #perform integration
      seurat_object.anchors <- FindIntegrationAnchors(object.list = seurat_object.list, anchor.features = features, reduction = "rpca")
      
      # this command creates an 'integrated' data assay
      seurat_object.combined <- IntegrateData(anchorset = seurat_object.anchors)
      
      # specify that we will perform downstream analysis on the corrected data note that the
      # original unmodified data still resides in the 'RNA' assay
      DefaultAssay(seurat_object.combined) <- "integrated"
      
      # Run the standard workflow for visualization and clustering
      seurat_object.combined <- ScaleData(seurat_object.combined, verbose = FALSE)
      seurat_object.combined <- RunPCA(seurat_object.combined, npcs = dim, verbose = FALSE)
      seurat_object.combined <- RunUMAP(seurat_object.combined, reduction = "pca", dims = 1:dim)
      print(paste(use, dim))
      seurat_object.combined <- FindNeighbors(seurat_object.combined, reduction = "pca", dims = 1:dim)
      seurat_object.combined <- FindClusters(seurat_object.combined, resolution = res)
      print(paste(use, res))
      
      # Save clustering results to metadata matrix
      clustering_results <- seurat_object.combined$seurat_clusters
      
      clustering_results_list[[column_name]] <- clustering_results
    }
  }
}

cluster_results_df <- as.data.frame(clustering_results_list)
print('converted list to dataframe')

#save the data as a dataframe and seurat object
write.csv(cluster_results_df, file = "/path/to/different/clustering/parameter/cluster/assignment.csv")
seurat_object.combined <- AddMetaData(object = seurat_object.combined, metadata = cluster_results_df)
saveRDS(seurat_object.combined, file = "/path/to/seurat_object/cluster_parameters.rds")

#---------compute the ARI-----------#

#use Dune library
library(Dune)

df <- read.csv("/path/to/different/clustering/parameter/cluster/assignment.csv")

#remove "X"
df <- df[, -1]

#get column names 
columns <- colnames(df)

merger <- Dune(clusMat = df %>% select(columns), parallel = TRUE, metric = "ARI")

#save the output
df_merge <- merger$currentMat

#save the ARI results
write.csv(df_merge, file = "/path/to/ari_results.csv", row.names = T)

#compute mean and median ARI scores

#ARI values 
plot <- plotARIs(clusMat = df_merge)

ari_values <- data.frame(plot[["plot_env"]][["Mat"]])

#calculate mean ARI
ari_values$mean_ari <- rowMeans(ari_values)

#top 10 mean ari 
mean_ari <- as.data.frame(ari_values$mean_ari)
rownames(mean_ari) <- rownames(ari_values)
top_n(mean_ari, 10)

#overall mean ari
summarise(ari_values, overall_mean = mean(mean_ari)) #0.8658195

#remove mean_ari from ari_values - remove last column - changes based on number of parameter combinations tested
ari_values <- ari_values[, -c(37)] 

#median ARI

ari_values$row_median = apply(ari_values, 1, median, na.rm=TRUE)

overall_median <- median(ari_values$row_median)  #0.8612033

#top 10 median ARI scores 
median_ari <- as.data.frame(ari_values$row_median)
rownames(median_ari) <- rownames(ari_values)
top_n(median_ari, 10)

#assess top 3 ranked ARI results (based on mean or median)

###########Step 4: ensure clustering repeatability ############

#use the selected clustering results 

#--------Jaccard Index----------#

#perform clustering stability assessment with different clustering resolutions to compare results 
#annotated cell types 

so_low <- readRDS(file = '/path/to/clustering/stability/assessment_low_resolution_seurat_object.rds')

so_high <- readRDS(file = '/path/to/clustering/stability/assessment_high_resolution_seurat_object.rds')


Idents(so_low) <- so_low$cell_type
Idents(so_high) <- ileum_high$cell_type

#get cell types 
clusters_low <- unique(Idents(so_low))
clusters_high <- unique(Idents(so_high))


#calculate jaccard index 
jaccard_fun <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

# Initialize a matrix to store the results

jaccard_matrix_barcode_high_low <- matrix(NA, nrow = length(clusters_low), ncol = length(clusters_high),
                                          dimnames = list(paste0("Low ", clusters_low), paste0("High ", clusters_high)))

# Loop over each cluster in the reference and query sets
for (low_cluster in clusters_low) {
  for (high_cluster in clusters_high) {
    # Get cells in each cluster
    low_cells <- WhichCells(so_low, idents = as.character(low_cluster))
    high_cells <- WhichCells(so_high, idents = as.character(high_cluster))
    
    # Calculate Jaccard index for the cluster pair
    jaccard_matrix_barcode_high_low[paste0("Low ", low_cluster), paste0("High ", high_cluster)] <- jaccard_fun(low_cells, high_cells)
  }
}

jaccard_plot <- pheatmap(jaccard_matrix_barcode_high_low,
                         color = colorRampPalette(c("white", "red"))(50), # Color gradient from white to blue
                         cluster_rows = TRUE,  # Disable clustering if you don't want it
                         cluster_cols = TRUE,
                         display_numbers = TRUE,
                         fontsize = 16,
                         fontsize_number = 12)  # Optionally display the actual Jaccard index values

#--------Spearman Correlation------------#

#use FindAllMarkers function in seurat to compute DEG 

#load in low resolution DEG results for each annotated cell type  

ambig_t_low <- read.csv(file = "/path/to/low_res/DEG/ambiguous_t_cell.csv")
cd4_t_low <- read.csv(file = "/path/to/low_res/DEG/cd4_t_cell.csv")
cd8_t_low <- read.csv(file = "/path/to/low_res/DEG/cd8_t_cell.csv")
doublet_low <- read.csv(file = "/path/to/low_res/DEG/doublet_cluster.csv")
endo_low <- read.csv(file = "/path/to/low_res/DEG/endothelial_cell.csv")
ent_low <- read.csv(file = "/path/to/low_res/DEG/enterocyte.csv")
gob_low <- read.csv(file = "/path/to/low_res/DEG/goblet_cell.csv")
ilc3_low <- read.csv(file = "/path/to/low_res/DEG/ILC3.csv")
inflam_mono_low <- read.csv(file = "/path/to/low_res/DEG/inflammatory_monocyte.csv")
macro_low <- read.csv(file = "/path/to/low_res/DEG/macrophage.csv")
mast_low <- read.csv(file = "/path/to/low_res/DEG/mast_cell.csv")
mem_b_low <- read.csv(file = "/path/to/low_res/DEG/memory_b_cell.csv")
mes_low <- read.csv(file = "/path/to/low_res/DEG/mesenchymal_cell.csv")
mono_low <- read.csv(file = "/path/to/low_res/DEG/monocyte.csv")
naive_b_low <- read.csv(file = "/path/to/low_res/DEG/naive_b_cell.csv")
nk_t_low <- read.csv(file = "/path/to/low_res/DEG/nk_t_cell.csv")
plasma_low <- read.csv(file = "/path/to/low_res/DEG/plasma_cell.csv")
plasmablast_low <- read.csv(file = "/path/to/low_res/DEG/plasmablast.csv")
stem_pan_low <- read.csv(file = "/path/to/low_res/DEG/stem_paneth_cell.csv")
ta_low <- read.csv(file = "/path/to/low_res/DEG/ta_cell.csv")
tuft_low <- read.csv(file = "/path/to/low_res/DEG/tuft_cell.csv")


#low res DEG list
deg_list_celltype_low <- list(ambig_t_low = ambig_t_low, cd4_t_low = cd4_t_low, cd8_t_low = cd8_t_low, doublet_low = doublet_low, endo_low = endo_low, ent_low = ent_low, gob_low = gob_low, ilc3_low = ilc3_low, inflam_mono_low = inflam_mono_low, macro_low = macro_low, mast_low = mast_low, mem_b_low = mem_b_low, mes_low = mes_low, mono_low = mono_low, naive_b_low = naive_b_low, nk_t_low = nk_t_low, plasma_low = plasma_low, plasmablast_low = plasmablast_low, stem_pan_low = stem_pan_low, ta_low = ta_low, tuft_low = tuft_low)


#load in high resolution DEG results for each annotated cell type

ambig_t_high <- read.csv(file = "/path/to/high_res/DEG/ambiguous_t_cell.csv")
cd4_t_high <- read.csv(file = "/path/to/high_res/DEG/cd4_t_cell.csv")
cd8_t_high <- read.csv(file = "/path/to/high_res/DEG/cd8_t_cell.csv")
cyc_nk_t_high <- read.csv(file = "/path/to/high_res/DEG/cycling_nk_t_cell.csv")
doublet_high <- read.csv(file = "/path/to/high_res/DEG/doublet_cluster.csv")
endo_high <- read.csv(file = "/path/to/high_res/DEG/endothelial_cell.csv")
ent_high <- read.csv(file = "/path/to/high_res/DEG/enterocyte.csv")
entero_high <- read.csv(file = "/path/to/high_res/DEG/enteroendocrine.csv")
gc_b_high <- read.csv(file = "/path/to/high_res/DEG/gc_b_cell.csv")
gob_high <- read.csv(file = "/path/to/high_res/DEG/goblet_cell.csv")
ilc3_high <- read.csv(file = "/path/to/high_res/DEG/ilc3.csv")
inflam_macro_high <- read.csv(file = "/path/to/high_res/DEG/inflammatory_macrophage.csv")
inflam_mono_high <- read.csv(file = "/path/to/high_res/DEG/inflammatory_monocyte.csv")
macro_high <- read.csv(file = "/path/to/high_res/DEG/macrophage.csv")
mast_high <- read.csv(file = "/path/to/high_res/DEG/mast_cell.csv")
mem_b_high <- read.csv(file = "/path/to/high_res/DEG/memory_b_cell.csv")
mes_high <- read.csv(file = "/path/to/high_res/DEG/mesenchymal_cell.csv")
mono_high <- read.csv(file = "/path/to/high_res/DEG/monocyte.csv")
naive_b_high <- read.csv(file = "/path/to/high_res/DEG/naive_b_cell.csv")
nk_t_high <- read.csv(file = "/path/to/high_res/DEG/nk_t_cell.csv")
plasma_high <- read.csv(file = "/path/to/high_res/DEG/plasma_cell.csv")
plasmablast_high <- read.csv(file = "/path/to/high_res/DEG/plasmablast.csv")
treg_high <- read.csv(file = "/path/to/high_res/DEG/regulatory_t_cell.csv")
ribo_high <- read.csv(file = "/path/to/high_res/DEG/ribosomal_cluster.csv")
stem_high <- read.csv(file = "/path/to/high_res/DEG/stem_cell.csv")
stem_pan_high <- read.csv(file = "/path/to/high_res/DEG/stem_paneth_cell.csv")
ta_high <- read.csv(file = "/path/to/high_res/DEG/ta_cell.csv")
tuft_high <- read.csv(file = "/path/to/high_res/DEG/tuft_cell.csv")


#high res 
deg_list_celltype_high <- list(ambig_t_high = ambig_t_high, cd4_t_high = cd4_t_high, cd8_t_high = cd8_t_high, cyc_nk_t_high = cyc_nk_t_high, doublet_high = doublet_high, endo_high = endo_high, ent_high = ent_high, entero_high = entero_high, gc_b_high = gc_b_high, gob_high = gob_high, ilc3_high = ilc3_high, inflam_macro_high = inflam_macro_high, inflam_mono_high = inflam_mono_high, macro_high = macro_high, mast_high = mast_high, mem_b_high = mem_b_high, mes_high = mes_high, mono_high = mono_high, naive_b_high = naive_b_high, nk_t_high = nk_t_high, plasma_high = plasma_high, plasmablast_high = plasmablast_high, treg_high = treg_high, ribo_high = ribo_high, stem_high = stem_high, stem_pan_high = stem_pan_high, tuft_high = tuft_high, ta_high = ta_high)

clusters_low <- c('ambig_t_low', 'cd4_t_low', 'cd8_t_low', 'doublet_low', 'endo_low', 'ent_low', 'gob_low', 'ilc3_low', 'inflam_mono_low', 'macro_low', 'mast_low', 'mem_b_low', 'mes_low', 'mono_low', 'naive_b_low', 'nk_t_low', 'plasma_low', 'plasmablast_low', 'stem_pan_low', 'ta_low', 'tuft_low')
clusters_high <- c('ambig_t_high', 'cd4_t_high', 'cd8_t_high', 'cyc_nk_t_high', 'doublet_high', 'endo_high', 'ent_high', 'entero_high', 'gc_b_high', 'gob_high', 'ilc3_high', 'inflam_macro_high', 'inflam_mono_high', 'macro_high', 'mast_high', 'mem_b_high', 'mes_high', 'mono_high', 'naive_b_high', 'nk_t_high', 'plasma_high', 'plasmablast_high', 'treg_high', 'ribo_high', 'stem_high', 'stem_pan_high', 'tuft_high', 'ta_high')


# Initialize a matrix to store the results - high res vs. low res 
spearman_matrix_high_low_log2fc <- matrix(NA, nrow = length(clusters_low), ncol = length(clusters_high),
                                          dimnames = list(paste0("Low ", clusters_low), paste0("High ", clusters_high)))


# Loop through each pair of clusters
for (low_cluster in clusters_low) {
  for (high_cluster in clusters_high) {
    
    # Extract DEG results for each cluster
    low_data <- deg_list_celltype_low[[low_cluster]]
    high_data <- deg_list_celltype_high[[high_cluster]]
    
    
    # Ensure that the data frames contain only the genes present in both clusters
    common_genes <- intersect(low_data$gene, high_data$gene)
    
    # Subset the dataframes to only include the common genes
    low_subset <- low_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    high_subset <- high_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    
    #check if genes are in the same order 
    print(all(low_subset$gene == high_subset$gene))
    # Check if there are any common genes
    if (length(common_genes) > 0) {
      # Calculate Spearman correlation for log2FC values of the common genes
      spearman_corr <- cor(low_subset$avg_log2FC, high_subset$avg_log2FC, method = "spearman")
      
      # Store the correlation in the matrix
      spearman_matrix_high_low_log2fc[paste0("Low ", low_cluster), paste0("High ", high_cluster)] <- spearman_corr
    } else {
      # If no common genes, set correlation to NA
      spearman_matrix_high_low_log2fc[paste0("Low ", low_cluster), paste0("High ", high_cluster)] <- NA
    }
  }
}

spearman_cor_plot <- pheatmap(spearman_matrix_high_low_log2fc,
                              color = colorRampPalette(c("purple", "white", "red"))(50),  
                              cluster_rows = TRUE,  
                              cluster_cols = TRUE,
                              display_numbers = TRUE,
                              number_color = "black", 
                              fontsize = 16,
                              fontsize_number = 12)  # Optionally display the actual Jaccard index values

