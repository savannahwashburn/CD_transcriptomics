#Helmsley Batch 1-13 ILEUM Data Analysis 

#10/17/2024

#load libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(SeuratDisk)
library(SeuratData)
library(RColorBrewer)
library(ArchR)
library(pheatmap)
library(Dune)


#--------Clustering results - no batch effect correction------------#

#load in clustered SO
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_clustered_10_17_24.rds")

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#QC metrics of ileum samples only 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident", cols = c("helm_sam1" = "#BBCCEE", "helm_sam4" = "#BBCCEE", "helm_sam7" = "#BBCCEE", "helm_sam24" = "#99DDFF", "helm_sam29" = "#99DDFF", "helm_sam32" = "#99DDFF", "helm_sam35" = "#99DDFF", "helm_sam44" = "#44BB99", "helm_sam47" = "#44BB99", "helm_sam50" = "#EE8866", "helm_sam53" = "#EE8866", "helm_sam56" = "#EE8866", "helm_sam59" = "#EE8866", "helm_sam62" = "#BBCC33", "helm_sam77" = "#882255", "helm_sam83" = "#882255", "helm_sam88" = "#EECC66", "helm_sam91" = "#EECC66", "helm_sam98" = "#EECC66", "helm_sam101" ="#EECC66", "helm_sam112" = "#FE036A", 'helm_sam119' = "#A8CED2", 'helm_sam134' = '#B2DE7C', 'helm_sam147' = '#B2DE7C', 'helm_sam153' = '#E19390', 'helm_sam162' = '#E19390', 'helm_sam171' = '#BA6BD4')) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident", cols = c("helm_sam1" = "#BBCCEE", "helm_sam4" = "#BBCCEE", "helm_sam7" = "#BBCCEE", "helm_sam24" = "#99DDFF", "helm_sam29" = "#99DDFF", "helm_sam32" = "#99DDFF", "helm_sam35" = "#99DDFF", "helm_sam44" = "#44BB99", "helm_sam47" = "#44BB99", "helm_sam50" = "#EE8866", "helm_sam53" = "#EE8866", "helm_sam56" = "#EE8866", "helm_sam59" = "#EE8866", "helm_sam62" = "#BBCC33", "helm_sam77" = "#882255", "helm_sam83" = "#882255", "helm_sam88" = "#EECC66", "helm_sam91" = "#EECC66", "helm_sam98" = "#EECC66", "helm_sam101" ="#EECC66", "helm_sam112" = "#FE036A", 'helm_sam119' = "#A8CED2", 'helm_sam134' = '#B2DE7C', 'helm_sam147' = '#B2DE7C', 'helm_sam153' = '#E19390', 'helm_sam162' = '#E19390', 'helm_sam171' = '#BA6BD4')) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
VlnPlot(ileum, features = "percent.mt", pt.size = 0, group.by = "orig.ident", cols = c("helm_sam1" = "#BBCCEE", "helm_sam4" = "#BBCCEE", "helm_sam7" = "#BBCCEE", "helm_sam24" = "#99DDFF", "helm_sam29" = "#99DDFF", "helm_sam32" = "#99DDFF", "helm_sam35" = "#99DDFF", "helm_sam44" = "#44BB99", "helm_sam47" = "#44BB99", "helm_sam50" = "#EE8866", "helm_sam53" = "#EE8866", "helm_sam56" = "#EE8866", "helm_sam59" = "#EE8866", "helm_sam62" = "#BBCC33", "helm_sam77" = "#882255", "helm_sam83" = "#882255", "helm_sam88" = "#EECC66", "helm_sam91" = "#EECC66", "helm_sam98" = "#EECC66", "helm_sam101" ="#EECC66", "helm_sam112" = "#FE036A", 'helm_sam119' = "#A8CED2", 'helm_sam134' = '#B2DE7C', 'helm_sam147' = '#B2DE7C', 'helm_sam153' = '#E19390', 'helm_sam162' = '#E19390', 'helm_sam171' = '#BA6BD4')) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"))

library(randomcoloR)
no_of_colors <- 15

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#-----------rPCA preliminary cluster results-------------#
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_rPCA_10_17_24.rds")


ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))


#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

library(randomcoloR)
no_of_colors <- 16

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#rPCA high resolution clustering results 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_rpca_res1_10_21_24.rds")


#--------------Harmony preliminary clustering results----------#

ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_harmony_10_17_24.rds")


ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))


#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

library(randomcoloR)
no_of_colors <- 14

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#looked at marker genes - seems like enterocytes and paneth cells are clustering together. 

DefaultAssay(ileum) <- "RNA"
Idents(ileum) <- ileum$seurat_clusters
#perform find markers on cluster 6 
cluster6 <- FindMarkers(ileum, ident.1 = '6')

#----------CCA preliminary clustering results------------#

ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cca_10_17_24.rds")


ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))


#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

library(randomcoloR)
no_of_colors <- 16

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))


#------------reference ARI results----------#


#merged ARI scores
df <- read.csv('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/df_ri_helm_batch1_13_ileum_ref_dune_merge_10_25_24.csv', row.names = NULL, sep = ",")

df <- df[,-1]
names(df) <- gsub("^paste0\\.", "", names(df))

columns <- colnames(df)

#plot the ARI - without paste0
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/RI_helm_batch1_13_ileum_ref_all_param_merge_10_26_24.pdf", width = 15, height = 15)
plot(plotARIs(df %>% select(columns)) + theme(text=element_text(size=8)) + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) + theme(axis.text.y = element_text(size = 10)))
dev.off()


#----mean and median----#
#ARI values 
plot <- plotARIs(clusMat = df)

ari_values <- data.frame(plot[["plot_env"]][["Mat"]])

#calculate mean ARI
ari_values$mean_ari <- rowMeans(ari_values)

#top 10 mean ari 
mean_ari <- as.data.frame(ari_values$mean_ari)
rownames(mean_ari) <- rownames(ari_values)
top_n(mean_ari, 10)

#overall mean ari
summarise(ari_values, overall_mean = mean(mean_ari)) #0.8213568

#remove mean_ari from ari_values 
ari_values <- ari_values[, -c(37)]

#median ARI

ari_values$row_median = apply(ari_values, 1, median, na.rm=TRUE)

overall_median <- median(ari_values$row_median)  #0.8355671


#take average of median 
summarise(ari_values, overall_mean = mean(row_median)) #0.8235083

#top 10 median ARI scores 
median_ari <- as.data.frame(ari_values$row_median)
rownames(median_ari) <- rownames(ari_values)
top_n(median_ari, 10) 


#----------Low resolution stability assessment----------#

#-----Query 1------#

cluster1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub_res20_1_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #9813

#cluster 1
Idents(cluster1) <- cluster1$seurat_clusters
sub1 <- WhichCells(cluster1, idents = c('1'))

Idents(cluster1) <- cluster1$predicted.id
pred1 <- WhichCells(cluster1, idents = '1')

int1 <- intersect(sub1, pred1) #7548

#cluster 2
Idents(cluster1) <- cluster1$seurat_clusters
sub2 <- WhichCells(cluster1, idents = c('2', '9'))

Idents(cluster1) <- cluster1$predicted.id
pred2 <- WhichCells(cluster1, idents = '2')

int2 <- intersect(sub2, pred2) #6162

#cluster 3
Idents(cluster1) <- cluster1$seurat_clusters
sub3 <- WhichCells(cluster1, idents = c('3'))

Idents(cluster1) <- cluster1$predicted.id
pred3 <- WhichCells(cluster1, idents = '3')

int3 <- intersect(sub3, pred3) #4348

#cluster 4
Idents(cluster1) <- cluster1$seurat_clusters
sub4 <- WhichCells(cluster1, idents = c('4'))

Idents(cluster1) <- cluster1$predicted.id
pred4 <- WhichCells(cluster1, idents = '4')

int4 <- intersect(sub4, pred4) #3203

#cluster 5
Idents(cluster1) <- cluster1$seurat_clusters
sub5 <- WhichCells(cluster1, idents = c('5'))

Idents(cluster1) <- cluster1$predicted.id
pred5 <- WhichCells(cluster1, idents = '5')

int5 <- intersect(sub5, pred5) #2583

#cluster 6
Idents(cluster1) <- cluster1$seurat_clusters
sub6 <- WhichCells(cluster1, idents = c('6'))

Idents(cluster1) <- cluster1$predicted.id
pred6 <- WhichCells(cluster1, idents = '6')

int6 <- intersect(sub6, pred6) #2427

#cluster 7
Idents(cluster1) <- cluster1$seurat_clusters
sub7 <- WhichCells(cluster1, idents = c('1', '7'))

Idents(cluster1) <- cluster1$predicted.id
pred7<- WhichCells(cluster1, idents = '7')

int7 <- intersect(sub7, pred7) #2215

#cluster 8 
Idents(cluster1) <- cluster1$seurat_clusters
sub8 <- WhichCells(cluster1, idents = c('8'))

Idents(cluster1) <- cluster1$predicted.id
pred8 <- WhichCells(cluster1, idents = '8')

int8 <- intersect(sub8, pred8) #1411

#cluster 9
Idents(cluster1) <- cluster1$seurat_clusters
sub9 <- WhichCells(cluster1, idents = c('10'))

Idents(cluster1) <- cluster1$predicted.id
pred9 <- WhichCells(cluster1, idents = '9')

int9 <- intersect(sub9, pred9) #992

#cluster 10
Idents(cluster1) <- cluster1$seurat_clusters
sub10 <- WhichCells(cluster1, idents = c('11'))

Idents(cluster1) <- cluster1$predicted.id
pred10 <- WhichCells(cluster1, idents = '10')

int10 <- intersect(sub10, pred10) #829

#cluster 11
Idents(cluster1) <- cluster1$seurat_clusters
sub11 <- WhichCells(cluster1, idents = c('12'))

Idents(cluster1) <- cluster1$predicted.id
pred11 <- WhichCells(cluster1, idents = '11')

int11 <- intersect(sub11, pred11) #294

#cluster 12
Idents(cluster1) <- cluster1$seurat_clusters
sub12 <- WhichCells(cluster1, idents = c('14'))

Idents(cluster1) <- cluster1$predicted.id
pred12 <- WhichCells(cluster1, idents = '12')

int12 <- intersect(sub12, pred12) #180

#cluster 13
Idents(cluster1) <- cluster1$seurat_clusters
sub13 <- WhichCells(cluster1, idents = c('13'))

Idents(cluster1) <- cluster1$predicted.id
pred13 <- WhichCells(cluster1, idents = '13')

int13 <- intersect(sub13, pred13) #206

#cluster 14
Idents(cluster1) <- cluster1$seurat_clusters
sub14 <- WhichCells(cluster1, idents = c('15'))

Idents(cluster1) <- cluster1$predicted.id
pred14 <- WhichCells(cluster1, idents = '14')

int14 <- intersect(sub14, pred14) #104

#cluster 15
Idents(cluster1) <- cluster1$seurat_clusters
sub15 <- WhichCells(cluster1, idents = c('10'))

Idents(cluster1) <- cluster1$predicted.id
pred15 <- WhichCells(cluster1, idents = '15')

int15 <- intersect(sub15, pred15) #33

#total 
#pre:43892
#42348
total1 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total1,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query1_low_res_barcodes_10_23_24.csv', row.names=F)

#-----Query 2------#

cluster2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub2_res20_1_low_res_cluster_pred.rds")

cluster2$predicted.id <- factor(cluster2$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster2, reduction = 'umap', label = TRUE) + DimPlot(cluster2, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table2 <- table(cluster2$predicted.id, cluster2$seurat_clusters)
table2 <- as.data.frame.matrix(table2)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table2_prop <- as.data.frame.matrix(apply(table2, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster2) <- cluster2$seurat_clusters
sub0 <- WhichCells(cluster2, idents = c('0'))

Idents(cluster2) <- cluster2$predicted.id
pred0 <- WhichCells(cluster2, idents = '0')

int0 <- intersect(sub0, pred0) #9396

#cluster 1
Idents(cluster2) <- cluster2$seurat_clusters
sub1 <- WhichCells(cluster2, idents = c('1', '2'))

Idents(cluster2) <- cluster2$predicted.id
pred1 <- WhichCells(cluster2, idents = '1')

int1 <- intersect(sub1, pred1) #7992

#cluster 2
Idents(cluster2) <- cluster2$seurat_clusters
sub2 <- WhichCells(cluster2, idents = c('2', '9'))

Idents(cluster2) <- cluster2$predicted.id
pred2 <- WhichCells(cluster2, idents = '2')

int2 <- intersect(sub2, pred2) #5927

#cluster 3
Idents(cluster2) <- cluster2$seurat_clusters
sub3 <- WhichCells(cluster2, idents = c('3'))

Idents(cluster2) <- cluster2$predicted.id
pred3 <- WhichCells(cluster2, idents = '3')

int3 <- intersect(sub3, pred3) #4335

#cluster 4
Idents(cluster2) <- cluster2$seurat_clusters
sub4 <- WhichCells(cluster2, idents = c('4'))

Idents(cluster2) <- cluster2$predicted.id
pred4 <- WhichCells(cluster2, idents = '4')

int4 <- intersect(sub4, pred4) #3166

#cluster 5
Idents(cluster2) <- cluster2$seurat_clusters
sub5 <- WhichCells(cluster2, idents = c('6'))

Idents(cluster2) <- cluster2$predicted.id
pred5 <- WhichCells(cluster2, idents = '5')

int5 <- intersect(sub5, pred5) #2590

#cluster 6
Idents(cluster2) <- cluster2$seurat_clusters
sub6 <- WhichCells(cluster2, idents = c('5'))

Idents(cluster2) <- cluster2$predicted.id
pred6 <- WhichCells(cluster2, idents = '6')

int6 <- intersect(sub6, pred6) #2531

#cluster 7
Idents(cluster2) <- cluster2$seurat_clusters
sub7 <- WhichCells(cluster2, idents = c('7'))

Idents(cluster2) <- cluster2$predicted.id
pred7 <- WhichCells(cluster2, idents = '7')

int7 <- intersect(sub7, pred7) #2211

#cluster 8
Idents(cluster2) <- cluster2$seurat_clusters
sub8 <- WhichCells(cluster2, idents = c('8'))

Idents(cluster2) <- cluster2$predicted.id
pred8 <- WhichCells(cluster2, idents = '8')

int8 <- intersect(sub8, pred8) #1432

#cluster 9 
Idents(cluster2) <- cluster2$seurat_clusters
sub9 <- WhichCells(cluster2, idents = c('10'))

Idents(cluster2) <- cluster2$predicted.id
pred9 <- WhichCells(cluster2, idents = '9')

int9 <- intersect(sub9, pred9) #1017

#cluster 10 
Idents(cluster2) <- cluster2$seurat_clusters
sub10 <- WhichCells(cluster2, idents = c('11'))

Idents(cluster2) <- cluster2$predicted.id
pred10 <- WhichCells(cluster2, idents = '10')

int10 <- intersect(sub10, pred10) #770

#cluster 11
Idents(cluster2) <- cluster2$seurat_clusters
sub11 <- WhichCells(cluster2, idents = c('12'))

Idents(cluster2) <- cluster2$predicted.id
pred11 <- WhichCells(cluster2, idents = '11')

int11 <- intersect(sub11, pred11) #352

#cluster 12
Idents(cluster2) <- cluster2$seurat_clusters
sub12 <- WhichCells(cluster2, idents = c('13'))

Idents(cluster2) <- cluster2$predicted.id
pred12 <- WhichCells(cluster2, idents = '12')

int12 <- intersect(sub12, pred12) #214

#cluster 13
Idents(cluster2) <- cluster2$seurat_clusters
sub13 <- WhichCells(cluster2, idents = c('14'))

Idents(cluster2) <- cluster2$predicted.id
pred13 <- WhichCells(cluster2, idents = '13')

int13 <- intersect(sub13, pred13) #158

#cluster 14
Idents(cluster2) <- cluster2$seurat_clusters
sub14 <- WhichCells(cluster2, idents = c('15'))

Idents(cluster2) <- cluster2$predicted.id
pred14 <- WhichCells(cluster2, idents = '14')

int14 <- intersect(sub14, pred14) #68

#cluster 15
Idents(cluster2) <- cluster2$seurat_clusters
sub15 <- WhichCells(cluster2, idents = c('4'))

Idents(cluster2) <- cluster2$predicted.id
pred15 <- WhichCells(cluster2, idents = '15')

int15 <- intersect(sub15, pred15) #41

#total 
#pre:43893
#42200
total2 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total2,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query2_low_res_barcodes_10_23_24.csv', row.names=F)

#-----Query 3------#

cluster3 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub_res20_2_low_res_cluster_pred.rds")

cluster3$predicted.id <- factor(cluster3$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster3, reduction = 'umap', label = TRUE) + DimPlot(cluster3, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table3 <- table(cluster3$predicted.id, cluster3$seurat_clusters)
table3 <- as.data.frame.matrix(table3)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table3_prop <- as.data.frame.matrix(apply(table3, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster3) <- cluster3$seurat_clusters
sub0 <- WhichCells(cluster3, idents = c('0'))

Idents(cluster3) <- cluster3$predicted.id
pred0 <- WhichCells(cluster3, idents = '0')

int0 <- intersect(sub0, pred0) #8802

#cluster 1 
Idents(cluster3) <- cluster3$seurat_clusters
sub1 <- WhichCells(cluster3, idents = c('1', '2'))

Idents(cluster3) <- cluster3$predicted.id
pred1 <- WhichCells(cluster3, idents = '1')

int1 <- intersect(sub1, pred1) #8129

#cluster 2
Idents(cluster3) <- cluster3$seurat_clusters
sub2 <- WhichCells(cluster3, idents = c('2', '9'))

Idents(cluster3) <- cluster3$predicted.id
pred2 <- WhichCells(cluster3, idents = '2')

int2 <- intersect(sub2, pred2) #5858

#cluster 3
Idents(cluster3) <- cluster3$seurat_clusters
sub3 <- WhichCells(cluster3, idents = c('3'))

Idents(cluster3) <- cluster3$predicted.id
pred3 <- WhichCells(cluster3, idents = '3')

int3 <- intersect(sub3, pred3) #4234

#cluster 4
Idents(cluster3) <- cluster3$seurat_clusters
sub4 <- WhichCells(cluster3, idents = c('4'))

Idents(cluster3) <- cluster3$predicted.id
pred4 <- WhichCells(cluster3, idents = '4')

int4 <- intersect(sub4, pred4) #3133

#cluster 5
Idents(cluster3) <- cluster3$seurat_clusters
sub5 <- WhichCells(cluster3, idents = c('5'))

Idents(cluster3) <- cluster3$predicted.id
pred5 <- WhichCells(cluster3, idents = '5')

int5 <- intersect(sub5, pred5) #2607

#cluster 6
Idents(cluster3) <- cluster3$seurat_clusters
sub6 <- WhichCells(cluster3, idents = c('6'))

Idents(cluster3) <- cluster3$predicted.id
pred6 <- WhichCells(cluster3, idents = '6')

int6 <- intersect(sub6, pred6) #2529

#cluster 7
Idents(cluster3) <- cluster3$seurat_clusters
sub7 <- WhichCells(cluster3, idents = c('7'))

Idents(cluster3) <- cluster3$predicted.id
pred7 <- WhichCells(cluster3, idents = '7')

int7 <- intersect(sub7, pred7) #2256

#cluster 8
Idents(cluster3) <- cluster3$seurat_clusters
sub8 <- WhichCells(cluster3, idents = c('8'))

Idents(cluster3) <- cluster3$predicted.id
pred8 <- WhichCells(cluster3, idents = '8')

int8 <- intersect(sub8, pred8) #1362

#cluster 9 
Idents(cluster3) <- cluster3$seurat_clusters
sub9 <- WhichCells(cluster3, idents = c('10'))

Idents(cluster3) <- cluster3$predicted.id
pred9 <- WhichCells(cluster3, idents = '9')

int9 <- intersect(sub9, pred9) #1094

#cluster 10 
Idents(cluster3) <- cluster3$seurat_clusters
sub10 <- WhichCells(cluster3, idents = c('12'))

Idents(cluster3) <- cluster3$predicted.id
pred10 <- WhichCells(cluster3, idents = '10')

int10 <- intersect(sub10, pred10) #816

#cluster 11
Idents(cluster3) <- cluster3$seurat_clusters
sub11 <- WhichCells(cluster3, idents = c('12'))

Idents(cluster3) <- cluster3$predicted.id
pred11 <- WhichCells(cluster3, idents = '11')

int11 <- intersect(sub11, pred11) #308

#cluster 12
Idents(cluster3) <- cluster3$seurat_clusters
sub12 <- WhichCells(cluster3, idents = c('13'))

Idents(cluster3) <- cluster3$predicted.id
pred12 <- WhichCells(cluster3, idents = '12')

int12 <- intersect(sub12, pred12) #213

#cluster 13
Idents(cluster3) <- cluster3$seurat_clusters
sub13 <- WhichCells(cluster3, idents = c('14'))

Idents(cluster3) <- cluster3$predicted.id
pred13 <- WhichCells(cluster3, idents = '13')

int13 <- intersect(sub13, pred13) #184

#cluster 14
Idents(cluster3) <- cluster3$seurat_clusters
sub14 <- WhichCells(cluster3, idents = c('3'))

Idents(cluster3) <- cluster3$predicted.id
pred14 <- WhichCells(cluster3, idents = '14')

int14 <- intersect(sub14, pred14) #77

#cluster 15
Idents(cluster3) <- cluster3$seurat_clusters
sub15 <- WhichCells(cluster3, idents = c('4', '10'))

Idents(cluster3) <- cluster3$predicted.id
pred15 <- WhichCells(cluster3, idents = '15')

int15 <- intersect(sub15, pred15) #34

#total 
#pre:43892
#41636
total3 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total3,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query3_low_res_barcodes_10_23_24.csv', row.names=F)


#-----Query 4------#

cluster4 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub2_res20_2_low_res_cluster_pred.rds")

cluster4$predicted.id <- factor(cluster4$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster4, reduction = 'umap', label = TRUE) + DimPlot(cluster4, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table4 <- table(cluster4$predicted.id, cluster4$seurat_clusters)
table4 <- as.data.frame.matrix(table4)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table4_prop <- as.data.frame.matrix(apply(table4, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster4) <- cluster4$seurat_clusters
sub0 <- WhichCells(cluster4, idents = c('0'))

Idents(cluster4) <- cluster4$predicted.id
pred0 <- WhichCells(cluster4, idents = '0')

int0 <- intersect(sub0, pred0) #9148

#cluster 1
Idents(cluster4) <- cluster4$seurat_clusters
sub1 <- WhichCells(cluster4, idents = c('1', '2'))

Idents(cluster4) <- cluster4$predicted.id
pred1 <- WhichCells(cluster4, idents = '1')

int1 <- intersect(sub1, pred1) #8121

#cluster 2
Idents(cluster4) <- cluster4$seurat_clusters
sub2 <- WhichCells(cluster4, idents = c('2', '8'))

Idents(cluster4) <- cluster4$predicted.id
pred2 <- WhichCells(cluster4, idents = '2')

int2 <- intersect(sub2, pred2) #6467

#cluster 3
Idents(cluster4) <- cluster4$seurat_clusters
sub3 <- WhichCells(cluster4, idents = c('3'))

Idents(cluster4) <- cluster4$predicted.id
pred3 <- WhichCells(cluster4, idents = '3')

int3 <- intersect(sub3, pred3) #4392

#cluster 4
Idents(cluster4) <- cluster4$seurat_clusters
sub4 <- WhichCells(cluster4, idents = c('4'))

Idents(cluster4) <- cluster4$predicted.id
pred4 <- WhichCells(cluster4, idents = '4')

int4 <- intersect(sub4, pred4) #3210

#cluster 5
Idents(cluster4) <- cluster4$seurat_clusters
sub5 <- WhichCells(cluster4, idents = c('5'))

Idents(cluster4) <- cluster4$predicted.id
pred5 <- WhichCells(cluster4, idents = '5')

int5 <- intersect(sub5, pred5) #2534

#cluster 6
Idents(cluster4) <- cluster4$seurat_clusters
sub6 <- WhichCells(cluster4, idents = c('6'))

Idents(cluster4) <- cluster4$predicted.id
pred6 <- WhichCells(cluster4, idents = '6')

int6 <- intersect(sub6, pred6) #2465

#cluster 7
Idents(cluster4) <- cluster4$seurat_clusters
sub7 <- WhichCells(cluster4, idents = c('7'))

Idents(cluster4) <- cluster4$predicted.id
pred7 <- WhichCells(cluster4, idents = '7')

int7 <- intersect(sub7, pred7) #2246

#cluster 8
Idents(cluster4) <- cluster4$seurat_clusters
sub8 <- WhichCells(cluster4, idents = c('9'))

Idents(cluster4) <- cluster4$predicted.id
pred8 <- WhichCells(cluster4, idents = '8')

int8 <- intersect(sub8, pred8) #1470

#cluster 9
Idents(cluster4) <- cluster4$seurat_clusters
sub9 <- WhichCells(cluster4, idents = c('1', '11'))

Idents(cluster4) <- cluster4$predicted.id
pred9 <- WhichCells(cluster4, idents = '9')

int9 <- intersect(sub9, pred9) #1083

#cluster 10 
Idents(cluster4) <- cluster4$seurat_clusters
sub10 <- WhichCells(cluster4, idents = c('10'))

Idents(cluster4) <- cluster4$predicted.id
pred10 <- WhichCells(cluster4, idents = '10')

int10 <- intersect(sub10, pred10) #794

#cluster 11
Idents(cluster4) <- cluster4$seurat_clusters
sub11 <- WhichCells(cluster4, idents = c('12'))

Idents(cluster4) <- cluster4$predicted.id
pred11 <- WhichCells(cluster4, idents = '11')

int11 <- intersect(sub11, pred11) #320

#cluster 12
Idents(cluster4) <- cluster4$seurat_clusters
sub12 <- WhichCells(cluster4, idents = c('14'))

Idents(cluster4) <- cluster4$predicted.id
pred12 <- WhichCells(cluster4, idents = '12')

int12 <- intersect(sub12, pred12) #161

#cluster 13
Idents(cluster4) <- cluster4$seurat_clusters
sub13 <- WhichCells(cluster4, idents = c('13'))

Idents(cluster4) <- cluster4$predicted.id
pred13 <- WhichCells(cluster4, idents = '13')

int13 <- intersect(sub13, pred13) #172

#cluster 14 
Idents(cluster4) <- cluster4$seurat_clusters
sub14 <- WhichCells(cluster4, idents = c('15'))

Idents(cluster4) <- cluster4$predicted.id
pred14 <- WhichCells(cluster4, idents = '14')

int14 <- intersect(sub14, pred14) #80

#cluster 15
Idents(cluster4) <- cluster4$seurat_clusters
sub15 <- WhichCells(cluster4, idents = c('4'))

Idents(cluster4) <- cluster4$predicted.id
pred15 <- WhichCells(cluster4, idents = '15')

int15 <- intersect(sub15, pred15) #39

#total 
#pre:43892
#42702
total4 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total4,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query4_low_res_barcodes_10_23_24.csv', row.names=F)

#-------Query 5--------#

cluster5 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub_res20_3_low_res_cluster_pred.rds")

cluster5$predicted.id <- factor(cluster5$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster5, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster5, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table5 <- table(cluster5$predicted.id, cluster5$seurat_clusters)
table5 <- as.data.frame.matrix(table5)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table5_prop <- as.data.frame.matrix(apply(table5, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster5) <- cluster5$seurat_clusters
sub0 <- WhichCells(cluster5, idents = c('0'))

Idents(cluster5) <- cluster5$predicted.id
pred0 <- WhichCells(cluster5, idents = '0')

int0 <- intersect(sub0, pred0) #9118

#cluster 1
Idents(cluster5) <- cluster5$seurat_clusters
sub1 <- WhichCells(cluster5, idents = c('1', '2'))

Idents(cluster5) <- cluster5$predicted.id
pred1 <- WhichCells(cluster5, idents = '1')

int1 <- intersect(sub1, pred1) #7559

#cluster 2 
Idents(cluster5) <- cluster5$seurat_clusters
sub2 <- WhichCells(cluster5, idents = c('1', '8'))

Idents(cluster5) <- cluster5$predicted.id
pred2 <- WhichCells(cluster5, idents = '2')

int2 <- intersect(sub2, pred2) #6497

#cluster 3
Idents(cluster5) <- cluster5$seurat_clusters
sub3 <- WhichCells(cluster5, idents = c('3'))

Idents(cluster5) <- cluster5$predicted.id
pred3 <- WhichCells(cluster5, idents = '3')

int3 <- intersect(sub3, pred3) #4326

#cluster 4
Idents(cluster5) <- cluster5$seurat_clusters
sub4 <- WhichCells(cluster5, idents = c('4'))

Idents(cluster5) <- cluster5$predicted.id
pred4 <- WhichCells(cluster5, idents = '4')

int4 <- intersect(sub4, pred4) #3239

#cluster 5
Idents(cluster5) <- cluster5$seurat_clusters
sub5 <- WhichCells(cluster5, idents = c('5'))

Idents(cluster5) <- cluster5$predicted.id
pred5 <- WhichCells(cluster5, idents = '5')

int5 <- intersect(sub5, pred5) #2587

#cluster 6
Idents(cluster5) <- cluster5$seurat_clusters
sub6 <- WhichCells(cluster5, idents = c('6'))

Idents(cluster5) <- cluster5$predicted.id
pred6 <- WhichCells(cluster5, idents = '6')

int6 <- intersect(sub6, pred6) #2475

#cluster 7
Idents(cluster5) <- cluster5$seurat_clusters
sub7 <- WhichCells(cluster5, idents = c('7'))

Idents(cluster5) <- cluster5$predicted.id
pred7 <- WhichCells(cluster5, idents = '7')

int7 <- intersect(sub7, pred7) #2182

#cluster 8
Idents(cluster5) <- cluster5$seurat_clusters
sub8 <- WhichCells(cluster5, idents = c('9'))

Idents(cluster5) <- cluster5$predicted.id
pred8 <- WhichCells(cluster5, idents = '8')

int8 <- intersect(sub8, pred8) #1434

#cluster 9
Idents(cluster5) <- cluster5$seurat_clusters
sub9 <- WhichCells(cluster5, idents = c('10'))

Idents(cluster5) <- cluster5$predicted.id
pred9 <- WhichCells(cluster5, idents = '9')

int9 <- intersect(sub9, pred9) #1068

#cluster 10
Idents(cluster5) <- cluster5$seurat_clusters
sub10 <- WhichCells(cluster5, idents = c('11'))

Idents(cluster5) <- cluster5$predicted.id
pred10 <- WhichCells(cluster5, idents = '10')

int10 <- intersect(sub10, pred10) #825

#cluster 11
Idents(cluster5) <- cluster5$seurat_clusters
sub11 <- WhichCells(cluster5, idents = c('12'))

Idents(cluster5) <- cluster5$predicted.id
pred11 <- WhichCells(cluster5, idents = '11')

int11 <- intersect(sub11, pred11) #354

#cluster 12
Idents(cluster5) <- cluster5$seurat_clusters
sub12 <- WhichCells(cluster5, idents = c('14'))

Idents(cluster5) <- cluster5$predicted.id
pred12 <- WhichCells(cluster5, idents = '12')

int12 <- intersect(sub12, pred12) #205

#cluster 13
Idents(cluster5) <- cluster5$seurat_clusters
sub13 <- WhichCells(cluster5, idents = c('13'))

Idents(cluster5) <- cluster5$predicted.id
pred13 <- WhichCells(cluster5, idents = '13')

int13 <- intersect(sub13, pred13) #183

#cluster 14
Idents(cluster5) <- cluster5$seurat_clusters
sub14 <- WhichCells(cluster5, idents = c('15'))

Idents(cluster5) <- cluster5$predicted.id
pred14 <- WhichCells(cluster5, idents = '14')

int14 <- intersect(sub14, pred14) #92

#cluster 15
Idents(cluster5) <- cluster5$seurat_clusters
sub15 <- WhichCells(cluster5, idents = c('10'))

Idents(cluster5) <- cluster5$predicted.id
pred15 <- WhichCells(cluster5, idents = '15')

int15 <- intersect(sub15, pred15) #37

#total 
#pre:43892
#42234
total5 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total5,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query5_low_res_barcodes_10_23_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster5$query_label <- paste(cluster5$seurat_clusters, "query", sep = "_")
cluster5$ref_label <- paste(cluster5$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster5$ref_label), paste0(cluster5$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-------Query 6--------#

cluster6 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub2_res20_3_low_res_cluster_pred.rds")

cluster6$predicted.id <- factor(cluster6$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster6, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster6, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table6 <- table(cluster6$predicted.id, cluster6$seurat_clusters)
table6 <- as.data.frame.matrix(table6)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table6_prop <- as.data.frame.matrix(apply(table6, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster6) <- cluster6$seurat_clusters
sub0 <- WhichCells(cluster6, idents = c('0'))

Idents(cluster6) <- cluster6$predicted.id
pred0 <- WhichCells(cluster6, idents = '0')

int0 <- intersect(sub0, pred0) #9823

#cluster 1
Idents(cluster6) <- cluster6$seurat_clusters
sub1 <- WhichCells(cluster6, idents = c('1'))

Idents(cluster6) <- cluster6$predicted.id
pred1 <- WhichCells(cluster6, idents = '1')

int1 <- intersect(sub1, pred1) #8166

#cluster 2
Idents(cluster6) <- cluster6$seurat_clusters
sub2 <- WhichCells(cluster6, idents = c('1', '2'))

Idents(cluster6) <- cluster6$predicted.id
pred2 <- WhichCells(cluster6, idents = '2')

int2 <- intersect(sub2, pred2) #6401

#cluster 3
Idents(cluster6) <- cluster6$seurat_clusters
sub3 <- WhichCells(cluster6, idents = c('3'))

Idents(cluster6) <- cluster6$predicted.id
pred3 <- WhichCells(cluster6, idents = '3')

int3 <- intersect(sub3, pred3) #4407

#cluster 4
Idents(cluster6) <- cluster6$seurat_clusters
sub4 <- WhichCells(cluster6, idents = c('4'))

Idents(cluster6) <- cluster6$predicted.id
pred4 <- WhichCells(cluster6, idents = '4')

int4 <- intersect(sub4, pred4) #3095

#cluster 5
Idents(cluster6) <- cluster6$seurat_clusters
sub5 <- WhichCells(cluster6, idents = c('5'))

Idents(cluster6) <- cluster6$predicted.id
pred5 <- WhichCells(cluster6, idents = '5')

int5 <- intersect(sub5, pred5) #2573

#cluster 6
Idents(cluster6) <- cluster6$seurat_clusters
sub6 <- WhichCells(cluster6, idents = c('6'))

Idents(cluster6) <- cluster6$predicted.id
pred6 <- WhichCells(cluster6, idents = '6')

int6 <- intersect(sub6, pred6) #2452

#cluster 7
Idents(cluster6) <- cluster6$seurat_clusters
sub7 <- WhichCells(cluster6, idents = c('7'))

Idents(cluster6) <- cluster6$predicted.id
pred7 <- WhichCells(cluster6, idents = '7')

int7 <- intersect(sub7, pred7) #2254

#cluster 8
Idents(cluster6) <- cluster6$seurat_clusters
sub8 <- WhichCells(cluster6, idents = c('8'))

Idents(cluster6) <- cluster6$predicted.id
pred8 <- WhichCells(cluster6, idents = '8')

int8 <- intersect(sub8, pred8) #1415

#cluster 9 
Idents(cluster6) <- cluster6$seurat_clusters
sub9 <- WhichCells(cluster6, idents = c('9'))

Idents(cluster6) <- cluster6$predicted.id
pred9 <- WhichCells(cluster6, idents = '9')

int9 <- intersect(sub9, pred9) #1093

#cluster 10
Idents(cluster6) <- cluster6$seurat_clusters
sub10 <- WhichCells(cluster6, idents = c('10'))

Idents(cluster6) <- cluster6$predicted.id
pred10 <- WhichCells(cluster6, idents = '10')

int10 <- intersect(sub10, pred10) #784

#cluster 11
Idents(cluster6) <- cluster6$seurat_clusters
sub11 <- WhichCells(cluster6, idents = c('10'))

Idents(cluster6) <- cluster6$predicted.id
pred11 <- WhichCells(cluster6, idents = '11')

int11 <- intersect(sub11, pred11) #287

#cluster 12
Idents(cluster6) <- cluster6$seurat_clusters
sub12 <- WhichCells(cluster6, idents = c('12'))

Idents(cluster6) <- cluster6$predicted.id
pred12 <- WhichCells(cluster6, idents = '12')

int12 <- intersect(sub12, pred12) #195

#cluster 13
Idents(cluster6) <- cluster6$seurat_clusters
sub13 <- WhichCells(cluster6, idents = c('11'))

Idents(cluster6) <- cluster6$predicted.id
pred13 <- WhichCells(cluster6, idents = '13')

int13 <- intersect(sub13, pred13) #195

#cluster 14
Idents(cluster6) <- cluster6$seurat_clusters
sub14 <- WhichCells(cluster6, idents = c('13'))

Idents(cluster6) <- cluster6$predicted.id
pred14 <- WhichCells(cluster6, idents = '14')

int14 <- intersect(sub14, pred14) #79

#cluster 15 
Idents(cluster6) <- cluster6$seurat_clusters
sub15 <- WhichCells(cluster6, idents = c('4'))

Idents(cluster6) <- cluster6$predicted.id
pred15 <- WhichCells(cluster6, idents = '15')

int15 <- intersect(sub15, pred15) #41

#total 
#pre:43893
#43260
total6 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total6,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query6_low_res_barcodes_10_23_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster6$query_label <- paste(cluster6$seurat_clusters, "query", sep = "_")
cluster6$ref_label <- paste(cluster6$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster6$ref_label), paste0(cluster6$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)


#-------Query 7--------#

cluster7 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub_res20_4_low_res_cluster_pred.rds")

cluster7$predicted.id <- factor(cluster7$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster7, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster7, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table7 <- table(cluster7$predicted.id, cluster7$seurat_clusters)
table7 <- as.data.frame.matrix(table7)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table7_prop <- as.data.frame.matrix(apply(table7, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster7) <- cluster7$seurat_clusters
sub0 <- WhichCells(cluster7, idents = c('0'))

Idents(cluster7) <- cluster7$predicted.id
pred0 <- WhichCells(cluster7, idents = '0')

int0 <- intersect(sub0, pred0) #9063

#cluster 1
Idents(cluster7) <- cluster7$seurat_clusters
sub1 <- WhichCells(cluster7, idents = c('1', '2'))

Idents(cluster7) <- cluster7$predicted.id
pred1 <- WhichCells(cluster7, idents = '1')

int1 <- intersect(sub1, pred1) #7939

#cluster 2
Idents(cluster7) <- cluster7$seurat_clusters
sub2 <- WhichCells(cluster7, idents = c('2', '8'))

Idents(cluster7) <- cluster7$predicted.id
pred2 <- WhichCells(cluster7, idents = '2')

int2 <- intersect(sub2, pred2) #5984

#cluster 3
Idents(cluster7) <- cluster7$seurat_clusters
sub3 <- WhichCells(cluster7, idents = c('3'))

Idents(cluster7) <- cluster7$predicted.id
pred3 <- WhichCells(cluster7, idents = '3')

int3 <- intersect(sub3, pred3) #4406

#cluster 4
Idents(cluster7) <- cluster7$seurat_clusters
sub4 <- WhichCells(cluster7, idents = c('4'))

Idents(cluster7) <- cluster7$predicted.id
pred4 <- WhichCells(cluster7, idents = '4')

int4 <- intersect(sub4, pred4) #3146

#cluster 5
Idents(cluster7) <- cluster7$seurat_clusters
sub5 <- WhichCells(cluster7, idents = c('6'))

Idents(cluster7) <- cluster7$predicted.id
pred5 <- WhichCells(cluster7, idents = '5')

int5 <- intersect(sub5, pred5) #2575

#cluster 6
Idents(cluster7) <- cluster7$seurat_clusters
sub6 <- WhichCells(cluster7, idents = c('5'))

Idents(cluster7) <- cluster7$predicted.id
pred6 <- WhichCells(cluster7, idents = '6')

int6 <- intersect(sub6, pred6) #2463

#cluster 7
Idents(cluster7) <- cluster7$seurat_clusters
sub7 <- WhichCells(cluster7, idents = c('7'))

Idents(cluster7) <- cluster7$predicted.id
pred7 <- WhichCells(cluster7, idents = '7')

int7 <- intersect(sub7, pred7) #2216

#cluster 8
Idents(cluster7) <- cluster7$seurat_clusters
sub8 <- WhichCells(cluster7, idents = c('9'))

Idents(cluster7) <- cluster7$predicted.id
pred8 <- WhichCells(cluster7, idents = '8')

int8 <- intersect(sub8, pred8) #1387

#cluster 9 
Idents(cluster7) <- cluster7$seurat_clusters
sub9 <- WhichCells(cluster7, idents = c('11', '12'))

Idents(cluster7) <- cluster7$predicted.id
pred9 <- WhichCells(cluster7, idents = '9')

int9 <- intersect(sub9, pred9) #1092

#cluster 10
Idents(cluster7) <- cluster7$seurat_clusters
sub10 <- WhichCells(cluster7, idents = c('10'))

Idents(cluster7) <- cluster7$predicted.id
pred10 <- WhichCells(cluster7, idents = '10')

int10 <- intersect(sub10, pred10) #808

#cluster 11
Idents(cluster7) <- cluster7$seurat_clusters
sub11 <- WhichCells(cluster7, idents = c('13'))

Idents(cluster7) <- cluster7$predicted.id
pred11 <- WhichCells(cluster7, idents = '11')

int11 <- intersect(sub11, pred11) #324

#cluster 12
Idents(cluster7) <- cluster7$seurat_clusters
sub12 <- WhichCells(cluster7, idents = c('14'))

Idents(cluster7) <- cluster7$predicted.id
pred12 <- WhichCells(cluster7, idents = '12')

int12 <- intersect(sub12, pred12) #211

#cluster 13
Idents(cluster7) <- cluster7$seurat_clusters
sub13 <- WhichCells(cluster7, idents = c('15'))

Idents(cluster7) <- cluster7$predicted.id
pred13 <- WhichCells(cluster7, idents = '13')

int13 <- intersect(sub13, pred13) #137

#cluster 14
Idents(cluster7) <- cluster7$seurat_clusters
sub14 <- WhichCells(cluster7, idents = c('16'))

Idents(cluster7) <- cluster7$predicted.id
pred14 <- WhichCells(cluster7, idents = '14')

int14 <- intersect(sub14, pred14) #92

#cluster 15
Idents(cluster7) <- cluster7$seurat_clusters
sub15 <- WhichCells(cluster7, idents = c('4'))

Idents(cluster7) <- cluster7$predicted.id
pred15 <- WhichCells(cluster7, idents = '15')

int15 <- intersect(sub15, pred15) #32

#total 
#pre:43892
#41875
total7 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total7,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query7_low_res_barcodes_10_23_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster7$query_label <- paste(cluster7$seurat_clusters, "query", sep = "_")
cluster7$ref_label <- paste(cluster7$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster7$ref_label), paste0(cluster7$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)


#-------Query 8--------#

cluster8 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub2_res20_4_low_res_cluster_pred.rds")

cluster8$predicted.id <- factor(cluster8$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster8, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster8, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table8 <- table(cluster8$predicted.id, cluster8$seurat_clusters)
table8 <- as.data.frame.matrix(table8)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table8_prop <- as.data.frame.matrix(apply(table8, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster8) <- cluster8$seurat_clusters
sub0 <- WhichCells(cluster8, idents = c('0'))

Idents(cluster8) <- cluster8$predicted.id
pred0 <- WhichCells(cluster8, idents = '0')

int0 <- intersect(sub0, pred0) #9296

#cluster 1
Idents(cluster8) <- cluster8$seurat_clusters
sub1 <- WhichCells(cluster8, idents = c('1', '2'))

Idents(cluster8) <- cluster8$predicted.id
pred1 <- WhichCells(cluster8, idents = '1')

int1 <- intersect(sub1, pred1) #8256

#cluster 2
Idents(cluster8) <- cluster8$seurat_clusters
sub2 <- WhichCells(cluster8, idents = c('2', '8'))

Idents(cluster8) <- cluster8$predicted.id
pred2 <- WhichCells(cluster8, idents = '2')

int2 <- intersect(sub2, pred2) #6247

#cluster 3
Idents(cluster8) <- cluster8$seurat_clusters
sub3 <- WhichCells(cluster8, idents = c('3'))

Idents(cluster8) <- cluster8$predicted.id
pred3 <- WhichCells(cluster8, idents = '3')

int3 <- intersect(sub3, pred3) #4235

#cluster 4
Idents(cluster8) <- cluster8$seurat_clusters
sub4 <- WhichCells(cluster8, idents = c('4'))

Idents(cluster8) <- cluster8$predicted.id
pred4 <- WhichCells(cluster8, idents = '4')

int4 <- intersect(sub4, pred4) #3206

#cluster 5
Idents(cluster8) <- cluster8$seurat_clusters
sub5 <- WhichCells(cluster8, idents = c('5'))

Idents(cluster8) <- cluster8$predicted.id
pred5 <- WhichCells(cluster8, idents = '5')

int5 <- intersect(sub5, pred5) #2596

#cluster 6
Idents(cluster8) <- cluster8$seurat_clusters
sub6 <- WhichCells(cluster8, idents = c('6'))

Idents(cluster8) <- cluster8$predicted.id
pred6 <- WhichCells(cluster8, idents = '6')

int6 <- intersect(sub6, pred6) #2463

#cluster 7
Idents(cluster8) <- cluster8$seurat_clusters
sub7 <- WhichCells(cluster8, idents = c('7'))

Idents(cluster8) <- cluster8$predicted.id
pred7 <- WhichCells(cluster8, idents = '7')

int7 <- intersect(sub7, pred7) #2122

#cluster 8
Idents(cluster8) <- cluster8$seurat_clusters
sub8 <- WhichCells(cluster8, idents = c('9'))

Idents(cluster8) <- cluster8$predicted.id
pred8 <- WhichCells(cluster8, idents = '8')

int8 <- intersect(sub8, pred8) #1463

#cluster 9 
Idents(cluster8) <- cluster8$seurat_clusters
sub9 <- WhichCells(cluster8, idents = c('10'))

Idents(cluster8) <- cluster8$predicted.id
pred9 <- WhichCells(cluster8, idents = '9')

int9 <- intersect(sub9, pred9) #935

#cluster 10 
Idents(cluster8) <- cluster8$seurat_clusters
sub10 <- WhichCells(cluster8, idents = c('11'))

Idents(cluster8) <- cluster8$predicted.id
pred10 <- WhichCells(cluster8, idents = '10')

int10 <- intersect(sub10, pred10) #800

#cluster 11
Idents(cluster8) <- cluster8$seurat_clusters
sub11 <- WhichCells(cluster8, idents = c('12'))

Idents(cluster8) <- cluster8$predicted.id
pred11 <- WhichCells(cluster8, idents = '11')

int11 <- intersect(sub11, pred11) #309

#cluster 12
Idents(cluster8) <- cluster8$seurat_clusters
sub12 <- WhichCells(cluster8, idents = c('13'))

Idents(cluster8) <- cluster8$predicted.id
pred12 <- WhichCells(cluster8, idents = '12')

int12 <- intersect(sub12, pred12) #190

#cluster 13
Idents(cluster8) <- cluster8$seurat_clusters
sub13 <- WhichCells(cluster8, idents = c('14'))

Idents(cluster8) <- cluster8$predicted.id
pred13 <- WhichCells(cluster8, idents = '13')

int13 <- intersect(sub13, pred13) #186

#cluster 14
Idents(cluster8) <- cluster8$seurat_clusters
sub14 <- WhichCells(cluster8, idents = c('15'))

Idents(cluster8) <- cluster8$predicted.id
pred14 <- WhichCells(cluster8, idents = '14')

int14 <- intersect(sub14, pred14) #83

#cluster 15
Idents(cluster8) <- cluster8$seurat_clusters
sub15 <- WhichCells(cluster8, idents = c('4'))

Idents(cluster8) <- cluster8$predicted.id
pred15 <- WhichCells(cluster8, idents = '15')

int15 <- intersect(sub15, pred15) #47

#total 
#pre:43893
#42434
total8 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total8,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query8_low_res_barcodes_10_23_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster8$query_label <- paste(cluster8$seurat_clusters, "query", sep = "_")
cluster8$ref_label <- paste(cluster8$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster8$ref_label), paste0(cluster8$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-------Query 9--------#

cluster9 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub_res20_5_low_res_cluster_pred.rds")

cluster9$predicted.id <- factor(cluster9$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster9, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster9, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table9 <- table(cluster9$predicted.id, cluster9$seurat_clusters)
table9 <- as.data.frame.matrix(table9)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table9_prop <- as.data.frame.matrix(apply(table9, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster9) <- cluster9$seurat_clusters
sub0 <- WhichCells(cluster9, idents = c('0'))

Idents(cluster9) <- cluster9$predicted.id
pred0 <- WhichCells(cluster9, idents = '0')

int0 <- intersect(sub0, pred0) #9419

#cluster 1
Idents(cluster9) <- cluster9$seurat_clusters
sub1 <- WhichCells(cluster9, idents = c('1'))

Idents(cluster9) <- cluster9$predicted.id
pred1 <- WhichCells(cluster9, idents = '1')

int1 <- intersect(sub1, pred1) #7804

#cluster 2 
Idents(cluster9) <- cluster9$seurat_clusters
sub2 <- WhichCells(cluster9, idents = c('2'))

Idents(cluster9) <- cluster9$predicted.id
pred2 <- WhichCells(cluster9, idents = '2')

int2 <- intersect(sub2, pred2) #6060

#cluster 3
Idents(cluster9) <- cluster9$seurat_clusters
sub3 <- WhichCells(cluster9, idents = c('3'))

Idents(cluster9) <- cluster9$predicted.id
pred3 <- WhichCells(cluster9, idents = '3')

int3 <- intersect(sub3, pred3) #4395

#cluster 4
Idents(cluster9) <- cluster9$seurat_clusters
sub4 <- WhichCells(cluster9, idents = c('4'))

Idents(cluster9) <- cluster9$predicted.id
pred4 <- WhichCells(cluster9, idents = '4')

int4 <- intersect(sub4, pred4) #3182

#cluster 5
Idents(cluster9) <- cluster9$seurat_clusters
sub5 <- WhichCells(cluster9, idents = c('5'))

Idents(cluster9) <- cluster9$predicted.id
pred5 <- WhichCells(cluster9, idents = '5')

int5 <- intersect(sub5, pred5) #2569

#cluster 6
Idents(cluster9) <- cluster9$seurat_clusters
sub6 <- WhichCells(cluster9, idents = c('6'))

Idents(cluster9) <- cluster9$predicted.id
pred6 <- WhichCells(cluster9, idents = '6')

int6 <- intersect(sub6, pred6) #2449

#cluster 7
Idents(cluster9) <- cluster9$seurat_clusters
sub7 <- WhichCells(cluster9, idents = c('7'))

Idents(cluster9) <- cluster9$predicted.id
pred7 <- WhichCells(cluster9, idents = '7')

int7 <- intersect(sub7, pred7) #2214

#cluster 8
Idents(cluster9) <- cluster9$seurat_clusters
sub8 <- WhichCells(cluster9, idents = c('8'))

Idents(cluster9) <- cluster9$predicted.id
pred8 <- WhichCells(cluster9, idents = '8')

int8 <- intersect(sub8, pred8) #1413

#cluster 9 
Idents(cluster9) <- cluster9$seurat_clusters
sub9 <- WhichCells(cluster9, idents = c('9'))

Idents(cluster9) <- cluster9$predicted.id
pred9 <- WhichCells(cluster9, idents = '9')

int9 <- intersect(sub9, pred9) #1127

#cluster 10
Idents(cluster9) <- cluster9$seurat_clusters
sub10 <- WhichCells(cluster9, idents = c('10'))

Idents(cluster9) <- cluster9$predicted.id
pred10 <- WhichCells(cluster9, idents = '10')

int10 <- intersect(sub10, pred10) #792

#cluster 11
Idents(cluster9) <- cluster9$seurat_clusters
sub11 <- WhichCells(cluster9, idents = c('11'))

Idents(cluster9) <- cluster9$predicted.id
pred11 <- WhichCells(cluster9, idents = '11')

int11 <- intersect(sub11, pred11) #309

#cluster 12
Idents(cluster9) <- cluster9$seurat_clusters
sub12 <- WhichCells(cluster9, idents = c('13'))

Idents(cluster9) <- cluster9$predicted.id
pred12 <- WhichCells(cluster9, idents = '12')

int12 <- intersect(sub12, pred12) #188

#cluster 13
Idents(cluster9) <- cluster9$seurat_clusters
sub13 <- WhichCells(cluster9, idents = c('12'))

Idents(cluster9) <- cluster9$predicted.id
pred13 <- WhichCells(cluster9, idents = '13')

int13 <- intersect(sub13, pred13) #192

#cluster 14
Idents(cluster9) <- cluster9$seurat_clusters
sub14 <- WhichCells(cluster9, idents = c('14'))

Idents(cluster9) <- cluster9$predicted.id
pred14 <- WhichCells(cluster9, idents = '14')

int14 <- intersect(sub14, pred14) #83

#cluster 15
Idents(cluster9) <- cluster9$seurat_clusters
sub15 <- WhichCells(cluster9, idents = c('4'))

Idents(cluster9) <- cluster9$predicted.id
pred15 <- WhichCells(cluster9, idents = '15')

int15 <- intersect(sub15, pred15) #36

#total 
#pre:43892
#42232
total9 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total9,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query9_low_res_barcodes_10_23_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster9$query_label <- paste(cluster9$seurat_clusters, "query", sep = "_")
cluster9$ref_label <- paste(cluster9$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster9$ref_label), paste0(cluster9$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-------Query 10--------#

cluster10 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub2_res20_5_low_res_cluster_pred.rds")

cluster10$predicted.id <- factor(cluster10$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster10, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster10, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table10 <- table(cluster10$predicted.id, cluster10$seurat_clusters)
table10 <- as.data.frame.matrix(table10)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table10_prop <- as.data.frame.matrix(apply(table10, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster10) <- cluster10$seurat_clusters
sub0 <- WhichCells(cluster10, idents = c('1'))

Idents(cluster10) <- cluster10$predicted.id
pred0 <- WhichCells(cluster10, idents = '0')

int0 <- intersect(sub0, pred0) #9090

#cluster 1
Idents(cluster10) <- cluster10$seurat_clusters
sub1 <- WhichCells(cluster10, idents = c('0'))

Idents(cluster10) <- cluster10$predicted.id
pred1 <- WhichCells(cluster10, idents = '1')

int1 <- intersect(sub1, pred1) #7896

#cluster 2
Idents(cluster10) <- cluster10$seurat_clusters
sub2 <- WhichCells(cluster10, idents = c('2', '8'))

Idents(cluster10) <- cluster10$predicted.id
pred2 <- WhichCells(cluster10, idents = '2')

int2 <- intersect(sub2, pred2) #6324

#cluster 3
Idents(cluster10) <- cluster10$seurat_clusters
sub3 <- WhichCells(cluster10, idents = c('3'))

Idents(cluster10) <- cluster10$predicted.id
pred3 <- WhichCells(cluster10, idents = '3')

int3 <- intersect(sub3, pred3) #4296

#cluster 4
Idents(cluster10) <- cluster10$seurat_clusters
sub4 <- WhichCells(cluster10, idents = c('4'))

Idents(cluster10) <- cluster10$predicted.id
pred4 <- WhichCells(cluster10, idents = '4')

int4 <- intersect(sub4, pred4) #3186

#cluster 5
Idents(cluster10) <- cluster10$seurat_clusters
sub5 <- WhichCells(cluster10, idents = c('6'))

Idents(cluster10) <- cluster10$predicted.id
pred5 <- WhichCells(cluster10, idents = '5')

int5 <- intersect(sub5, pred5) #2584

#cluster 6
Idents(cluster10) <- cluster10$seurat_clusters
sub6 <- WhichCells(cluster10, idents = c('5'))

Idents(cluster10) <- cluster10$predicted.id
pred6 <- WhichCells(cluster10, idents = '6')

int6 <- intersect(sub6, pred6) #2509

#cluster 7
Idents(cluster10) <- cluster10$seurat_clusters
sub7 <- WhichCells(cluster10, idents = c('7'))

Idents(cluster10) <- cluster10$predicted.id
pred7 <- WhichCells(cluster10, idents = '7')

int7 <- intersect(sub7, pred7) #2158

#cluster 8 
Idents(cluster10) <- cluster10$seurat_clusters
sub8 <- WhichCells(cluster10, idents = c('9'))

Idents(cluster10) <- cluster10$predicted.id
pred8 <- WhichCells(cluster10, idents = '8')

int8 <- intersect(sub8, pred8) #1443

#cluster 9 
Idents(cluster10) <- cluster10$seurat_clusters
sub9 <- WhichCells(cluster10, idents = c('10'))

Idents(cluster10) <- cluster10$predicted.id
pred9 <- WhichCells(cluster10, idents = '9')

int9 <- intersect(sub9, pred9) #979

#cluster 10 
Idents(cluster10) <- cluster10$seurat_clusters
sub10 <- WhichCells(cluster10, idents = c('11'))

Idents(cluster10) <- cluster10$predicted.id
pred10 <- WhichCells(cluster10, idents = '10')

int10 <- intersect(sub10, pred10) #805

#cluster 11
Idents(cluster10) <- cluster10$seurat_clusters
sub11 <- WhichCells(cluster10, idents = c('12'))

Idents(cluster10) <- cluster10$predicted.id
pred11 <- WhichCells(cluster10, idents = '11')

int11 <- intersect(sub11, pred11) #323

#cluster 12
Idents(cluster10) <- cluster10$seurat_clusters
sub12 <- WhichCells(cluster10, idents = c('13'))

Idents(cluster10) <- cluster10$predicted.id
pred12 <- WhichCells(cluster10, idents = '12')

int12 <- intersect(sub12, pred12) #205

#cluster 13
Idents(cluster10) <- cluster10$seurat_clusters
sub13 <- WhichCells(cluster10, idents = c('14'))

Idents(cluster10) <- cluster10$predicted.id
pred13 <- WhichCells(cluster10, idents = '13')

int13 <- intersect(sub13, pred13) #164

#cluster 14
Idents(cluster10) <- cluster10$seurat_clusters
sub14 <- WhichCells(cluster10, idents = c('15'))

Idents(cluster10) <- cluster10$predicted.id
pred14 <- WhichCells(cluster10, idents = '14')

int14 <- intersect(sub14, pred14) #83

#cluster 15
Idents(cluster10) <- cluster10$seurat_clusters
sub15 <- WhichCells(cluster10, idents = c('10'))

Idents(cluster10) <- cluster10$predicted.id
pred15 <- WhichCells(cluster10, idents = '15')

int15 <- intersect(sub15, pred15) #43

#total 
#pre:43893
#42088
total10 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total10,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query10_low_res_barcodes_10_23_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster10$query_label <- paste(cluster10$seurat_clusters, "query", sep = "_")
cluster10$ref_label <- paste(cluster10$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster10$ref_label), paste0(cluster10$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-----load in the stable barcodes----------#
#load in barcodes 
total1 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query1_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)
total2 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query2_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)
total3 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query3_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)
total4 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query4_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)
total5 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query5_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)
total6 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query6_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)
total7 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query7_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)
total8 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query8_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)
total9 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query9_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)
total10 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/query10_low_res_barcodes_10_23_24.csv', what = "", sep = ",", skip = 1)

#merge the barcodes into a list 
#423009 barcodes 
total_merge <- c(total1, total2, total3, total4, total5, total6, total7, total8, total9, total10)

#counts the number of times each barcode is present
value_counts <- table(total_merge)
table(value_counts)

#1: 837
#2: 1108
#3: 1791
#4: 3927
#5: 79775

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
write.csv(filtered_values, file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/filtered_low_res_barcodes_10_24_24.csv', row.names = F)

#subset the ileum seurat object on the server 

#load SO
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_low_res_ref_10_21_24.rds")

#load in stably assigned cells 
filtered_values <- scan('/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/filtered_low_res_barcodes_10_24_24.csv', what = "", sep = ",", skip = 1)

#subset ileum based on stable barcodes 
ileum_sub <- subset(ileum, cells = filtered_values)

#save 
saveRDS(ileum_sub, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_low_res_stable_cell_10_14_24.rds")

#--------Visualize the unstably assigned cells----------#

#load in ileum SO 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_rPCA_10_17_24.rds")

#ileum barcodes 
barcodes <- colnames(ileum)

unstable <- setdiff(barcodes, filtered_values)

DimPlot(ileum, reduction = "umap", cells.highlight = unstable, sizes.highlight = 0.1) + scale_color_manual(labels = c("Stably Assigned Cells", "Unstably Assigned Cells"), values = c("grey", "red"))


#----------ARI-----------#

#10/25/2024

#merged ARI scores
df <- read.csv('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/df_ri_helm_batch1_13_ileum_low_res_dune_merge_10_25_24.csv', row.names = NULL, sep = ",")

df <- df[,-1]
names(df) <- gsub("^paste0\\.", "", names(df))

columns <- colnames(df)

#plot the ARI - without paste0
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/figures/RI_helm_batch1_13_ileum_low_res_all_param_merge_10_25_24.pdf", width = 15, height = 15)
plot(plotARIs(df %>% select(columns)) + theme(text=element_text(size=8)) + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) + theme(axis.text.y = element_text(size = 10)))
dev.off()


#----mean and median----#
#ARI values 
plot <- plotARIs(clusMat = df)

ari_values <- data.frame(plot[["plot_env"]][["Mat"]])

#calculate mean ARI
ari_values$mean_ari <- rowMeans(ari_values)

#top 10 mean ari 
mean_ari <- as.data.frame(ari_values$mean_ari)
rownames(mean_ari) <- rownames(ari_values)
top_n(mean_ari, 10)

#overall mean ari
summarise(ari_values, overall_mean = mean(mean_ari)) #0.8658195

#remove mean_ari from ari_values 
ari_values <- ari_values[, -c(37)]

#median ARI

ari_values$row_median = apply(ari_values, 1, median, na.rm=TRUE)

overall_median <- median(ari_values$row_median)  #0.8612033


#take average of median 
summarise(ari_values, overall_mean = mean(row_median)) #0.7654572

#top 10 median ARI scores 
median_ari <- as.data.frame(ari_values$row_median)
rownames(median_ari) <- rownames(ari_values)
top_n(median_ari, 10)


#--------------High Resolution Stability Assessment--------------#

#10/25/2024

#-------Query 1---------#

cluster1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub_res15_1_high_res_cluster_pred.rds")

cluster1$predicted.id <- factor(cluster1$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster1, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster1, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

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

int0 <- intersect(sub0, pred0) #4760

#cluster 1
Idents(cluster1) <- cluster1$seurat_clusters
sub1 <- WhichCells(cluster1, idents = c('4', '7'))

Idents(cluster1) <- cluster1$predicted.id
pred1 <- WhichCells(cluster1, idents = '1')

int1 <- intersect(sub1, pred1) #3254

#cluster 2
Idents(cluster1) <- cluster1$seurat_clusters
sub2 <- WhichCells(cluster1, idents = c('1'))

Idents(cluster1) <- cluster1$predicted.id
pred2 <- WhichCells(cluster1, idents = '2')

int2 <- intersect(sub2, pred2) #2395

#cluster 3
Idents(cluster1) <- cluster1$seurat_clusters
sub3 <- WhichCells(cluster1, idents = c('5', '18'))

Idents(cluster1) <- cluster1$predicted.id
pred3 <- WhichCells(cluster1, idents = '3')

int3 <- intersect(sub3, pred3) #2182

#cluster 4
Idents(cluster1) <- cluster1$seurat_clusters
sub4 <- WhichCells(cluster1, idents = c('3'))

Idents(cluster1) <- cluster1$predicted.id
pred4 <- WhichCells(cluster1, idents = '4')

int4 <- intersect(sub4, pred4) #2109

#cluster 5
Idents(cluster1) <- cluster1$seurat_clusters
sub5 <- WhichCells(cluster1, idents = c('2'))

Idents(cluster1) <- cluster1$predicted.id
pred5 <- WhichCells(cluster1, idents = '5')

int5 <- intersect(sub5, pred5) #2190

#cluster 6
Idents(cluster1) <- cluster1$seurat_clusters
sub6 <- WhichCells(cluster1, idents = c('1', '8'))

Idents(cluster1) <- cluster1$predicted.id
pred6 <- WhichCells(cluster1, idents = '6')

int6 <- intersect(sub6, pred6) #1928

#cluster 7
Idents(cluster1) <- cluster1$seurat_clusters
sub7 <- WhichCells(cluster1, idents = c('5', '6'))

Idents(cluster1) <- cluster1$predicted.id
pred7 <- WhichCells(cluster1, idents = '7')

int7 <- intersect(sub7, pred7) #1857

#cluster 8 
Idents(cluster1) <- cluster1$seurat_clusters
sub8 <- WhichCells(cluster1, idents = c('15', '28'))

Idents(cluster1) <- cluster1$predicted.id
pred8 <- WhichCells(cluster1, idents = '8')

int8 <- intersect(sub8, pred8) #1753

#cluster 9 
Idents(cluster1) <- cluster1$seurat_clusters
sub9 <- WhichCells(cluster1, idents = c('12', '27'))

Idents(cluster1) <- cluster1$predicted.id
pred9 <- WhichCells(cluster1, idents = '9')

int9 <- intersect(sub9, pred9) #1889

#cluster 10 
Idents(cluster1) <- cluster1$seurat_clusters
sub10 <- WhichCells(cluster1, idents = c('9'))

Idents(cluster1) <- cluster1$predicted.id
pred10 <- WhichCells(cluster1, idents = '10')

int10 <- intersect(sub10, pred10) #1411

#cluster 11
Idents(cluster1) <- cluster1$seurat_clusters
sub11 <- WhichCells(cluster1, idents = c('10'))

Idents(cluster1) <- cluster1$predicted.id
pred11 <- WhichCells(cluster1, idents = '11')

int11 <- intersect(sub11, pred11) #1105

#cluster 12
Idents(cluster1) <- cluster1$seurat_clusters
sub12 <- WhichCells(cluster1, idents = c('14', '19'))

Idents(cluster1) <- cluster1$predicted.id
pred12 <- WhichCells(cluster1, idents = '12')

int12 <- intersect(sub12, pred12) #887

#cluster 13
Idents(cluster1) <- cluster1$seurat_clusters
sub13 <- WhichCells(cluster1, idents = c('14'))

Idents(cluster1) <- cluster1$predicted.id
pred13 <- WhichCells(cluster1, idents = '13')

int13 <- intersect(sub13, pred13) #1168

#cluster 14
Idents(cluster1) <- cluster1$seurat_clusters
sub14 <- WhichCells(cluster1, idents = c('13'))

Idents(cluster1) <- cluster1$predicted.id
pred14 <- WhichCells(cluster1, idents = '14')

int14 <- intersect(sub14, pred14) #1184

#cluster 15
Idents(cluster1) <- cluster1$seurat_clusters
sub15 <- WhichCells(cluster1, idents = c('21', '25'))

Idents(cluster1) <- cluster1$predicted.id
pred15 <- WhichCells(cluster1, idents = '15')

int15 <- intersect(sub15, pred15) #1190

#cluster 16
Idents(cluster1) <- cluster1$seurat_clusters
sub16 <- WhichCells(cluster1, idents = c('16', '22'))

Idents(cluster1) <- cluster1$predicted.id
pred16 <- WhichCells(cluster1, idents = '16')

int16 <- intersect(sub16, pred16) #882

#cluster 17
Idents(cluster1) <- cluster1$seurat_clusters
sub17 <- WhichCells(cluster1, idents = c('11'))

Idents(cluster1) <- cluster1$predicted.id
pred17 <- WhichCells(cluster1, idents = '17')

int17 <- intersect(sub17, pred17) #1004

#cluster 18 
Idents(cluster1) <- cluster1$seurat_clusters
sub18 <- WhichCells(cluster1, idents = c('20', '25'))

Idents(cluster1) <- cluster1$predicted.id
pred18 <- WhichCells(cluster1, idents = '18')

int18 <- intersect(sub18, pred18) #979

#cluster 19 
Idents(cluster1) <- cluster1$seurat_clusters
sub19 <- WhichCells(cluster1, idents = c('2', '24'))

Idents(cluster1) <- cluster1$predicted.id
pred19 <- WhichCells(cluster1, idents = '19')

int19 <- intersect(sub19, pred19) #814

#cluster 20 
Idents(cluster1) <- cluster1$seurat_clusters
sub20 <- WhichCells(cluster1, idents = c('16'))

Idents(cluster1) <- cluster1$predicted.id
pred20 <- WhichCells(cluster1, idents = '20')

int20 <- intersect(sub20, pred20) #954

#cluster 21
Idents(cluster1) <- cluster1$seurat_clusters
sub21 <- WhichCells(cluster1, idents = c('23'))

Idents(cluster1) <- cluster1$predicted.id
pred21 <- WhichCells(cluster1, idents = '21')

int21 <- intersect(sub21, pred21) #629

#cluster 22
Idents(cluster1) <- cluster1$seurat_clusters
sub22 <- WhichCells(cluster1, idents = c('26'))

Idents(cluster1) <- cluster1$predicted.id
pred22 <- WhichCells(cluster1, idents = '22')

int22 <- intersect(sub22, pred22) #575

#cluster 23
Idents(cluster1) <- cluster1$seurat_clusters
sub23 <- WhichCells(cluster1, idents = c('31'))

Idents(cluster1) <- cluster1$predicted.id
pred23 <- WhichCells(cluster1, idents = '23')

int23 <- intersect(sub23, pred23) #408

#cluster 24
Idents(cluster1) <- cluster1$seurat_clusters
sub24 <- WhichCells(cluster1, idents = c('29'))

Idents(cluster1) <- cluster1$predicted.id
pred24 <- WhichCells(cluster1, idents = '24')

int24 <- intersect(sub24, pred24) #448

#cluster 25
Idents(cluster1) <- cluster1$seurat_clusters
sub25 <- WhichCells(cluster1, idents = c('17'))

Idents(cluster1) <- cluster1$predicted.id
pred25 <- WhichCells(cluster1, idents = '25')

int25 <- intersect(sub25, pred25) #309

#cluster 26
Idents(cluster1) <- cluster1$seurat_clusters
sub26 <- WhichCells(cluster1, idents = c('30'))

Idents(cluster1) <- cluster1$predicted.id
pred26 <- WhichCells(cluster1, idents = '26')

int26 <- intersect(sub26, pred26) #352

#cluster 27
Idents(cluster1) <- cluster1$seurat_clusters
sub27 <- WhichCells(cluster1, idents = c('32'))

Idents(cluster1) <- cluster1$predicted.id
pred27 <- WhichCells(cluster1, idents = '27')

int27 <- intersect(sub27, pred27) #294

#cluster 28
Idents(cluster1) <- cluster1$seurat_clusters
sub28 <- WhichCells(cluster1, idents = c('17'))

Idents(cluster1) <- cluster1$predicted.id
pred28 <- WhichCells(cluster1, idents = '28')

int28 <- intersect(sub28, pred28) #274

#cluster 29
Idents(cluster1) <- cluster1$seurat_clusters
sub29 <- WhichCells(cluster1, idents = c('33'))

Idents(cluster1) <- cluster1$predicted.id
pred29 <- WhichCells(cluster1, idents = '29')

int29 <- intersect(sub29, pred29) #221

#cluster 30 
Idents(cluster1) <- cluster1$seurat_clusters
sub30 <- WhichCells(cluster1, idents = c('34'))

Idents(cluster1) <- cluster1$predicted.id
pred30 <- WhichCells(cluster1, idents = '30')

int30 <- intersect(sub30, pred30) #180

#cluster 31
Idents(cluster1) <- cluster1$seurat_clusters
sub31 <- WhichCells(cluster1, idents = c('17', '36'))

Idents(cluster1) <- cluster1$predicted.id
pred31 <- WhichCells(cluster1, idents = '31')

int31 <- intersect(sub31, pred31) #92

#cluster 32
Idents(cluster1) <- cluster1$seurat_clusters
sub32 <- WhichCells(cluster1, idents = c('35'))

Idents(cluster1) <- cluster1$predicted.id
pred32 <- WhichCells(cluster1, idents = '32')

int32 <- intersect(sub32, pred32) #101

#cluster 33
Idents(cluster1) <- cluster1$seurat_clusters
sub33 <- WhichCells(cluster1, idents = c('18'))

Idents(cluster1) <- cluster1$predicted.id
pred33 <- WhichCells(cluster1, idents = '33')

int33 <- intersect(sub33, pred33) #79

#cluster 34
Idents(cluster1) <- cluster1$seurat_clusters
sub34 <- WhichCells(cluster1, idents = c('37'))

Idents(cluster1) <- cluster1$predicted.id
pred34 <- WhichCells(cluster1, idents = '34')

int34 <- intersect(sub34, pred34) #55

#cluster 35
Idents(cluster1) <- cluster1$seurat_clusters
sub35 <- WhichCells(cluster1, idents = c('17'))

Idents(cluster1) <- cluster1$predicted.id
pred35 <- WhichCells(cluster1, idents = '35')

int35 <- intersect(sub35, pred35) #34

#cluster 36
Idents(cluster1) <- cluster1$seurat_clusters
sub36 <- WhichCells(cluster1, idents = c('39'))

Idents(cluster1) <- cluster1$predicted.id
pred36 <- WhichCells(cluster1, idents = '36')

int36 <- intersect(sub36, pred36) #53

#cluster 37
Idents(cluster1) <- cluster1$seurat_clusters
sub37 <- WhichCells(cluster1, idents = c('38'))

Idents(cluster1) <- cluster1$predicted.id
pred37 <- WhichCells(cluster1, idents = '37')

int37 <- intersect(sub37, pred37) #54

#cluster 38
Idents(cluster1) <- cluster1$seurat_clusters
sub38 <- WhichCells(cluster1, idents = c('4', '16'))

Idents(cluster1) <- cluster1$predicted.id
pred38 <- WhichCells(cluster1, idents = '38')

int38 <- intersect(sub38, pred38) #16

#total 
#pre:43892
#40069
total1 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total1,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query1_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster1$query_label <- paste(cluster1$seurat_clusters, "query", sep = "_")
cluster1$ref_label <- paste(cluster1$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster1$ref_label), paste0(cluster1$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)


#-------Query 2---------#

cluster2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub2_res15_1_high_res_cluster_pred.rds")

cluster2$predicted.id <- factor(cluster2$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster2, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster2, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table2 <- table(cluster2$predicted.id, cluster2$seurat_clusters)
table2 <- as.data.frame.matrix(table2)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table2_prop <- as.data.frame.matrix(apply(table2, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster2) <- cluster2$seurat_clusters
sub0 <- WhichCells(cluster2, idents = c('0'))

Idents(cluster2) <- cluster2$predicted.id
pred0 <- WhichCells(cluster2, idents = '0')

int0 <- intersect(sub0, pred0) #4751

#cluster 1
Idents(cluster2) <- cluster2$seurat_clusters
sub1 <- WhichCells(cluster2, idents = c('1'))

Idents(cluster2) <- cluster2$predicted.id
pred1 <- WhichCells(cluster2, idents = '1')

int1 <- intersect(sub1, pred1) #3070

#cluster 2
Idents(cluster2) <- cluster2$seurat_clusters
sub2 <- WhichCells(cluster2, idents = c('2'))

Idents(cluster2) <- cluster2$predicted.id
pred2 <- WhichCells(cluster2, idents = '2')

int2 <- intersect(sub2, pred2) #2201

#cluster 3
Idents(cluster2) <- cluster2$seurat_clusters
sub3 <- WhichCells(cluster2, idents = c('3', '4'))

Idents(cluster2) <- cluster2$predicted.id
pred3 <- WhichCells(cluster2, idents = '3')

int3 <- intersect(sub3, pred3) #2045

#cluster 4
Idents(cluster2) <- cluster2$seurat_clusters
sub4 <- WhichCells(cluster2, idents = c('4'))

Idents(cluster2) <- cluster2$predicted.id
pred4 <- WhichCells(cluster2, idents = '4')

int4 <- intersect(sub4, pred4) #1986

#cluster 5
Idents(cluster2) <- cluster2$seurat_clusters
sub5 <- WhichCells(cluster2, idents = c('6', '29'))

Idents(cluster2) <- cluster2$predicted.id
pred5 <- WhichCells(cluster2, idents = '5')

int5 <- intersect(sub5, pred5) #2142

#cluster 6
Idents(cluster2) <- cluster2$seurat_clusters
sub6 <- WhichCells(cluster2, idents = c('7', '9', '25'))

Idents(cluster2) <- cluster2$predicted.id
pred6 <- WhichCells(cluster2, idents = '6')

int6 <- intersect(sub6, pred6) #2167

#cluster 7 
Idents(cluster2) <- cluster2$seurat_clusters
sub7 <- WhichCells(cluster2, idents = c('3', '16'))

Idents(cluster2) <- cluster2$predicted.id
pred7 <- WhichCells(cluster2, idents = '7')

int7 <- intersect(sub7, pred7) #1825

#cluster 8 
Idents(cluster2) <- cluster2$seurat_clusters
sub8 <- WhichCells(cluster2, idents = c('13', '22'))

Idents(cluster2) <- cluster2$predicted.id
pred8 <- WhichCells(cluster2, idents = '8')

int8 <- intersect(sub8, pred8) #1872

#cluster 9 
Idents(cluster2) <- cluster2$seurat_clusters
sub9 <- WhichCells(cluster2, idents = c('5'))

Idents(cluster2) <- cluster2$predicted.id
pred9 <- WhichCells(cluster2, idents = '9')

int9 <- intersect(sub9, pred9) #1901

#cluster 10 
Idents(cluster2) <- cluster2$seurat_clusters
sub10 <- WhichCells(cluster2, idents = c('10'))

Idents(cluster2) <- cluster2$predicted.id
pred10 <- WhichCells(cluster2, idents = '10')

int10 <- intersect(sub10, pred10) #1445

#cluster 11
Idents(cluster2) <- cluster2$seurat_clusters
sub11 <- WhichCells(cluster2, idents = c('9', '26'))

Idents(cluster2) <- cluster2$predicted.id
pred11 <- WhichCells(cluster2, idents = '11')

int11 <- intersect(sub11, pred11) #1120

#cluster 12
Idents(cluster2) <- cluster2$seurat_clusters
sub12 <- WhichCells(cluster2, idents = c('8'))

Idents(cluster2) <- cluster2$predicted.id
pred12 <- WhichCells(cluster2, idents = '12')

int12 <- intersect(sub12, pred12) #1116

#cluster 13
Idents(cluster2) <- cluster2$seurat_clusters
sub13 <- WhichCells(cluster2, idents = c('17', '26'))

Idents(cluster2) <- cluster2$predicted.id
pred13 <- WhichCells(cluster2, idents = '13')

int13 <- intersect(sub13, pred13) #1233

#cluster 14
Idents(cluster2) <- cluster2$seurat_clusters
sub14 <- WhichCells(cluster2, idents = c('12'))

Idents(cluster2) <- cluster2$predicted.id
pred14 <- WhichCells(cluster2, idents = '14')

int14 <- intersect(sub14, pred14) #1195

#cluster 15
Idents(cluster2) <- cluster2$seurat_clusters
sub15 <- WhichCells(cluster2, idents = c('15'))

Idents(cluster2) <- cluster2$predicted.id
pred15 <- WhichCells(cluster2, idents = '15')

int15 <- intersect(sub15, pred15) #1123

#cluster 16
Idents(cluster2) <- cluster2$seurat_clusters
sub16 <- WhichCells(cluster2, idents = c('18'))

Idents(cluster2) <- cluster2$predicted.id
pred16 <- WhichCells(cluster2, idents = '16')

int16 <- intersect(sub16, pred16) #940

#cluster 17 
Idents(cluster2) <- cluster2$seurat_clusters
sub17 <- WhichCells(cluster2, idents = c('11'))

Idents(cluster2) <- cluster2$predicted.id
pred17 <- WhichCells(cluster2, idents = '17')

int17 <- intersect(sub17, pred17) #1050

#cluster 18 
Idents(cluster2) <- cluster2$seurat_clusters
sub18 <- WhichCells(cluster2, idents = c('14'))

Idents(cluster2) <- cluster2$predicted.id
pred18 <- WhichCells(cluster2, idents = '18')

int18 <- intersect(sub18, pred18) #972

#cluster 19 
Idents(cluster2) <- cluster2$seurat_clusters
sub19 <- WhichCells(cluster2, idents = c('20'))

Idents(cluster2) <- cluster2$predicted.id
pred19 <- WhichCells(cluster2, idents = '19')

int19 <- intersect(sub19, pred19) #714

#cluster 20 
Idents(cluster2) <- cluster2$seurat_clusters
sub20 <- WhichCells(cluster2, idents = c('19'))

Idents(cluster2) <- cluster2$predicted.id
pred20 <- WhichCells(cluster2, idents = '20')

int20 <- intersect(sub20, pred20) #910

#cluster 21
Idents(cluster2) <- cluster2$seurat_clusters
sub21 <- WhichCells(cluster2, idents = c('21'))

Idents(cluster2) <- cluster2$predicted.id
pred21 <- WhichCells(cluster2, idents = '21')

int21 <- intersect(sub21, pred21) #635

#cluster 22
Idents(cluster2) <- cluster2$seurat_clusters
sub22 <- WhichCells(cluster2, idents = c('23'))

Idents(cluster2) <- cluster2$predicted.id
pred22 <- WhichCells(cluster2, idents = '22')

int22 <- intersect(sub22, pred22) #554

#cluster 23
Idents(cluster2) <- cluster2$seurat_clusters
sub23 <- WhichCells(cluster2, idents = c('28', '31'))

Idents(cluster2) <- cluster2$predicted.id
pred23 <- WhichCells(cluster2, idents = '23')

int23 <- intersect(sub23, pred23) #425

#cluster 24
Idents(cluster2) <- cluster2$seurat_clusters
sub24 <- WhichCells(cluster2, idents = c('27'))

Idents(cluster2) <- cluster2$predicted.id
pred24 <- WhichCells(cluster2, idents = '24')

int24 <- intersect(sub24, pred24) #466

#cluster 25
Idents(cluster2) <- cluster2$seurat_clusters
sub25 <- WhichCells(cluster2, idents = c('32'))

Idents(cluster2) <- cluster2$predicted.id
pred25 <- WhichCells(cluster2, idents = '25')

int25 <- intersect(sub25, pred25) #310

#cluster 26
Idents(cluster2) <- cluster2$seurat_clusters
sub26 <- WhichCells(cluster2, idents = c('28'))

Idents(cluster2) <- cluster2$predicted.id
pred26 <- WhichCells(cluster2, idents = '26')

int26 <- intersect(sub26, pred26) #344

#cluster 27
Idents(cluster2) <- cluster2$seurat_clusters
sub27 <- WhichCells(cluster2, idents = c('30'))

Idents(cluster2) <- cluster2$predicted.id
pred27 <- WhichCells(cluster2, idents = '27')

int27 <- intersect(sub27, pred27) #352

#cluster 28
Idents(cluster2) <- cluster2$seurat_clusters
sub28 <- WhichCells(cluster2, idents = c('25'))

Idents(cluster2) <- cluster2$predicted.id
pred28 <- WhichCells(cluster2, idents = '28')

int28 <- intersect(sub28, pred28) #250

#cluster 29 
Idents(cluster2) <- cluster2$seurat_clusters
sub29 <- WhichCells(cluster2, idents = c('34'))

Idents(cluster2) <- cluster2$predicted.id
pred29 <- WhichCells(cluster2, idents = '29')

int29 <- intersect(sub29, pred29) #173

#cluster 30 
Idents(cluster2) <- cluster2$seurat_clusters
sub30 <- WhichCells(cluster2, idents = c('33'))

Idents(cluster2) <- cluster2$predicted.id
pred30 <- WhichCells(cluster2, idents = '30')

int30 <- intersect(sub30, pred30) #211

#cluster 31
Idents(cluster2) <- cluster2$seurat_clusters
sub31 <- WhichCells(cluster2, idents = c('9', '35'))

Idents(cluster2) <- cluster2$predicted.id
pred31 <- WhichCells(cluster2, idents = '31')

int31 <- intersect(sub31, pred31) #121

#cluster 32
Idents(cluster2) <- cluster2$seurat_clusters
sub32 <- WhichCells(cluster2, idents = c('36'))

Idents(cluster2) <- cluster2$predicted.id
pred32 <- WhichCells(cluster2, idents = '32')

int32 <- intersect(sub32, pred32) #72

#cluster 33
Idents(cluster2) <- cluster2$seurat_clusters
sub33 <- WhichCells(cluster2, idents = c('37'))

Idents(cluster2) <- cluster2$predicted.id
pred33 <- WhichCells(cluster2, idents = '33')

int33 <- intersect(sub33, pred33) #68

#cluster 34
Idents(cluster2) <- cluster2$seurat_clusters
sub34 <- WhichCells(cluster2, idents = c('5'))

Idents(cluster2) <- cluster2$predicted.id
pred34 <- WhichCells(cluster2, idents = '34')

int34 <- intersect(sub34, pred34) #49

#cluster 35
Idents(cluster2) <- cluster2$seurat_clusters
sub35 <- WhichCells(cluster2, idents = c('40'))

Idents(cluster2) <- cluster2$predicted.id
pred35 <- WhichCells(cluster2, idents = '35')

int35 <- intersect(sub35, pred35) #40

#cluster 36
Idents(cluster2) <- cluster2$seurat_clusters
sub36 <- WhichCells(cluster2, idents = c('39'))

Idents(cluster2) <- cluster2$predicted.id
pred36 <- WhichCells(cluster2, idents = '36')

int36 <- intersect(sub36, pred36) #44

#cluster 37
Idents(cluster2) <- cluster2$seurat_clusters
sub37 <- WhichCells(cluster2, idents = c('38'))

Idents(cluster2) <- cluster2$predicted.id
pred37 <- WhichCells(cluster2, idents = '37')

int37 <- intersect(sub37, pred37) #38

#cluster 38 
Idents(cluster2) <- cluster2$seurat_clusters
sub38 <- WhichCells(cluster2, idents = c('18'))

Idents(cluster2) <- cluster2$predicted.id
pred38 <- WhichCells(cluster2, idents = '38')

int38 <- intersect(sub38, pred38) #27

#total 
#pre:43893
#39957
total2 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total2,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query2_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster2$query_label <- paste(cluster2$seurat_clusters, "query", sep = "_")
cluster2$ref_label <- paste(cluster2$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster2$ref_label), paste0(cluster2$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-------Query 3---------#

cluster3 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub_res15_2_high_res_cluster_pred.rds")

cluster3$predicted.id <- factor(cluster3$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster3, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster3, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table3 <- table(cluster3$predicted.id, cluster3$seurat_clusters)
table3 <- as.data.frame.matrix(table3)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table3_prop <- as.data.frame.matrix(apply(table3, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster3) <- cluster3$seurat_clusters
sub0 <- WhichCells(cluster3, idents = c('0', '15'))

Idents(cluster3) <- cluster3$predicted.id
pred0 <- WhichCells(cluster3, idents = '0')

int0 <- intersect(sub0, pred0) #4751

#cluster 1
Idents(cluster3) <- cluster3$seurat_clusters
sub1 <- WhichCells(cluster3, idents = c('2', '5'))

Idents(cluster3) <- cluster3$predicted.id
pred1 <- WhichCells(cluster3, idents = '1')

int1 <- intersect(sub1, pred1) #3667

#cluster 2
Idents(cluster3) <- cluster3$seurat_clusters
sub2 <- WhichCells(cluster3, idents = c('1'))

Idents(cluster3) <- cluster3$predicted.id
pred2 <- WhichCells(cluster3, idents = '2')

int2 <- intersect(sub2, pred2) #2402

#cluster 3
Idents(cluster3) <- cluster3$seurat_clusters
sub3 <- WhichCells(cluster3, idents = c('6', '16'))

Idents(cluster3) <- cluster3$predicted.id
pred3 <- WhichCells(cluster3, idents = '3')

int3 <- intersect(sub3, pred3) #2014

#cluster 4
Idents(cluster3) <- cluster3$seurat_clusters
sub4 <- WhichCells(cluster3, idents = c('16', '19', '21'))

Idents(cluster3) <- cluster3$predicted.id
pred4 <- WhichCells(cluster3, idents = '4')

int4 <- intersect(sub4, pred4) #2345

#cluster 5
Idents(cluster3) <- cluster3$seurat_clusters
sub5 <- WhichCells(cluster3, idents = c('3'))

Idents(cluster3) <- cluster3$predicted.id
pred5 <- WhichCells(cluster3, idents = '5')

int5 <- intersect(sub5, pred5) #2195

#cluster 6
Idents(cluster3) <- cluster3$seurat_clusters
sub6 <- WhichCells(cluster3, idents = c('4','13'))

Idents(cluster3) <- cluster3$predicted.id
pred6 <- WhichCells(cluster3, idents = '6')

int6 <- intersect(sub6, pred6) #1843

#cluster 7
Idents(cluster3) <- cluster3$seurat_clusters
sub7 <- WhichCells(cluster3, idents = c('1', '6', '10'))

Idents(cluster3) <- cluster3$predicted.id
pred7 <- WhichCells(cluster3, idents = '7')

int7 <- intersect(sub7, pred7) #2074

#cluster 8 
Idents(cluster3) <- cluster3$seurat_clusters
sub8 <- WhichCells(cluster3, idents = c('9', '25'))

Idents(cluster3) <- cluster3$predicted.id
pred8 <- WhichCells(cluster3, idents = '8')

int8 <- intersect(sub8, pred8) #1961

#cluster 9 
Idents(cluster3) <- cluster3$seurat_clusters
sub9 <- WhichCells(cluster3, idents = c('7'))

Idents(cluster3) <- cluster3$predicted.id
pred9 <- WhichCells(cluster3, idents = '9')

int9 <- intersect(sub9, pred9) #1696

#cluster 10 
Idents(cluster3) <- cluster3$seurat_clusters
sub10 <- WhichCells(cluster3, idents = c('14'))

Idents(cluster3) <- cluster3$predicted.id
pred10 <- WhichCells(cluster3, idents = '10')

int10 <- intersect(sub10, pred10) #1366

#cluster 11
Idents(cluster3) <- cluster3$seurat_clusters
sub11 <- WhichCells(cluster3, idents = c('4', '26'))

Idents(cluster3) <- cluster3$predicted.id
pred11 <- WhichCells(cluster3, idents = '11')

int11 <- intersect(sub11, pred11) #1259

#cluster 12
Idents(cluster3) <- cluster3$seurat_clusters
sub12 <- WhichCells(cluster3, idents = c('18'))

Idents(cluster3) <- cluster3$predicted.id
pred12 <- WhichCells(cluster3, idents = '12')

int12 <- intersect(sub12, pred12) #1083

#cluster 13
Idents(cluster3) <- cluster3$seurat_clusters
sub13 <- WhichCells(cluster3, idents = c('20'))

Idents(cluster3) <- cluster3$predicted.id
pred13 <- WhichCells(cluster3, idents = '13')

int13 <- intersect(sub13, pred13) #1010

#cluster 14
Idents(cluster3) <- cluster3$seurat_clusters
sub14 <- WhichCells(cluster3, idents = c('11'))

Idents(cluster3) <- cluster3$predicted.id
pred14 <- WhichCells(cluster3, idents = '14')

int14 <- intersect(sub14, pred14) #1110

#cluster 15
Idents(cluster3) <- cluster3$seurat_clusters
sub15 <- WhichCells(cluster3, idents = c('12', '22'))

Idents(cluster3) <- cluster3$predicted.id
pred15 <- WhichCells(cluster3, idents = '15')

int15 <- intersect(sub15, pred15) #1196

#cluster 16
Idents(cluster3) <- cluster3$seurat_clusters
sub16 <- WhichCells(cluster3, idents = c('8'))

Idents(cluster3) <- cluster3$predicted.id
pred16 <- WhichCells(cluster3, idents = '16')

int16 <- intersect(sub16, pred16) #917

#cluster 17 
Idents(cluster3) <- cluster3$seurat_clusters
sub17 <- WhichCells(cluster3, idents = c('17'))

Idents(cluster3) <- cluster3$predicted.id
pred17 <- WhichCells(cluster3, idents = '17')

int17 <- intersect(sub17, pred17) #1157

#cluster 18 
Idents(cluster3) <- cluster3$seurat_clusters
sub18 <- WhichCells(cluster3, idents = c('12'))

Idents(cluster3) <- cluster3$predicted.id
pred18 <- WhichCells(cluster3, idents = '18')

int18 <- intersect(sub18, pred18) #1021

#cluster 19 
Idents(cluster3) <- cluster3$seurat_clusters
sub19 <- WhichCells(cluster3, idents = c('23'))

Idents(cluster3) <- cluster3$predicted.id
pred19 <- WhichCells(cluster3, idents = '19')

int19 <- intersect(sub19, pred19) #719

#cluster 20 
Idents(cluster3) <- cluster3$seurat_clusters
sub20 <- WhichCells(cluster3, idents = c('8', '27'))

Idents(cluster3) <- cluster3$predicted.id
pred20 <- WhichCells(cluster3, idents = '20')

int20 <- intersect(sub20, pred20) #1059

#cluster 21
Idents(cluster3) <- cluster3$seurat_clusters
sub21 <- WhichCells(cluster3, idents = c('28'))

Idents(cluster3) <- cluster3$predicted.id
pred21 <- WhichCells(cluster3, idents = '21')

int21 <- intersect(sub21, pred21) #517

#cluster 22
Idents(cluster3) <- cluster3$seurat_clusters
sub22 <- WhichCells(cluster3, idents = c('7', '31'))

Idents(cluster3) <- cluster3$predicted.id
pred22 <- WhichCells(cluster3, idents = '22')

int22 <- intersect(sub22, pred22) #398

#cluster 23
Idents(cluster3) <- cluster3$seurat_clusters
sub23 <- WhichCells(cluster3, idents = c('30'))

Idents(cluster3) <- cluster3$predicted.id
pred23 <- WhichCells(cluster3, idents = '23')

int23 <- intersect(sub23, pred23) #434

#cluster 24
Idents(cluster3) <- cluster3$seurat_clusters
sub24 <- WhichCells(cluster3, idents = c('29'))

Idents(cluster3) <- cluster3$predicted.id
pred24 <- WhichCells(cluster3, idents = '24')

int24 <- intersect(sub24, pred24) #449

#cluster 25
Idents(cluster3) <- cluster3$seurat_clusters
sub25 <- WhichCells(cluster3, idents = c('24'))

Idents(cluster3) <- cluster3$predicted.id
pred25 <- WhichCells(cluster3, idents = '25')

int25 <- intersect(sub25, pred25) #306

#cluster 26
Idents(cluster3) <- cluster3$seurat_clusters
sub26 <- WhichCells(cluster3, idents = c('32'))

Idents(cluster3) <- cluster3$predicted.id
pred26 <- WhichCells(cluster3, idents = '26')

int26 <- intersect(sub26, pred26) #347

#cluster 27
Idents(cluster3) <- cluster3$seurat_clusters
sub27 <- WhichCells(cluster3, idents = c('33'))

Idents(cluster3) <- cluster3$predicted.id
pred27 <- WhichCells(cluster3, idents = '27')

int27 <- intersect(sub27, pred27) #309

#cluster 28
Idents(cluster3) <- cluster3$seurat_clusters
sub28 <- WhichCells(cluster3, idents = c('24'))

Idents(cluster3) <- cluster3$predicted.id
pred28 <- WhichCells(cluster3, idents = '28')

int28 <- intersect(sub28, pred28) #279

#cluster 29
Idents(cluster3) <- cluster3$seurat_clusters
sub29 <- WhichCells(cluster3, idents = c('35'))

Idents(cluster3) <- cluster3$predicted.id
pred29 <- WhichCells(cluster3, idents = '29')

int29 <- intersect(sub29, pred29) #205

#cluster 30
Idents(cluster3) <- cluster3$seurat_clusters
sub30 <- WhichCells(cluster3, idents = c('34'))

Idents(cluster3) <- cluster3$predicted.id
pred30 <- WhichCells(cluster3, idents = '30')

int30 <- intersect(sub30, pred30) #213

#cluster 31
Idents(cluster3) <- cluster3$seurat_clusters
sub31 <- WhichCells(cluster3, idents = c('4', '24', '36'))

Idents(cluster3) <- cluster3$predicted.id
pred31 <- WhichCells(cluster3, idents = '31')

int31 <- intersect(sub31, pred31) #71

#cluster 32
Idents(cluster3) <- cluster3$seurat_clusters
sub32 <- WhichCells(cluster3, idents = c('37'))

Idents(cluster3) <- cluster3$predicted.id
pred32 <- WhichCells(cluster3, idents = '32')

int32 <- intersect(sub32, pred32) #76

#cluster 33
Idents(cluster3) <- cluster3$seurat_clusters
sub33 <- WhichCells(cluster3, idents = c('19', '26'))

Idents(cluster3) <- cluster3$predicted.id
pred33 <- WhichCells(cluster3, idents = '33')

int33 <- intersect(sub33, pred33) #46

#cluster 34
Idents(cluster3) <- cluster3$seurat_clusters
sub34 <- WhichCells(cluster3, idents = c('11'))

Idents(cluster3) <- cluster3$predicted.id
pred34 <- WhichCells(cluster3, idents = '34')

int34 <- intersect(sub34, pred34) #58

#cluster 35
Idents(cluster3) <- cluster3$seurat_clusters
sub35 <- WhichCells(cluster3, idents = c('39'))

Idents(cluster3) <- cluster3$predicted.id
pred35 <- WhichCells(cluster3, idents = '35')

int35 <- intersect(sub35, pred35) #34

#cluster 36
Idents(cluster3) <- cluster3$seurat_clusters
sub36 <- WhichCells(cluster3, idents = c('9'))

Idents(cluster3) <- cluster3$predicted.id
pred36 <- WhichCells(cluster3, idents = '36')

int36 <- intersect(sub36, pred36) #49

#cluster 37
Idents(cluster3) <- cluster3$seurat_clusters
sub37 <- WhichCells(cluster3, idents = c('38'))

Idents(cluster3) <- cluster3$predicted.id
pred37 <- WhichCells(cluster3, idents = '37')

int37 <- intersect(sub37, pred37) #49

#cluster 38 
Idents(cluster3) <- cluster3$seurat_clusters
sub38 <- WhichCells(cluster3, idents = c('5', '8'))

Idents(cluster3) <- cluster3$predicted.id
pred38 <- WhichCells(cluster3, idents = '38')

int38 <- intersect(sub38, pred38) #15

#total 
#pre:43892
#39927
total3 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total3,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query3_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster3$query_label <- paste(cluster3$seurat_clusters, "query", sep = "_")
cluster3$ref_label <- paste(cluster3$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster3$ref_label), paste0(cluster3$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-------Query 4---------#

cluster4 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub2_res15_2_high_res_cluster_pred.rds")

cluster4$predicted.id <- factor(cluster4$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster4, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster4, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table4 <- table(cluster4$predicted.id, cluster4$seurat_clusters)
table4 <- as.data.frame.matrix(table4)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table4_prop <- as.data.frame.matrix(apply(table4, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster4) <- cluster4$seurat_clusters
sub0 <- WhichCells(cluster4, idents = c('0'))

Idents(cluster4) <- cluster4$predicted.id
pred0 <- WhichCells(cluster4, idents = '0')

int0 <- intersect(sub0, pred0) #4570

#cluster 1
Idents(cluster4) <- cluster4$seurat_clusters
sub1 <- WhichCells(cluster4, idents = c('1'))

Idents(cluster4) <- cluster4$predicted.id
pred1 <- WhichCells(cluster4, idents = '1')

int1 <- intersect(sub1, pred1) #3206

#cluster 2
Idents(cluster4) <- cluster4$seurat_clusters
sub2 <- WhichCells(cluster4, idents = c('2'))

Idents(cluster4) <- cluster4$predicted.id
pred2 <- WhichCells(cluster4, idents = '2')

int2 <- intersect(sub2, pred2) #2231

#cluster 3
Idents(cluster4) <- cluster4$seurat_clusters
sub3 <- WhichCells(cluster4, idents = c('3', '6'))

Idents(cluster4) <- cluster4$predicted.id
pred3 <- WhichCells(cluster4, idents = '3')

int3 <- intersect(sub3, pred3) #2103

#cluster 4
Idents(cluster4) <- cluster4$seurat_clusters
sub4 <- WhichCells(cluster4, idents = c('3', '6'))

Idents(cluster4) <- cluster4$predicted.id
pred4 <- WhichCells(cluster4, idents = '4')

int4 <- intersect(sub4, pred4) #2083

#cluster 5
Idents(cluster4) <- cluster4$seurat_clusters
sub5 <- WhichCells(cluster4, idents = c('4'))

Idents(cluster4) <- cluster4$predicted.id
pred5 <- WhichCells(cluster4, idents = '5')

int5 <- intersect(sub5, pred5) #2024

#cluster 6
Idents(cluster4) <- cluster4$seurat_clusters
sub6 <- WhichCells(cluster4, idents = c('5', '16'))

Idents(cluster4) <- cluster4$predicted.id
pred6 <- WhichCells(cluster4, idents = '6')

int6 <- intersect(sub6, pred6) #2017

#cluster 7
Idents(cluster4) <- cluster4$seurat_clusters
sub7 <- WhichCells(cluster4, idents = c('5', '7'))

Idents(cluster4) <- cluster4$predicted.id
pred7 <- WhichCells(cluster4, idents = '7')

int7 <- intersect(sub7, pred7) #1830

#cluster 8
Idents(cluster4) <- cluster4$seurat_clusters
sub8 <- WhichCells(cluster4, idents = c('9', '26'))

Idents(cluster4) <- cluster4$predicted.id
pred8 <- WhichCells(cluster4, idents = '8')

int8 <- intersect(sub8, pred8) #1868

#cluster 9 
Idents(cluster4) <- cluster4$seurat_clusters
sub9 <- WhichCells(cluster4, idents = c('18', '21'))

Idents(cluster4) <- cluster4$predicted.id
pred9 <- WhichCells(cluster4, idents = '9')

int9 <- intersect(sub9, pred9) #1871

#cluster 10 
Idents(cluster4) <- cluster4$seurat_clusters
sub10 <- WhichCells(cluster4, idents = c('8'))

Idents(cluster4) <- cluster4$predicted.id
pred10 <- WhichCells(cluster4, idents = '10')

int10 <- intersect(sub10, pred10) #1470

#cluster 11
Idents(cluster4) <- cluster4$seurat_clusters
sub11 <- WhichCells(cluster4, idents = c('14', '17'))

Idents(cluster4) <- cluster4$predicted.id
pred11 <- WhichCells(cluster4, idents = '11')

int11 <- intersect(sub11, pred11) #1111

#cluster 12
Idents(cluster4) <- cluster4$seurat_clusters
sub12 <- WhichCells(cluster4, idents = c('13'))

Idents(cluster4) <- cluster4$predicted.id
pred12 <- WhichCells(cluster4, idents = '12')

int12 <- intersect(sub12, pred12) #998

#cluster 13
Idents(cluster4) <- cluster4$seurat_clusters
sub13 <- WhichCells(cluster4, idents = c('14', '20'))

Idents(cluster4) <- cluster4$predicted.id
pred13 <- WhichCells(cluster4, idents = '13')

int13 <- intersect(sub13, pred13) #1355

#cluster 14
Idents(cluster4) <- cluster4$seurat_clusters
sub14 <- WhichCells(cluster4, idents = c('11'))

Idents(cluster4) <- cluster4$predicted.id
pred14 <- WhichCells(cluster4, idents = '14')

int14 <- intersect(sub14, pred14) #1206

#cluster 15
Idents(cluster4) <- cluster4$seurat_clusters
sub15 <- WhichCells(cluster4, idents = c('12'))

Idents(cluster4) <- cluster4$predicted.id
pred15 <- WhichCells(cluster4, idents = '15')

int15 <- intersect(sub15, pred15) #1236

#cluster 16
Idents(cluster4) <- cluster4$seurat_clusters
sub16 <- WhichCells(cluster4, idents = c('19'))

Idents(cluster4) <- cluster4$predicted.id
pred16 <- WhichCells(cluster4, idents = '16')

int16 <- intersect(sub16, pred16) #954

#cluster 17
Idents(cluster4) <- cluster4$seurat_clusters
sub17 <- WhichCells(cluster4, idents = c('10'))

Idents(cluster4) <- cluster4$predicted.id
pred17 <- WhichCells(cluster4, idents = '17')

int17 <- intersect(sub17, pred17) #1038

#cluster 18 
Idents(cluster4) <- cluster4$seurat_clusters
sub18 <- WhichCells(cluster4, idents = c('15'))

Idents(cluster4) <- cluster4$predicted.id
pred18 <- WhichCells(cluster4, idents = '18')

int18 <- intersect(sub18, pred18) #971

#cluster 19
Idents(cluster4) <- cluster4$seurat_clusters
sub19 <- WhichCells(cluster4, idents = c('23'))

Idents(cluster4) <- cluster4$predicted.id
pred19 <- WhichCells(cluster4, idents = '19')

int19 <- intersect(sub19, pred19) #829

#cluster 20
Idents(cluster4) <- cluster4$seurat_clusters
sub20 <- WhichCells(cluster4, idents = c('22'))

Idents(cluster4) <- cluster4$predicted.id
pred20 <- WhichCells(cluster4, idents = '20')

int20 <- intersect(sub20, pred20) #869

#cluster 21
Idents(cluster4) <- cluster4$seurat_clusters
sub21 <- WhichCells(cluster4, idents = c('25'))

Idents(cluster4) <- cluster4$predicted.id
pred21 <- WhichCells(cluster4, idents = '21')

int21 <- intersect(sub21, pred21) #627

#cluster 22
Idents(cluster4) <- cluster4$seurat_clusters
sub22 <- WhichCells(cluster4, idents = c('9', '28'))

Idents(cluster4) <- cluster4$predicted.id
pred22 <- WhichCells(cluster4, idents = '22')

int22 <- intersect(sub22, pred22) #585

#cluster 23
Idents(cluster4) <- cluster4$seurat_clusters
sub23 <- WhichCells(cluster4, idents = c('24'))

Idents(cluster4) <- cluster4$predicted.id
pred23 <- WhichCells(cluster4, idents = '23')

int23 <- intersect(sub23, pred23) #457

#cluster 24
Idents(cluster4) <- cluster4$seurat_clusters
sub24 <- WhichCells(cluster4, idents = c('27'))

Idents(cluster4) <- cluster4$predicted.id
pred24 <- WhichCells(cluster4, idents = '24')

int24 <- intersect(sub24, pred24) #459

#cluster 25
Idents(cluster4) <- cluster4$seurat_clusters
sub25 <- WhichCells(cluster4, idents = c('16'))

Idents(cluster4) <- cluster4$predicted.id
pred25 <- WhichCells(cluster4, idents = '25')

int25 <- intersect(sub25, pred25) #339

#cluster 26
Idents(cluster4) <- cluster4$seurat_clusters
sub26 <- WhichCells(cluster4, idents = c('24'))

Idents(cluster4) <- cluster4$predicted.id
pred26 <- WhichCells(cluster4, idents = '26')

int26 <- intersect(sub26, pred26) #338

#cluster 27
Idents(cluster4) <- cluster4$seurat_clusters
sub27 <- WhichCells(cluster4, idents = c('29'))

Idents(cluster4) <- cluster4$predicted.id
pred27 <- WhichCells(cluster4, idents = '27')

int27 <- intersect(sub27, pred27) #320

#cluster 28 
Idents(cluster4) <- cluster4$seurat_clusters
sub28 <- WhichCells(cluster4, idents = c('16'))

Idents(cluster4) <- cluster4$predicted.id
pred28 <- WhichCells(cluster4, idents = '28')

int28 <- intersect(sub28, pred28) #233

#cluster 29 
Idents(cluster4) <- cluster4$seurat_clusters
sub29 <- WhichCells(cluster4, idents = c('30'))

Idents(cluster4) <- cluster4$predicted.id
pred29 <- WhichCells(cluster4, idents = '29')

int29 <- intersect(sub29, pred29) #180

#cluster 30 
Idents(cluster4) <- cluster4$seurat_clusters
sub30 <- WhichCells(cluster4, idents = c('31'))

Idents(cluster4) <- cluster4$predicted.id
pred30 <- WhichCells(cluster4, idents = '30')

int30 <- intersect(sub30, pred30) #161

#cluster 31
Idents(cluster4) <- cluster4$seurat_clusters
sub31 <- WhichCells(cluster4, idents = c('17', '27', '33'))

Idents(cluster4) <- cluster4$predicted.id
pred31 <- WhichCells(cluster4, idents = '31')

int31 <- intersect(sub31, pred31) #108

#cluster 32
Idents(cluster4) <- cluster4$seurat_clusters
sub32 <- WhichCells(cluster4, idents = c('32'))

Idents(cluster4) <- cluster4$predicted.id
pred32 <- WhichCells(cluster4, idents = '32')

int32 <- intersect(sub32, pred32) #81

#cluster 33
Idents(cluster4) <- cluster4$seurat_clusters
sub33 <- WhichCells(cluster4, idents = c('14'))

Idents(cluster4) <- cluster4$predicted.id
pred33 <- WhichCells(cluster4, idents = '33')

int33 <- intersect(sub33, pred33) #72

#cluster 34
Idents(cluster4) <- cluster4$seurat_clusters
sub34 <- WhichCells(cluster4, idents = c('21'))

Idents(cluster4) <- cluster4$predicted.id
pred34 <- WhichCells(cluster4, idents = '34')

int34 <- intersect(sub34, pred34) #43

#cluster 35
Idents(cluster4) <- cluster4$seurat_clusters
sub35 <- WhichCells(cluster4, idents = c('35'))

Idents(cluster4) <- cluster4$predicted.id
pred35 <- WhichCells(cluster4, idents = '35')

int35 <- intersect(sub35, pred35) #39

#cluster 36
Idents(cluster4) <- cluster4$seurat_clusters
sub36 <- WhichCells(cluster4, idents = c('34'))

Idents(cluster4) <- cluster4$predicted.id
pred36 <- WhichCells(cluster4, idents = '36')

int36 <- intersect(sub36, pred36) #45

#cluster 37
Idents(cluster4) <- cluster4$seurat_clusters
sub37 <- WhichCells(cluster4, idents = c('30'))

Idents(cluster4) <- cluster4$predicted.id
pred37 <- WhichCells(cluster4, idents = '37')

int37 <- intersect(sub37, pred37) #45

#cluster 38
Idents(cluster4) <- cluster4$seurat_clusters
sub38 <- WhichCells(cluster4, idents = c('19'))

Idents(cluster4) <- cluster4$predicted.id
pred38 <- WhichCells(cluster4, idents = '38')

int38 <- intersect(sub38, pred38) #21

#total 
#pre:43893
#39993
total4 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total4,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query4_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster4$query_label <- paste(cluster4$seurat_clusters, "query", sep = "_")
cluster4$ref_label <- paste(cluster4$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster4$ref_label), paste0(cluster4$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-------Query 5---------#

cluster5 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub_res15_3_high_res_cluster_pred.rds")

cluster5$predicted.id <- factor(cluster5$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster5, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster5, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table5 <- table(cluster5$predicted.id, cluster5$seurat_clusters)
table5 <- as.data.frame.matrix(table5)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table5_prop <- as.data.frame.matrix(apply(table5, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster5) <- cluster5$seurat_clusters
sub0 <- WhichCells(cluster5, idents = c('0'))

Idents(cluster5) <- cluster5$predicted.id
pred0 <- WhichCells(cluster5, idents = '0')

int0 <- intersect(sub0, pred0) #4092

#cluster 1
Idents(cluster5) <- cluster5$seurat_clusters
sub1 <- WhichCells(cluster5, idents = c('1'))

Idents(cluster5) <- cluster5$predicted.id
pred1 <- WhichCells(cluster5, idents = '1')

int1 <- intersect(sub1, pred1) #2794

#cluster 2
Idents(cluster5) <- cluster5$seurat_clusters
sub2 <- WhichCells(cluster5, idents = c('3'))

Idents(cluster5) <- cluster5$predicted.id
pred2 <- WhichCells(cluster5, idents = '2')

int2 <- intersect(sub2, pred2) #2060

#cluster 3
Idents(cluster5) <- cluster5$seurat_clusters
sub3 <- WhichCells(cluster5, idents = c('5', '11'))

Idents(cluster5) <- cluster5$predicted.id
pred3 <- WhichCells(cluster5, idents = '3')

int3 <- intersect(sub3, pred3) #2361

#cluster 4
Idents(cluster5) <- cluster5$seurat_clusters
sub4 <- WhichCells(cluster5, idents = c('8', '11'))

Idents(cluster5) <- cluster5$predicted.id
pred4 <- WhichCells(cluster5, idents = '4')

int4 <- intersect(sub4, pred4) #2206

#cluster 5
Idents(cluster5) <- cluster5$seurat_clusters
sub5 <- WhichCells(cluster5, idents = c('2'))

Idents(cluster5) <- cluster5$predicted.id
pred5 <- WhichCells(cluster5, idents = '5')

int5 <- intersect(sub5, pred5) #2336

#cluster 6
Idents(cluster5) <- cluster5$seurat_clusters
sub6 <- WhichCells(cluster5, idents = c('4', '15'))

Idents(cluster5) <- cluster5$predicted.id
pred6 <- WhichCells(cluster5, idents = '6')

int6 <- intersect(sub6, pred6) #2278

#cluster 7 
Idents(cluster5) <- cluster5$seurat_clusters
sub7 <- WhichCells(cluster5, idents = c('5', '9'))

Idents(cluster5) <- cluster5$predicted.id
pred7 <- WhichCells(cluster5, idents = '7')

int7 <- intersect(sub7, pred7) #1848

#cluster 8 
Idents(cluster5) <- cluster5$seurat_clusters
sub8 <- WhichCells(cluster5, idents = c('16', '24'))

Idents(cluster5) <- cluster5$predicted.id
pred8 <- WhichCells(cluster5, idents = '8')

int8 <- intersect(sub8, pred8) #1800

#cluster 9 
Idents(cluster5) <- cluster5$seurat_clusters
sub9 <- WhichCells(cluster5, idents = c('7'))

Idents(cluster5) <- cluster5$predicted.id
pred9 <- WhichCells(cluster5, idents = '9')

int9 <- intersect(sub9, pred9) #1876

#cluster 10 
Idents(cluster5) <- cluster5$seurat_clusters
sub10 <- WhichCells(cluster5, idents = c('12'))

Idents(cluster5) <- cluster5$predicted.id
pred10 <- WhichCells(cluster5, idents = '10')

int10 <- intersect(sub10, pred10) #1435

#cluster 11
Idents(cluster5) <- cluster5$seurat_clusters
sub11 <- WhichCells(cluster5, idents = c('5', '15', '19'))

Idents(cluster5) <- cluster5$predicted.id
pred11 <- WhichCells(cluster5, idents = '11')

int11 <- intersect(sub11, pred11) #1097

#cluster 12
Idents(cluster5) <- cluster5$seurat_clusters
sub12 <- WhichCells(cluster5, idents = c('4', '6'))

Idents(cluster5) <- cluster5$predicted.id
pred12 <- WhichCells(cluster5, idents = '12')

int12 <- intersect(sub12, pred12) #1211

#cluster 13
Idents(cluster5) <- cluster5$seurat_clusters
sub13 <- WhichCells(cluster5, idents = c('18', '19'))

Idents(cluster5) <- cluster5$predicted.id
pred13 <- WhichCells(cluster5, idents = '13')

int13 <- intersect(sub13, pred13) #1243

#cluster 14
Idents(cluster5) <- cluster5$seurat_clusters
sub14 <- WhichCells(cluster5, idents = c('14'))

Idents(cluster5) <- cluster5$predicted.id
pred14 <- WhichCells(cluster5, idents = '14')

int14 <- intersect(sub14, pred14) #1215

#cluster 15
Idents(cluster5) <- cluster5$seurat_clusters
sub15 <- WhichCells(cluster5, idents = c('13', '20'))

Idents(cluster5) <- cluster5$predicted.id
pred15 <- WhichCells(cluster5, idents = '15')

int15 <- intersect(sub15, pred15) #1158

#cluster 16
Idents(cluster5) <- cluster5$seurat_clusters
sub16 <- WhichCells(cluster5, idents = c('10'))

Idents(cluster5) <- cluster5$predicted.id
pred16 <- WhichCells(cluster5, idents = '16')

int16 <- intersect(sub16, pred16) #1056

#cluster 17 
Idents(cluster5) <- cluster5$seurat_clusters
sub17 <- WhichCells(cluster5, idents = c('17'))

Idents(cluster5) <- cluster5$predicted.id
pred17 <- WhichCells(cluster5, idents = '17')

int17 <- intersect(sub17, pred17) #1027

#cluster 18
Idents(cluster5) <- cluster5$seurat_clusters
sub18 <- WhichCells(cluster5, idents = c('13'))

Idents(cluster5) <- cluster5$predicted.id
pred18 <- WhichCells(cluster5, idents = '18')

int18 <- intersect(sub18, pred18) #1022

#cluster 19 
Idents(cluster5) <- cluster5$seurat_clusters
sub19 <- WhichCells(cluster5, idents = c('21'))

Idents(cluster5) <- cluster5$predicted.id
pred19 <- WhichCells(cluster5, idents = '19')

int19 <- intersect(sub19, pred19) #779

#cluster 20 
Idents(cluster5) <- cluster5$seurat_clusters
sub20 <- WhichCells(cluster5, idents = c('10', '27'))

Idents(cluster5) <- cluster5$predicted.id
pred20 <- WhichCells(cluster5, idents = '20')

int20 <- intersect(sub20, pred20) #927

#cluster 21
Idents(cluster5) <- cluster5$seurat_clusters
sub21 <- WhichCells(cluster5, idents = c('25'))

Idents(cluster5) <- cluster5$predicted.id
pred21 <- WhichCells(cluster5, idents = '21')

int21 <- intersect(sub21, pred21) #578

#cluster 22
Idents(cluster5) <- cluster5$seurat_clusters
sub22 <- WhichCells(cluster5, idents = c('23'))

Idents(cluster5) <- cluster5$predicted.id
pred22 <- WhichCells(cluster5, idents = '22')

int22 <- intersect(sub22, pred22) #584

#cluster 23
Idents(cluster5) <- cluster5$seurat_clusters
sub23 <- WhichCells(cluster5, idents = c('29'))

Idents(cluster5) <- cluster5$predicted.id
pred23 <- WhichCells(cluster5, idents = '23')

int23 <- intersect(sub23, pred23) #400

#cluster 24
Idents(cluster5) <- cluster5$seurat_clusters
sub24 <- WhichCells(cluster5, idents = c('26'))

Idents(cluster5) <- cluster5$predicted.id
pred24 <- WhichCells(cluster5, idents = '24')

int24 <- intersect(sub24, pred24) #487

#cluster 25
Idents(cluster5) <- cluster5$seurat_clusters
sub25 <- WhichCells(cluster5, idents = c('22'))

Idents(cluster5) <- cluster5$predicted.id
pred25 <- WhichCells(cluster5, idents = '25')

int25 <- intersect(sub25, pred25) #279

#cluster 26 
Idents(cluster5) <- cluster5$seurat_clusters
sub26 <- WhichCells(cluster5, idents = c('28'))

Idents(cluster5) <- cluster5$predicted.id
pred26 <- WhichCells(cluster5, idents = '26')

int26 <- intersect(sub26, pred26) #358

#cluster 27 
Idents(cluster5) <- cluster5$seurat_clusters
sub27 <- WhichCells(cluster5, idents = c('30'))

Idents(cluster5) <- cluster5$predicted.id
pred27 <- WhichCells(cluster5, idents = '27')

int27 <- intersect(sub27, pred27) #354

#cluster 28 
Idents(cluster5) <- cluster5$seurat_clusters
sub28 <- WhichCells(cluster5, idents = c('22'))

Idents(cluster5) <- cluster5$predicted.id
pred28 <- WhichCells(cluster5, idents = '28')

int28 <- intersect(sub28, pred28) #250

#cluster 29 
Idents(cluster5) <- cluster5$seurat_clusters
sub29 <- WhichCells(cluster5, idents = c('32'))

Idents(cluster5) <- cluster5$predicted.id
pred29 <- WhichCells(cluster5, idents = '29')

int29 <- intersect(sub29, pred29) #198

#cluster 30 
Idents(cluster5) <- cluster5$seurat_clusters
sub30 <- WhichCells(cluster5, idents = c('31'))

Idents(cluster5) <- cluster5$predicted.id
pred30 <- WhichCells(cluster5, idents = '30')

int30 <- intersect(sub30, pred30) #205

#cluster 31
Idents(cluster5) <- cluster5$seurat_clusters
sub31 <- WhichCells(cluster5, idents = c('15'))

Idents(cluster5) <- cluster5$predicted.id
pred31 <- WhichCells(cluster5, idents = '31')

int31 <- intersect(sub31, pred31) #109

#cluster 32
Idents(cluster5) <- cluster5$seurat_clusters
sub32 <- WhichCells(cluster5, idents = c('33'))

Idents(cluster5) <- cluster5$predicted.id
pred32 <- WhichCells(cluster5, idents = '32')

int32 <- intersect(sub32, pred32) #93

#cluster 33
Idents(cluster5) <- cluster5$seurat_clusters
sub33 <- WhichCells(cluster5, idents = c('5'))

Idents(cluster5) <- cluster5$predicted.id
pred33 <- WhichCells(cluster5, idents = '33')

int33 <- intersect(sub33, pred33) #60

#cluster 34
Idents(cluster5) <- cluster5$seurat_clusters
sub34 <- WhichCells(cluster5, idents = c('14'))

Idents(cluster5) <- cluster5$predicted.id
pred34 <- WhichCells(cluster5, idents = '34')

int34 <- intersect(sub34, pred34) #69

#cluster 35
Idents(cluster5) <- cluster5$seurat_clusters
sub35 <- WhichCells(cluster5, idents = c('22'))

Idents(cluster5) <- cluster5$predicted.id
pred35 <- WhichCells(cluster5, idents = '35')

int35 <- intersect(sub35, pred35) #43

#cluster 36
Idents(cluster5) <- cluster5$seurat_clusters
sub36 <- WhichCells(cluster5, idents = c('35'))

Idents(cluster5) <- cluster5$predicted.id
pred36 <- WhichCells(cluster5, idents = '36')

int36 <- intersect(sub36, pred36) #47

#cluster 37
Idents(cluster5) <- cluster5$seurat_clusters
sub37 <- WhichCells(cluster5, idents = c('34'))

Idents(cluster5) <- cluster5$predicted.id
pred37 <- WhichCells(cluster5, idents = '37')

int37 <- intersect(sub37, pred37) #51

#cluster 38 
Idents(cluster5) <- cluster5$seurat_clusters
sub38 <- WhichCells(cluster5, idents = c('36'))

Idents(cluster5) <- cluster5$predicted.id
pred38 <- WhichCells(cluster5, idents = '38')

int38 <- intersect(sub38, pred38) #25

#total 
#pre:43892
#40011
total5 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total5,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query5_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster5$query_label <- paste(cluster5$seurat_clusters, "query", sep = "_")
cluster5$ref_label <- paste(cluster5$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster5$ref_label), paste0(cluster5$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)


#-------Query 6---------#

cluster6 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub2_res15_3_high_res_cluster_pred.rds")

cluster6$predicted.id <- factor(cluster6$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster6, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster6, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table6 <- table(cluster6$predicted.id, cluster6$seurat_clusters)
table6 <- as.data.frame.matrix(table6)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table6_prop <- as.data.frame.matrix(apply(table6, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster6) <- cluster6$seurat_clusters
sub0 <- WhichCells(cluster6, idents = c('0'))

Idents(cluster6) <- cluster6$predicted.id
pred0 <- WhichCells(cluster6, idents = '0')

int0 <- intersect(sub0, pred0) #4286

#cluster 1
Idents(cluster6) <- cluster6$seurat_clusters
sub1 <- WhichCells(cluster6, idents = c('1'))

Idents(cluster6) <- cluster6$predicted.id
pred1 <- WhichCells(cluster6, idents = '1')

int1 <- intersect(sub1, pred1) #3136

#cluster 2
Idents(cluster6) <- cluster6$seurat_clusters
sub2 <- WhichCells(cluster6, idents = c('3'))

Idents(cluster6) <- cluster6$predicted.id
pred2 <- WhichCells(cluster6, idents = '2')

int2 <- intersect(sub2, pred2) #2323

#cluster 3
Idents(cluster6) <- cluster6$seurat_clusters
sub3 <- WhichCells(cluster6, idents = c('4', '6'))

Idents(cluster6) <- cluster6$predicted.id
pred3 <- WhichCells(cluster6, idents = '3')

int3 <- intersect(sub3, pred3) #2330

#cluster 4
Idents(cluster6) <- cluster6$seurat_clusters
sub4 <- WhichCells(cluster6, idents = c('4', '13'))

Idents(cluster6) <- cluster6$predicted.id
pred4 <- WhichCells(cluster6, idents = '4')

int4 <- intersect(sub4, pred4) #2078

#cluster 5
Idents(cluster6) <- cluster6$seurat_clusters
sub5 <- WhichCells(cluster6, idents = c('2'))

Idents(cluster6) <- cluster6$predicted.id
pred5 <- WhichCells(cluster6, idents = '5')

int5 <- intersect(sub5, pred5) #2365

#cluster 6 
Idents(cluster6) <- cluster6$seurat_clusters
sub6 <- WhichCells(cluster6, idents = c('8', '11', '17'))

Idents(cluster6) <- cluster6$predicted.id
pred6 <- WhichCells(cluster6, idents = '6')

int6 <- intersect(sub6, pred6) #2120

#cluster 7 
Idents(cluster6) <- cluster6$seurat_clusters
sub7 <- WhichCells(cluster6, idents = c('5'))

Idents(cluster6) <- cluster6$predicted.id
pred7 <- WhichCells(cluster6, idents = '7')

int7 <- intersect(sub7, pred7) #1761

#cluster 8 
Idents(cluster6) <- cluster6$seurat_clusters
sub8 <- WhichCells(cluster6, idents = c('10', '28'))

Idents(cluster6) <- cluster6$predicted.id
pred8 <- WhichCells(cluster6, idents = '8')

int8 <- intersect(sub8, pred8) #1846

#cluster 9 
Idents(cluster6) <- cluster6$seurat_clusters
sub9 <- WhichCells(cluster6, idents = c('14', '21'))

Idents(cluster6) <- cluster6$predicted.id
pred9 <- WhichCells(cluster6, idents = '9')

int9 <- intersect(sub9, pred9) #1856

#cluster 10 
Idents(cluster6) <- cluster6$seurat_clusters
sub10 <- WhichCells(cluster6, idents = c('9'))

Idents(cluster6) <- cluster6$predicted.id
pred10 <- WhichCells(cluster6, idents = '10')

int10 <- intersect(sub10, pred10) #1416

#cluster 11
Idents(cluster6) <- cluster6$seurat_clusters
sub11 <- WhichCells(cluster6, idents = c('11', '27'))

Idents(cluster6) <- cluster6$predicted.id
pred11 <- WhichCells(cluster6, idents = '11')

int11 <- intersect(sub11, pred11) #1200

#cluster 12
Idents(cluster6) <- cluster6$seurat_clusters
sub12 <- WhichCells(cluster6, idents = c('16'))

Idents(cluster6) <- cluster6$predicted.id
pred12 <- WhichCells(cluster6, idents = '12')

int12 <- intersect(sub12, pred12) #999

#cluster 13 
Idents(cluster6) <- cluster6$seurat_clusters
sub13 <- WhichCells(cluster6, idents = c('7'))

Idents(cluster6) <- cluster6$predicted.id
pred13 <- WhichCells(cluster6, idents = '13')

int13 <- intersect(sub13, pred13) #1253

#cluster 14
Idents(cluster6) <- cluster6$seurat_clusters
sub14 <- WhichCells(cluster6, idents = c('12'))

Idents(cluster6) <- cluster6$predicted.id
pred14 <- WhichCells(cluster6, idents = '14')

int14 <- intersect(sub14, pred14) #1131

#cluster 15
Idents(cluster6) <- cluster6$seurat_clusters
sub15 <- WhichCells(cluster6, idents = c('22', '23'))

Idents(cluster6) <- cluster6$predicted.id
pred15 <- WhichCells(cluster6, idents = '15')

int15 <- intersect(sub15, pred15) #1242

#cluster 16
Idents(cluster6) <- cluster6$seurat_clusters
sub16 <- WhichCells(cluster6, idents = c('19'))

Idents(cluster6) <- cluster6$predicted.id
pred16 <- WhichCells(cluster6, idents = '16')

int16 <- intersect(sub16, pred16) #895

#cluster 17
Idents(cluster6) <- cluster6$seurat_clusters
sub17 <- WhichCells(cluster6, idents = c('15'))

Idents(cluster6) <- cluster6$predicted.id
pred17 <- WhichCells(cluster6, idents = '17')

int17 <- intersect(sub17, pred17) #1021

#cluster 18 
Idents(cluster6) <- cluster6$seurat_clusters
sub18 <- WhichCells(cluster6, idents = c('23', '24'))

Idents(cluster6) <- cluster6$predicted.id
pred18 <- WhichCells(cluster6, idents = '18')

int18 <- intersect(sub18, pred18) #991

#cluster 19 
Idents(cluster6) <- cluster6$seurat_clusters
sub19 <- WhichCells(cluster6, idents = c('25'))

Idents(cluster6) <- cluster6$predicted.id
pred19 <- WhichCells(cluster6, idents = '19')

int19 <- intersect(sub19, pred19) #727

#cluster 20 
Idents(cluster6) <- cluster6$seurat_clusters
sub20 <- WhichCells(cluster6, idents = c('18'))

Idents(cluster6) <- cluster6$predicted.id
pred20 <- WhichCells(cluster6, idents = '20')

int20 <- intersect(sub20, pred20) #930

#cluster 21
Idents(cluster6) <- cluster6$seurat_clusters
sub21 <- WhichCells(cluster6, idents = c('26'))

Idents(cluster6) <- cluster6$predicted.id
pred21 <- WhichCells(cluster6, idents = '21')

int21 <- intersect(sub21, pred21) #660

#cluster 22
Idents(cluster6) <- cluster6$seurat_clusters
sub22 <- WhichCells(cluster6, idents = c('29'))

Idents(cluster6) <- cluster6$predicted.id
pred22 <- WhichCells(cluster6, idents = '22')

int22 <- intersect(sub22, pred22) #525

#cluster 23
Idents(cluster6) <- cluster6$seurat_clusters
sub23 <- WhichCells(cluster6, idents = c('31'))

Idents(cluster6) <- cluster6$predicted.id
pred23 <- WhichCells(cluster6, idents = '23')

int23 <- intersect(sub23, pred23) #435

#cluster 24
Idents(cluster6) <- cluster6$seurat_clusters
sub24 <- WhichCells(cluster6, idents = c('30'))

Idents(cluster6) <- cluster6$predicted.id
pred24 <- WhichCells(cluster6, idents = '24')

int24 <- intersect(sub24, pred24) #448

#cluster 25
Idents(cluster6) <- cluster6$seurat_clusters
sub25 <- WhichCells(cluster6, idents = c('17'))

Idents(cluster6) <- cluster6$predicted.id
pred25 <- WhichCells(cluster6, idents = '25')

int25 <- intersect(sub25, pred25) #362

#cluster 26
Idents(cluster6) <- cluster6$seurat_clusters
sub26 <- WhichCells(cluster6, idents = c('32'))

Idents(cluster6) <- cluster6$predicted.id
pred26 <- WhichCells(cluster6, idents = '26')

int26 <- intersect(sub26, pred26) #324

#cluster 27 
Idents(cluster6) <- cluster6$seurat_clusters
sub27 <- WhichCells(cluster6, idents = c('33'))

Idents(cluster6) <- cluster6$predicted.id
pred27 <- WhichCells(cluster6, idents = '27')

int27 <- intersect(sub27, pred27) #287

#cluster 28
Idents(cluster6) <- cluster6$seurat_clusters
sub28 <- WhichCells(cluster6, idents = c('17'))

Idents(cluster6) <- cluster6$predicted.id
pred28 <- WhichCells(cluster6, idents = '28')

int28 <- intersect(sub28, pred28) #281

#cluster 29
Idents(cluster6) <- cluster6$seurat_clusters
sub29 <- WhichCells(cluster6, idents = c('34'))

Idents(cluster6) <- cluster6$predicted.id
pred29 <- WhichCells(cluster6, idents = '29')

int29 <- intersect(sub29, pred29) #204

#cluster 30 
Idents(cluster6) <- cluster6$seurat_clusters
sub30 <- WhichCells(cluster6, idents = c('35'))

Idents(cluster6) <- cluster6$predicted.id
pred30 <- WhichCells(cluster6, idents = '30')

int30 <- intersect(sub30, pred30) #197

#cluster 31
Idents(cluster6) <- cluster6$seurat_clusters
sub31 <- WhichCells(cluster6, idents = c('11', '17', '20'))

Idents(cluster6) <- cluster6$predicted.id
pred31 <- WhichCells(cluster6, idents = '31')

int31 <- intersect(sub31, pred31) #111

#cluster 32
Idents(cluster6) <- cluster6$seurat_clusters
sub32 <- WhichCells(cluster6, idents = c('36'))

Idents(cluster6) <- cluster6$predicted.id
pred32 <- WhichCells(cluster6, idents = '32')

int32 <- intersect(sub32, pred32) #79

#cluster 33
Idents(cluster6) <- cluster6$seurat_clusters
sub33 <- WhichCells(cluster6, idents = c('7'))

Idents(cluster6) <- cluster6$predicted.id
pred33 <- WhichCells(cluster6, idents = '33')

int33 <- intersect(sub33, pred33) #66

#cluster 34
Idents(cluster6) <- cluster6$seurat_clusters
sub34 <- WhichCells(cluster6, idents = c('12', '14'))

Idents(cluster6) <- cluster6$predicted.id
pred34 <- WhichCells(cluster6, idents = '34')

int34 <- intersect(sub34, pred34) #38

#cluster 35
Idents(cluster6) <- cluster6$seurat_clusters
sub35 <- WhichCells(cluster6, idents = c('37'))

Idents(cluster6) <- cluster6$predicted.id
pred35 <- WhichCells(cluster6, idents = '35')

int35 <- intersect(sub35, pred35) #41

#cluster 36
Idents(cluster6) <- cluster6$seurat_clusters
sub36 <- WhichCells(cluster6, idents = c('10'))

Idents(cluster6) <- cluster6$predicted.id
pred36 <- WhichCells(cluster6, idents = '36')

int36 <- intersect(sub36, pred36) #50

#cluster 37
Idents(cluster6) <- cluster6$seurat_clusters
sub37 <- WhichCells(cluster6, idents = c('34'))

Idents(cluster6) <- cluster6$predicted.id
pred37 <- WhichCells(cluster6, idents = '37')

int37 <- intersect(sub37, pred37) #49

#cluster 38 
Idents(cluster6) <- cluster6$seurat_clusters
sub38 <- WhichCells(cluster6, idents = c('1', '18', '19'))

Idents(cluster6) <- cluster6$predicted.id
pred38 <- WhichCells(cluster6, idents = '38')

int38 <- intersect(sub38, pred38) #13

#total 
#pre:43893
#40076
total6 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total6,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query6_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster6$query_label <- paste(cluster6$seurat_clusters, "query", sep = "_")
cluster6$ref_label <- paste(cluster6$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster6$ref_label), paste0(cluster6$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-------Query 7---------#

cluster7 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub_res15_4_high_res_cluster_pred.rds")

cluster7$predicted.id <- factor(cluster7$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster7, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster7, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table7 <- table(cluster7$predicted.id, cluster7$seurat_clusters)
table7 <- as.data.frame.matrix(table7)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table7_prop <- as.data.frame.matrix(apply(table7, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster7) <- cluster7$seurat_clusters
sub0 <- WhichCells(cluster7, idents = c('0', '8'))

Idents(cluster7) <- cluster7$predicted.id
pred0 <- WhichCells(cluster7, idents = '0')

int0 <- intersect(sub0, pred0) #4738

#cluster 1
Idents(cluster7) <- cluster7$seurat_clusters
sub1 <- WhichCells(cluster7, idents = c('1'))

Idents(cluster7) <- cluster7$predicted.id
pred1 <- WhichCells(cluster7, idents = '1')

int1 <- intersect(sub1, pred1) #3171

#cluster 2
Idents(cluster7) <- cluster7$seurat_clusters
sub2 <- WhichCells(cluster7, idents = c('2'))

Idents(cluster7) <- cluster7$predicted.id
pred2 <- WhichCells(cluster7, idents = '2')

int2 <- intersect(sub2, pred2) #2378

#cluster 3 
Idents(cluster7) <- cluster7$seurat_clusters
sub3 <- WhichCells(cluster7, idents = c('2', '5'))

Idents(cluster7) <- cluster7$predicted.id
pred3 <- WhichCells(cluster7, idents = '3')

int3 <- intersect(sub3, pred3) #1740

#cluster 4
Idents(cluster7) <- cluster7$seurat_clusters
sub4 <- WhichCells(cluster7, idents = c('3'))

Idents(cluster7) <- cluster7$predicted.id
pred4 <- WhichCells(cluster7, idents = '4')

int4 <- intersect(sub4, pred4) #2178

#cluster 5
Idents(cluster7) <- cluster7$seurat_clusters
sub5 <- WhichCells(cluster7, idents = c('4'))

Idents(cluster7) <- cluster7$predicted.id
pred5 <- WhichCells(cluster7, idents = '5')

int5 <- intersect(sub5, pred5) #2412

#cluster 6
Idents(cluster7) <- cluster7$seurat_clusters
sub6 <- WhichCells(cluster7, idents = c('9', '18'))

Idents(cluster7) <- cluster7$predicted.id
pred6 <- WhichCells(cluster7, idents = '6')

int6 <- intersect(sub6, pred6) #1950

#cluster 7
Idents(cluster7) <- cluster7$seurat_clusters
sub7 <- WhichCells(cluster7, idents = c('6'))

Idents(cluster7) <- cluster7$predicted.id
pred7 <- WhichCells(cluster7, idents = '7')

int7 <- intersect(sub7, pred7) #1753

#cluster 8 
Idents(cluster7) <- cluster7$seurat_clusters
sub8 <- WhichCells(cluster7, idents = c('10', '25'))

Idents(cluster7) <- cluster7$predicted.id
pred8 <- WhichCells(cluster7, idents = '8')

int8 <- intersect(sub8, pred8) #1946

#cluster 9 
Idents(cluster7) <- cluster7$seurat_clusters
sub9 <- WhichCells(cluster7, idents = c('7'))

Idents(cluster7) <- cluster7$predicted.id
pred9 <- WhichCells(cluster7, idents = '9')

int9 <- intersect(sub9, pred9) #1816

#cluster 10 
Idents(cluster7) <- cluster7$seurat_clusters
sub10 <- WhichCells(cluster7, idents = c('16', '29'))

Idents(cluster7) <- cluster7$predicted.id
pred10 <- WhichCells(cluster7, idents = '10')

int10 <- intersect(sub10, pred10) #1390

#cluster 11
Idents(cluster7) <- cluster7$seurat_clusters
sub11 <- WhichCells(cluster7, idents = c('2', '20'))

Idents(cluster7) <- cluster7$predicted.id
pred11 <- WhichCells(cluster7, idents = '11')

int11 <- intersect(sub11, pred11) #1044

#cluster 12
Idents(cluster7) <- cluster7$seurat_clusters
sub12 <- WhichCells(cluster7, idents = c('8', '22'))

Idents(cluster7) <- cluster7$predicted.id
pred12 <- WhichCells(cluster7, idents = '12')

int12 <- intersect(sub12, pred12) #1078

#cluster 13
Idents(cluster7) <- cluster7$seurat_clusters
sub13 <- WhichCells(cluster7, idents = c('11'))

Idents(cluster7) <- cluster7$predicted.id
pred13 <- WhichCells(cluster7, idents = '13')

int13 <- intersect(sub13, pred13) #1180

#cluster 14
Idents(cluster7) <- cluster7$seurat_clusters
sub14 <- WhichCells(cluster7, idents = c('12'))

Idents(cluster7) <- cluster7$predicted.id
pred14 <- WhichCells(cluster7, idents = '14')

int14 <- intersect(sub14, pred14) #1172

#cluster 15
Idents(cluster7) <- cluster7$seurat_clusters
sub15 <- WhichCells(cluster7, idents = c('13'))

Idents(cluster7) <- cluster7$predicted.id
pred15 <- WhichCells(cluster7, idents = '15')

int15 <- intersect(sub15, pred15) #1224

#cluster 16
Idents(cluster7) <- cluster7$seurat_clusters
sub16 <- WhichCells(cluster7, idents = c('19'))

Idents(cluster7) <- cluster7$predicted.id
pred16 <- WhichCells(cluster7, idents = '16')

int16 <- intersect(sub16, pred16) #932

#cluster 17 
Idents(cluster7) <- cluster7$seurat_clusters
sub17 <- WhichCells(cluster7, idents = c('14'))

Idents(cluster7) <- cluster7$predicted.id
pred17 <- WhichCells(cluster7, idents = '17')

int17 <- intersect(sub17, pred17) #1038

#cluster 18 
Idents(cluster7) <- cluster7$seurat_clusters
sub18 <- WhichCells(cluster7, idents = c('15'))

Idents(cluster7) <- cluster7$predicted.id
pred18 <- WhichCells(cluster7, idents = '18')

int18 <- intersect(sub18, pred18) #948

#cluster 19 
Idents(cluster7) <- cluster7$seurat_clusters
sub19 <- WhichCells(cluster7, idents = c('21'))

Idents(cluster7) <- cluster7$predicted.id
pred19 <- WhichCells(cluster7, idents = '19')

int19 <- intersect(sub19, pred19) #797

#cluster 20
Idents(cluster7) <- cluster7$seurat_clusters
sub20 <- WhichCells(cluster7, idents = c('17'))

Idents(cluster7) <- cluster7$predicted.id
pred20 <- WhichCells(cluster7, idents = '20')

int20 <- intersect(sub20, pred20) #949

#cluster 21
Idents(cluster7) <- cluster7$seurat_clusters
sub21 <- WhichCells(cluster7, idents = c('24'))

Idents(cluster7) <- cluster7$predicted.id
pred21 <- WhichCells(cluster7, idents = '21')

int21 <- intersect(sub21, pred21) #604

#cluster 22
Idents(cluster7) <- cluster7$seurat_clusters
sub22 <- WhichCells(cluster7, idents = c('27'))

Idents(cluster7) <- cluster7$predicted.id
pred22 <- WhichCells(cluster7, idents = '22')

int22 <- intersect(sub22, pred22) #503

#cluster 23
Idents(cluster7) <- cluster7$seurat_clusters
sub23 <- WhichCells(cluster7, idents = c('26'))

Idents(cluster7) <- cluster7$predicted.id
pred23 <- WhichCells(cluster7, idents = '23')

int23 <- intersect(sub23, pred23) #430

#cluster 24
Idents(cluster7) <- cluster7$seurat_clusters
sub24 <- WhichCells(cluster7, idents = c('28'))

Idents(cluster7) <- cluster7$predicted.id
pred24 <- WhichCells(cluster7, idents = '24')

int24 <- intersect(sub24, pred24) #461

#cluster 25
Idents(cluster7) <- cluster7$seurat_clusters
sub25 <- WhichCells(cluster7, idents = c('23'))

Idents(cluster7) <- cluster7$predicted.id
pred25 <- WhichCells(cluster7, idents = '25')

int25 <- intersect(sub25, pred25) #299

#cluster 26
Idents(cluster7) <- cluster7$seurat_clusters
sub26 <- WhichCells(cluster7, idents = c('26', '31'))

Idents(cluster7) <- cluster7$predicted.id
pred26 <- WhichCells(cluster7, idents = '26')

int26 <- intersect(sub26, pred26) #381

#cluster 27 
Idents(cluster7) <- cluster7$seurat_clusters
sub27 <- WhichCells(cluster7, idents = c('30'))

Idents(cluster7) <- cluster7$predicted.id
pred27 <- WhichCells(cluster7, idents = '27')

int27 <- intersect(sub27, pred27) #325

#cluster 28 
Idents(cluster7) <- cluster7$seurat_clusters
sub28 <- WhichCells(cluster7, idents = c('23'))

Idents(cluster7) <- cluster7$predicted.id
pred28 <- WhichCells(cluster7, idents = '28')

int28 <- intersect(sub28, pred28) #293

#cluster 29 
Idents(cluster7) <- cluster7$seurat_clusters
sub29 <- WhichCells(cluster7, idents = c('32'))

Idents(cluster7) <- cluster7$predicted.id
pred29 <- WhichCells(cluster7, idents = '29')

int29 <- intersect(sub29, pred29) #178

#cluster 30 
Idents(cluster7) <- cluster7$seurat_clusters
sub30 <- WhichCells(cluster7, idents = c('33'))

Idents(cluster7) <- cluster7$predicted.id
pred30 <- WhichCells(cluster7, idents = '30')

int30 <- intersect(sub30, pred30) #211

#cluster 31
Idents(cluster7) <- cluster7$seurat_clusters
sub31 <- WhichCells(cluster7, idents = c('8', '18'))

Idents(cluster7) <- cluster7$predicted.id
pred31 <- WhichCells(cluster7, idents = '31')

int31 <- intersect(sub31, pred31) #108

#cluster 32
Idents(cluster7) <- cluster7$seurat_clusters
sub32 <- WhichCells(cluster7, idents = c('34'))

Idents(cluster7) <- cluster7$predicted.id
pred32 <- WhichCells(cluster7, idents = '32')

int32 <- intersect(sub32, pred32) #95

#cluster 33
Idents(cluster7) <- cluster7$seurat_clusters
sub33 <- WhichCells(cluster7, idents = c('5', '11'))

Idents(cluster7) <- cluster7$predicted.id
pred33 <- WhichCells(cluster7, idents = '33')

int33 <- intersect(sub33, pred33) #16

#cluster 34
Idents(cluster7) <- cluster7$seurat_clusters
sub34 <- WhichCells(cluster7, idents = c('12'))

Idents(cluster7) <- cluster7$predicted.id
pred34 <- WhichCells(cluster7, idents = '34')

int34 <- intersect(sub34, pred34) #56

#cluster 35
Idents(cluster7) <- cluster7$seurat_clusters
sub35 <- WhichCells(cluster7, idents = c('36'))

Idents(cluster7) <- cluster7$predicted.id
pred35 <- WhichCells(cluster7, idents = '35')

int35 <- intersect(sub35, pred35) #32

#cluster 36
Idents(cluster7) <- cluster7$seurat_clusters
sub36 <- WhichCells(cluster7, idents = c('35'))

Idents(cluster7) <- cluster7$predicted.id
pred36 <- WhichCells(cluster7, idents = '36')

int36 <- intersect(sub36, pred36) #37

#cluster 37 
Idents(cluster7) <- cluster7$seurat_clusters
sub37 <- WhichCells(cluster7, idents = c('32'))

Idents(cluster7) <- cluster7$predicted.id
pred37 <- WhichCells(cluster7, idents = '37')

int37 <- intersect(sub37, pred37) #46

#cluster 38 
Idents(cluster7) <- cluster7$seurat_clusters
sub38 <- WhichCells(cluster7, idents = c('1', '17'))

Idents(cluster7) <- cluster7$predicted.id
pred38 <- WhichCells(cluster7, idents = '38')

int38 <- intersect(sub38, pred38) #14

#total 
#pre:43892
#39923
total7 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total7,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query7_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster7$query_label <- paste(cluster7$seurat_clusters, "query", sep = "_")
cluster7$ref_label <- paste(cluster7$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster7$ref_label), paste0(cluster7$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)


#-------Query 8---------#

cluster8 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub2_res15_4_high_res_cluster_pred.rds")

cluster8$predicted.id <- factor(cluster8$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster8, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster8, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table8 <- table(cluster8$predicted.id, cluster8$seurat_clusters)
table8 <- as.data.frame.matrix(table8)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table8_prop <- as.data.frame.matrix(apply(table8, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster8) <- cluster8$seurat_clusters
sub0 <- WhichCells(cluster8, idents = c('0'))

Idents(cluster8) <- cluster8$predicted.id
pred0 <- WhichCells(cluster8, idents = '0')

int0 <- intersect(sub0, pred0) #4258

#cluster 1
Idents(cluster8) <- cluster8$seurat_clusters
sub1 <- WhichCells(cluster8, idents = c('1'))

Idents(cluster8) <- cluster8$predicted.id
pred1 <- WhichCells(cluster8, idents = '1')

int1 <- intersect(sub1, pred1) #3039

#cluster 2
Idents(cluster8) <- cluster8$seurat_clusters
sub2 <- WhichCells(cluster8, idents = c('2'))

Idents(cluster8) <- cluster8$predicted.id
pred2 <- WhichCells(cluster8, idents = '2')

int2 <- intersect(sub2, pred2) #2326

#cluster 3
Idents(cluster8) <- cluster8$seurat_clusters
sub3 <- WhichCells(cluster8, idents = c('4', '6', '13'))

Idents(cluster8) <- cluster8$predicted.id
pred3 <- WhichCells(cluster8, idents = '3')

int3 <- intersect(sub3, pred3) #2389

#cluster 4
Idents(cluster8) <- cluster8$seurat_clusters
sub4 <- WhichCells(cluster8, idents = c('6', '20'))

Idents(cluster8) <- cluster8$predicted.id
pred4 <- WhichCells(cluster8, idents = '4')

int4 <- intersect(sub4, pred4) #2046

#cluster 5 
Idents(cluster8) <- cluster8$seurat_clusters
sub5 <- WhichCells(cluster8, idents = c('3'))

Idents(cluster8) <- cluster8$predicted.id
pred5 <- WhichCells(cluster8, idents = '5')

int5 <- intersect(sub5, pred5) #2299

#cluster 6
Idents(cluster8) <- cluster8$seurat_clusters
sub6 <- WhichCells(cluster8, idents = c('5', '8', '18'))

Idents(cluster8) <- cluster8$predicted.id
pred6 <- WhichCells(cluster8, idents = '6')

int6 <- intersect(sub6, pred6) #1846

#cluster 7 
Idents(cluster8) <- cluster8$seurat_clusters
sub7 <- WhichCells(cluster8, idents = c('4', '8'))

Idents(cluster8) <- cluster8$predicted.id
pred7 <- WhichCells(cluster8, idents = '7')

int7 <- intersect(sub7, pred7) #2171

#cluster 8 
Idents(cluster8) <- cluster8$seurat_clusters
sub8 <- WhichCells(cluster8, idents = c('12', '26'))

Idents(cluster8) <- cluster8$predicted.id
pred8 <- WhichCells(cluster8, idents = '8')

int8 <- intersect(sub8, pred8) #1816

#cluster 9 
Idents(cluster8) <- cluster8$seurat_clusters
sub9 <- WhichCells(cluster8, idents = c('7'))

Idents(cluster8) <- cluster8$predicted.id
pred9 <- WhichCells(cluster8, idents = '9')

int9 <- intersect(sub9, pred9) #1714

#cluster 10 
Idents(cluster8) <- cluster8$seurat_clusters
sub10 <- WhichCells(cluster8, idents = c('11'))

Idents(cluster8) <- cluster8$predicted.id
pred10 <- WhichCells(cluster8, idents = '10')

int10 <- intersect(sub10, pred10) #1464

#cluster 11
Idents(cluster8) <- cluster8$seurat_clusters
sub11 <- WhichCells(cluster8, idents = c('5', '23'))

Idents(cluster8) <- cluster8$predicted.id
pred11 <- WhichCells(cluster8, idents = '11')

int11 <- intersect(sub11, pred11) #1252

#cluster 12
Idents(cluster8) <- cluster8$seurat_clusters
sub12 <- WhichCells(cluster8, idents = c('10', '18'))

Idents(cluster8) <- cluster8$predicted.id
pred12 <- WhichCells(cluster8, idents = '12')

int12 <- intersect(sub12, pred12) #1262

#cluster 13
Idents(cluster8) <- cluster8$seurat_clusters
sub13 <- WhichCells(cluster8, idents = c('19', '23'))

Idents(cluster8) <- cluster8$predicted.id
pred13 <- WhichCells(cluster8, idents = '13')

int13 <- intersect(sub13, pred13) #1251

#cluster 14
Idents(cluster8) <- cluster8$seurat_clusters
sub14 <- WhichCells(cluster8, idents = c('9'))

Idents(cluster8) <- cluster8$predicted.id
pred14 <- WhichCells(cluster8, idents = '14')

int14 <- intersect(sub14, pred14) #1237

#cluster 15
Idents(cluster8) <- cluster8$seurat_clusters
sub15 <- WhichCells(cluster8, idents = c('15'))

Idents(cluster8) <- cluster8$predicted.id
pred15 <- WhichCells(cluster8, idents = '15')

int15 <- intersect(sub15, pred15) #1038

#cluster 16
Idents(cluster8) <- cluster8$seurat_clusters
sub16 <- WhichCells(cluster8, idents = c('16'))

Idents(cluster8) <- cluster8$predicted.id
pred16 <- WhichCells(cluster8, idents = '16')

int16 <- intersect(sub16, pred16) #971

#cluster 17 
Idents(cluster8) <- cluster8$seurat_clusters
sub17 <- WhichCells(cluster8, idents = c('17'))

Idents(cluster8) <- cluster8$predicted.id
pred17 <- WhichCells(cluster8, idents = '17')

int17 <- intersect(sub17, pred17) #967

#cluster 18 
Idents(cluster8) <- cluster8$seurat_clusters
sub18 <- WhichCells(cluster8, idents = c('14'))

Idents(cluster8) <- cluster8$predicted.id
pred18 <- WhichCells(cluster8, idents = '18')

int18 <- intersect(sub18, pred18) #913

#cluster 19 
Idents(cluster8) <- cluster8$seurat_clusters
sub19 <- WhichCells(cluster8, idents = c('22'))

Idents(cluster8) <- cluster8$predicted.id
pred19 <- WhichCells(cluster8, idents = '19')

int19 <- intersect(sub19, pred19) #774

#cluster 20 
Idents(cluster8) <- cluster8$seurat_clusters
sub20 <- WhichCells(cluster8, idents = c('21'))

Idents(cluster8) <- cluster8$predicted.id
pred20 <- WhichCells(cluster8, idents = '20')

int20 <- intersect(sub20, pred20) #908

#cluster 21
Idents(cluster8) <- cluster8$seurat_clusters
sub21 <- WhichCells(cluster8, idents = c('24'))

Idents(cluster8) <- cluster8$predicted.id
pred21 <- WhichCells(cluster8, idents = '21')

int21 <- intersect(sub21, pred21) #653

#cluster 22
Idents(cluster8) <- cluster8$seurat_clusters
sub22 <- WhichCells(cluster8, idents = c('25'))

Idents(cluster8) <- cluster8$predicted.id
pred22 <- WhichCells(cluster8, idents = '22')

int22 <- intersect(sub22, pred22) #563

#cluster 23
Idents(cluster8) <- cluster8$seurat_clusters
sub23 <- WhichCells(cluster8, idents = c('30'))

Idents(cluster8) <- cluster8$predicted.id
pred23 <- WhichCells(cluster8, idents = '23')

int23 <- intersect(sub23, pred23) #370

#cluster 24
Idents(cluster8) <- cluster8$seurat_clusters
sub24 <- WhichCells(cluster8, idents = c('28'))

Idents(cluster8) <- cluster8$predicted.id
pred24 <- WhichCells(cluster8, idents = '24')

int24 <- intersect(sub24, pred24) #463

#cluster 25
Idents(cluster8) <- cluster8$seurat_clusters
sub25 <- WhichCells(cluster8, idents = c('5', '32'))

Idents(cluster8) <- cluster8$predicted.id
pred25 <- WhichCells(cluster8, idents = '25')

int25 <- intersect(sub25, pred25) #356

#cluster 26
Idents(cluster8) <- cluster8$seurat_clusters
sub26 <- WhichCells(cluster8, idents = c('29'))

Idents(cluster8) <- cluster8$predicted.id
pred26 <- WhichCells(cluster8, idents = '26')

int26 <- intersect(sub26, pred26) #340

#cluster 27
Idents(cluster8) <- cluster8$seurat_clusters
sub27 <- WhichCells(cluster8, idents = c('31'))

Idents(cluster8) <- cluster8$predicted.id
pred27 <- WhichCells(cluster8, idents = '27')

int27 <- intersect(sub27, pred27) #304

#cluster 28 
Idents(cluster8) <- cluster8$seurat_clusters
sub28 <- WhichCells(cluster8, idents = c('27'))

Idents(cluster8) <- cluster8$predicted.id
pred28 <- WhichCells(cluster8, idents = '28')

int28 <- intersect(sub28, pred28) #217

#cluster 29 
Idents(cluster8) <- cluster8$seurat_clusters
sub29 <- WhichCells(cluster8, idents = c('33'))

Idents(cluster8) <- cluster8$predicted.id
pred29 <- WhichCells(cluster8, idents = '29')

int29 <- intersect(sub29, pred29) #217

#cluster 30 
Idents(cluster8) <- cluster8$seurat_clusters
sub30 <- WhichCells(cluster8, idents = c('34'))

Idents(cluster8) <- cluster8$predicted.id
pred30 <- WhichCells(cluster8, idents = '30')

int30 <- intersect(sub30, pred30) #190

#cluster 31
Idents(cluster8) <- cluster8$seurat_clusters
sub31 <- WhichCells(cluster8, idents = c('27', '37'))

Idents(cluster8) <- cluster8$predicted.id
pred31 <- WhichCells(cluster8, idents = '31')

int31 <- intersect(sub31, pred31) #113

#cluster 32
Idents(cluster8) <- cluster8$seurat_clusters
sub32 <- WhichCells(cluster8, idents = c('36'))

Idents(cluster8) <- cluster8$predicted.id
pred32 <- WhichCells(cluster8, idents = '32')

int32 <- intersect(sub32, pred32) #85

#cluster 33
Idents(cluster8) <- cluster8$seurat_clusters
sub33 <- WhichCells(cluster8, idents = c('35'))

Idents(cluster8) <- cluster8$predicted.id
pred33 <- WhichCells(cluster8, idents = '33')

int33 <- intersect(sub33, pred33) #82

#cluster 34
Idents(cluster8) <- cluster8$seurat_clusters
sub34 <- WhichCells(cluster8, idents = c('9'))

Idents(cluster8) <- cluster8$predicted.id
pred34 <- WhichCells(cluster8, idents = '34')

int34 <- intersect(sub34, pred34) #48

#cluster 35
Idents(cluster8) <- cluster8$seurat_clusters
sub35 <- WhichCells(cluster8, idents = c('38'))

Idents(cluster8) <- cluster8$predicted.id
pred35 <- WhichCells(cluster8, idents = '35')

int35 <- intersect(sub35, pred35) #50

#cluster 36
Idents(cluster8) <- cluster8$seurat_clusters
sub36 <- WhichCells(cluster8, idents = c('12'))

Idents(cluster8) <- cluster8$predicted.id
pred36 <- WhichCells(cluster8, idents = '36')

int36 <- intersect(sub36, pred36) #57

#cluster 37
Idents(cluster8) <- cluster8$seurat_clusters
sub37 <- WhichCells(cluster8, idents = c('39'))

Idents(cluster8) <- cluster8$predicted.id
pred37 <- WhichCells(cluster8, idents = '37')

int37 <- intersect(sub37, pred37) #45

#cluster 38 
Idents(cluster8) <- cluster8$seurat_clusters
sub38 <- WhichCells(cluster8, idents = c('16'))

Idents(cluster8) <- cluster8$predicted.id
pred38 <- WhichCells(cluster8, idents = '38')

int38 <- intersect(sub38, pred38) #24

#total 
#pre:43893
#40118
total8 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total8,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query8_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster8$query_label <- paste(cluster8$seurat_clusters, "query", sep = "_")
cluster8$ref_label <- paste(cluster8$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster8$ref_label), paste0(cluster8$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-------Query 9---------#

cluster9 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub_res15_5_high_res_cluster_pred.rds")

cluster9$predicted.id <- factor(cluster9$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster9, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster9, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table9 <- table(cluster9$predicted.id, cluster9$seurat_clusters)
table9 <- as.data.frame.matrix(table9)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table9_prop <- as.data.frame.matrix(apply(table9, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster9) <- cluster9$seurat_clusters
sub0 <- WhichCells(cluster9, idents = c('0'))

Idents(cluster9) <- cluster9$predicted.id
pred0 <- WhichCells(cluster9, idents = '0')

int0 <- intersect(sub0, pred0) #4247

#cluster 1
Idents(cluster9) <- cluster9$seurat_clusters
sub1 <- WhichCells(cluster9, idents = c('1'))

Idents(cluster9) <- cluster9$predicted.id
pred1 <- WhichCells(cluster9, idents = '1')

int1 <- intersect(sub1, pred1) #3268

#cluster 2
Idents(cluster9) <- cluster9$seurat_clusters
sub2 <- WhichCells(cluster9, idents = c('2'))

Idents(cluster9) <- cluster9$predicted.id
pred2 <- WhichCells(cluster9, idents = '2')

int2 <- intersect(sub2, pred2) #2326

#cluster 3
Idents(cluster9) <- cluster9$seurat_clusters
sub3 <- WhichCells(cluster9, idents = c('4', '13', '15'))

Idents(cluster9) <- cluster9$predicted.id
pred3 <- WhichCells(cluster9, idents = '3')

int3 <- intersect(sub3, pred3) #2200

#cluster 4
Idents(cluster9) <- cluster9$seurat_clusters
sub4 <- WhichCells(cluster9, idents = c('5', '13'))

Idents(cluster9) <- cluster9$predicted.id
pred4 <- WhichCells(cluster9, idents = '4')

int4 <- intersect(sub4, pred4) #2068

#cluster 5
Idents(cluster9) <- cluster9$seurat_clusters
sub5 <- WhichCells(cluster9, idents = c('3'))

Idents(cluster9) <- cluster9$predicted.id
pred5 <- WhichCells(cluster9, idents = '5')

int5 <- intersect(sub5, pred5) #2352

#cluster 6
Idents(cluster9) <- cluster9$seurat_clusters
sub6 <- WhichCells(cluster9, idents = c('2', '6', '7'))

Idents(cluster9) <- cluster9$predicted.id
pred6 <- WhichCells(cluster9, idents = '6')

int6 <- intersect(sub6, pred6) #1818

#cluster 7
Idents(cluster9) <- cluster9$seurat_clusters
sub7 <- WhichCells(cluster9, idents = c('4'))

Idents(cluster9) <- cluster9$predicted.id
pred7 <- WhichCells(cluster9, idents = '7')

int7 <- intersect(sub7, pred7) #1757

#cluster 8
Idents(cluster9) <- cluster9$seurat_clusters
sub8 <- WhichCells(cluster9, idents = c('12', '27'))

Idents(cluster9) <- cluster9$predicted.id
pred8 <- WhichCells(cluster9, idents = '8')

int8 <- intersect(sub8, pred8) #1811

#cluster 9 
Idents(cluster9) <- cluster9$seurat_clusters
sub9 <- WhichCells(cluster9, idents = c('8'))

Idents(cluster9) <- cluster9$predicted.id
pred9 <- WhichCells(cluster9, idents = '9')

int9 <- intersect(sub9, pred9) #1568

#cluster 10 
Idents(cluster9) <- cluster9$seurat_clusters
sub10 <- WhichCells(cluster9, idents = c('10'))

Idents(cluster9) <- cluster9$predicted.id
pred10 <- WhichCells(cluster9, idents = '10')

int10 <- intersect(sub10, pred10) #1413

#cluster 11
Idents(cluster9) <- cluster9$seurat_clusters
sub11 <- WhichCells(cluster9, idents = c('16', '20'))

Idents(cluster9) <- cluster9$predicted.id
pred11 <- WhichCells(cluster9, idents = '11')

int11 <- intersect(sub11, pred11) #1207

#cluster 12
Idents(cluster9) <- cluster9$seurat_clusters
sub12 <- WhichCells(cluster9, idents = c('6'))

Idents(cluster9) <- cluster9$predicted.id
pred12 <- WhichCells(cluster9, idents = '12')

int12 <- intersect(sub12, pred12) #1136

#cluster 13
Idents(cluster9) <- cluster9$seurat_clusters
sub13 <- WhichCells(cluster9, idents = c('11'))

Idents(cluster9) <- cluster9$predicted.id
pred13 <- WhichCells(cluster9, idents = '13')

int13 <- intersect(sub13, pred13) #1173

#cluster 14
Idents(cluster9) <- cluster9$seurat_clusters
sub14 <- WhichCells(cluster9, idents = c('9'))

Idents(cluster9) <- cluster9$predicted.id
pred14 <- WhichCells(cluster9, idents = '14')

int14 <- intersect(sub14, pred14) #1205

#cluster 15
Idents(cluster9) <- cluster9$seurat_clusters
sub15 <- WhichCells(cluster9, idents = c('23', '24'))

Idents(cluster9) <- cluster9$predicted.id
pred15 <- WhichCells(cluster9, idents = '15')

int15 <- intersect(sub15, pred15) #1244

#cluster 16
Idents(cluster9) <- cluster9$seurat_clusters
sub16 <- WhichCells(cluster9, idents = c('19'))

Idents(cluster9) <- cluster9$predicted.id
pred16 <- WhichCells(cluster9, idents = '16')

int16 <- intersect(sub16, pred16) #902

#cluster 17
Idents(cluster9) <- cluster9$seurat_clusters
sub17 <- WhichCells(cluster9, idents = c('14'))

Idents(cluster9) <- cluster9$predicted.id
pred17 <- WhichCells(cluster9, idents = '17')

int17 <- intersect(sub17, pred17) #1002

#cluster 18
Idents(cluster9) <- cluster9$seurat_clusters
sub18 <- WhichCells(cluster9, idents = c('18', '24'))

Idents(cluster9) <- cluster9$predicted.id
pred18 <- WhichCells(cluster9, idents = '18')

int18 <- intersect(sub18, pred18) #965

#cluster 19 
Idents(cluster9) <- cluster9$seurat_clusters
sub19 <- WhichCells(cluster9, idents = c('21'))

Idents(cluster9) <- cluster9$predicted.id
pred19 <- WhichCells(cluster9, idents = '19')

int19 <- intersect(sub19, pred19) #842

#cluster 20 
Idents(cluster9) <- cluster9$seurat_clusters
sub20 <- WhichCells(cluster9, idents = c('17'))

Idents(cluster9) <- cluster9$predicted.id
pred20 <- WhichCells(cluster9, idents = '20')

int20 <- intersect(sub20, pred20) #895

#cluster 21
Idents(cluster9) <- cluster9$seurat_clusters
sub21 <- WhichCells(cluster9, idents = c('25'))

Idents(cluster9) <- cluster9$predicted.id
pred21 <- WhichCells(cluster9, idents = '21')

int21 <- intersect(sub21, pred21) #647

#cluster 22
Idents(cluster9) <- cluster9$seurat_clusters
sub22 <- WhichCells(cluster9, idents = c('26'))

Idents(cluster9) <- cluster9$predicted.id
pred22 <- WhichCells(cluster9, idents = '22')

int22 <- intersect(sub22, pred22) #571


#cluster 23
Idents(cluster9) <- cluster9$seurat_clusters
sub23 <- WhichCells(cluster9, idents = c('29'))

Idents(cluster9) <- cluster9$predicted.id
pred23 <- WhichCells(cluster9, idents = '23')

int23 <- intersect(sub23, pred23) #423

#cluster 24
Idents(cluster9) <- cluster9$seurat_clusters
sub24 <- WhichCells(cluster9, idents = c('28'))

Idents(cluster9) <- cluster9$predicted.id
pred24 <- WhichCells(cluster9, idents = '24')

int24 <- intersect(sub24, pred24) #486

#cluster 25
Idents(cluster9) <- cluster9$seurat_clusters
sub25 <- WhichCells(cluster9, idents = c('22'))

Idents(cluster9) <- cluster9$predicted.id
pred25 <- WhichCells(cluster9, idents = '25')

int25 <- intersect(sub25, pred25) #350

#cluster 26
Idents(cluster9) <- cluster9$seurat_clusters
sub26 <- WhichCells(cluster9, idents = c('30'))

Idents(cluster9) <- cluster9$predicted.id
pred26 <- WhichCells(cluster9, idents = '26')

int26 <- intersect(sub26, pred26) #330

#cluster 27 
Idents(cluster9) <- cluster9$seurat_clusters
sub27 <- WhichCells(cluster9, idents = c('31'))

Idents(cluster9) <- cluster9$predicted.id
pred27 <- WhichCells(cluster9, idents = '27')

int27 <- intersect(sub27, pred27) #309

#cluster 28 
Idents(cluster9) <- cluster9$seurat_clusters
sub28 <- WhichCells(cluster9, idents = c('22'))

Idents(cluster9) <- cluster9$predicted.id
pred28 <- WhichCells(cluster9, idents = '28')

int28 <- intersect(sub28, pred28) #252

#cluster 29 
Idents(cluster9) <- cluster9$seurat_clusters
sub29 <- WhichCells(cluster9, idents = c('32'))

Idents(cluster9) <- cluster9$predicted.id
pred29 <- WhichCells(cluster9, idents = '29')

int29 <- intersect(sub29, pred29) #218

#cluster 30 
Idents(cluster9) <- cluster9$seurat_clusters
sub30 <- WhichCells(cluster9, idents = c('33'))

Idents(cluster9) <- cluster9$predicted.id
pred30 <- WhichCells(cluster9, idents = '30')

int30 <- intersect(sub30, pred30) #190

#cluster 31
Idents(cluster9) <- cluster9$seurat_clusters
sub31 <- WhichCells(cluster9, idents = c('7', '34'))

Idents(cluster9) <- cluster9$predicted.id
pred31 <- WhichCells(cluster9, idents = '31')

int31 <- intersect(sub31, pred31) #101

#cluster 32
Idents(cluster9) <- cluster9$seurat_clusters
sub32 <- WhichCells(cluster9, idents = c('35'))

Idents(cluster9) <- cluster9$predicted.id
pred32 <- WhichCells(cluster9, idents = '32')

int32 <- intersect(sub32, pred32) #86

#cluster 33
Idents(cluster9) <- cluster9$seurat_clusters
sub33 <- WhichCells(cluster9, idents = c('16'))

Idents(cluster9) <- cluster9$predicted.id
pred33 <- WhichCells(cluster9, idents = '33')

int33 <- intersect(sub33, pred33) #66

#cluster 34
Idents(cluster9) <- cluster9$seurat_clusters
sub34 <- WhichCells(cluster9, idents = c('36'))

Idents(cluster9) <- cluster9$predicted.id
pred34 <- WhichCells(cluster9, idents = '34')

int34 <- intersect(sub34, pred34) #60

#cluster 35
Idents(cluster9) <- cluster9$seurat_clusters
sub35 <- WhichCells(cluster9, idents = c('38'))

Idents(cluster9) <- cluster9$predicted.id
pred35 <- WhichCells(cluster9, idents = '35')

int35 <- intersect(sub35, pred35) #36

#cluster 36
Idents(cluster9) <- cluster9$seurat_clusters
sub36 <- WhichCells(cluster9, idents = c('39'))

Idents(cluster9) <- cluster9$predicted.id
pred36 <- WhichCells(cluster9, idents = '36')

int36 <- intersect(sub36, pred36) #33

#cluster 37
Idents(cluster9) <- cluster9$seurat_clusters
sub37 <- WhichCells(cluster9, idents = c('37'))

Idents(cluster9) <- cluster9$predicted.id
pred37 <- WhichCells(cluster9, idents = '37')

int37 <- intersect(sub37, pred37) #54

#cluster 38
Idents(cluster9) <- cluster9$seurat_clusters
sub38 <- WhichCells(cluster9, idents = c('19'))

Idents(cluster9) <- cluster9$predicted.id
pred38 <- WhichCells(cluster9, idents = '38')

int38 <- intersect(sub38, pred38) #19

#total 
#pre:43892
#39630
total9 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total9,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query9_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster9$query_label <- paste(cluster9$seurat_clusters, "query", sep = "_")
cluster9$ref_label <- paste(cluster9$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster9$ref_label), paste0(cluster9$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-------Query 10---------#

cluster10 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/ileum_sub2_res15_5_high_res_cluster_pred.rds")

cluster10$predicted.id <- factor(cluster10$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38"))

#visualize the umap 
DimPlot(cluster10, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) + DimPlot(cluster10, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

#create table that has the number of cells in each cluster
table10 <- table(cluster10$predicted.id, cluster10$seurat_clusters)
table10 <- as.data.frame.matrix(table10)

#seurat cluster is the x axis and predicted id is the y axis 

#proportion of cells in each cluster 
table10_prop <- as.data.frame.matrix(apply(table10, 1, function(x) x/sum(x)))

#seurat cluster is the y axis and predicted id is the x axis - using this for identifying the cluster overlap

#rules - take the largest proportions that add up to at least 80% 
#cluster0 

Idents(cluster10) <- cluster10$seurat_clusters
sub0 <- WhichCells(cluster10, idents = c('0'))

Idents(cluster10) <- cluster10$predicted.id
pred0 <- WhichCells(cluster10, idents = '0')

int0 <- intersect(sub0, pred0) #4867

#cluster 1
Idents(cluster10) <- cluster10$seurat_clusters
sub1 <- WhichCells(cluster10, idents = c('4', '13'))

Idents(cluster10) <- cluster10$predicted.id
pred1 <- WhichCells(cluster10, idents = '1')

int1 <- intersect(sub1, pred1) #3226

#cluster 2
Idents(cluster10) <- cluster10$seurat_clusters
sub2 <- WhichCells(cluster10, idents = c('1'))

Idents(cluster10) <- cluster10$predicted.id
pred2 <- WhichCells(cluster10, idents = '2')

int2 <- intersect(sub2, pred2) #2548

#cluster 3
Idents(cluster10) <- cluster10$seurat_clusters
sub3 <- WhichCells(cluster10, idents = c('2', '6'))

Idents(cluster10) <- cluster10$predicted.id
pred3 <- WhichCells(cluster10, idents = '3')

int3 <- intersect(sub3, pred3) #2191

#cluster 4
Idents(cluster10) <- cluster10$seurat_clusters
sub4 <- WhichCells(cluster10, idents = c('2'))

Idents(cluster10) <- cluster10$predicted.id
pred4 <- WhichCells(cluster10, idents = '4')

int4 <- intersect(sub4, pred4) #2125

#cluster 5
Idents(cluster10) <- cluster10$seurat_clusters
sub5 <- WhichCells(cluster10, idents = c('3'))

Idents(cluster10) <- cluster10$predicted.id
pred5 <- WhichCells(cluster10, idents = '5')

int5 <- intersect(sub5, pred5) #2056

#cluster 6
Idents(cluster10) <- cluster10$seurat_clusters
sub6 <- WhichCells(cluster10, idents = c('9', '11'))

Idents(cluster10) <- cluster10$predicted.id
pred6 <- WhichCells(cluster10, idents = '6')

int6 <- intersect(sub6, pred6) #1896

#cluster 7 
Idents(cluster10) <- cluster10$seurat_clusters
sub7 <- WhichCells(cluster10, idents = c('7'))

Idents(cluster10) <- cluster10$predicted.id
pred7 <- WhichCells(cluster10, idents = '7')

int7 <- intersect(sub7, pred7) #1710

#cluster 8 
Idents(cluster10) <- cluster10$seurat_clusters
sub8 <- WhichCells(cluster10, idents = c('14', '22'))

Idents(cluster10) <- cluster10$predicted.id
pred8 <- WhichCells(cluster10, idents = '8')

int8 <- intersect(sub8, pred8) #1871

#cluster 9 
Idents(cluster10) <- cluster10$seurat_clusters
sub9 <- WhichCells(cluster10, idents = c('5'))

Idents(cluster10) <- cluster10$predicted.id
pred9 <- WhichCells(cluster10, idents = '9')

int9 <- intersect(sub9, pred9) #1916

#cluster 10 
Idents(cluster10) <- cluster10$seurat_clusters
sub10 <- WhichCells(cluster10, idents = c('10'))

Idents(cluster10) <- cluster10$predicted.id
pred10 <- WhichCells(cluster10, idents = '10')

int10 <- intersect(sub10, pred10) #1443

#cluster 11
Idents(cluster10) <- cluster10$seurat_clusters
sub11<- WhichCells(cluster10, idents = c('11', '21'))

Idents(cluster10) <- cluster10$predicted.id
pred11 <- WhichCells(cluster10, idents = '11')

int11 <- intersect(sub11, pred11) #1137

#cluster 12
Idents(cluster10) <- cluster10$seurat_clusters
sub12 <- WhichCells(cluster10, idents = c('4', '17'))

Idents(cluster10) <- cluster10$predicted.id
pred12 <- WhichCells(cluster10, idents = '12')

int12 <- intersect(sub12, pred12) #1104

#cluster 13
Idents(cluster10) <- cluster10$seurat_clusters
sub13 <- WhichCells(cluster10, idents = c('12'))

Idents(cluster10) <- cluster10$predicted.id
pred13 <- WhichCells(cluster10, idents = '13')

int13 <- intersect(sub13, pred13) #1171

#cluster 14
Idents(cluster10) <- cluster10$seurat_clusters
sub14 <- WhichCells(cluster10, idents = c('15'))

Idents(cluster10) <- cluster10$predicted.id
pred14 <- WhichCells(cluster10, idents = '14')

int14 <- intersect(sub14, pred14) #1180

#cluster 15
Idents(cluster10) <- cluster10$seurat_clusters
sub15 <- WhichCells(cluster10, idents = c('16'))

Idents(cluster10) <- cluster10$predicted.id
pred15 <- WhichCells(cluster10, idents = '15')

int15 <- intersect(sub15, pred15) #1142

#cluster 16 
Idents(cluster10) <- cluster10$seurat_clusters
sub16 <- WhichCells(cluster10, idents = c('8'))

Idents(cluster10) <- cluster10$predicted.id
pred16 <- WhichCells(cluster10, idents = '16')

int16 <- intersect(sub16, pred16) #965

#cluster 17
Idents(cluster10) <- cluster10$seurat_clusters
sub17 <- WhichCells(cluster10, idents = c('3', '25'))

Idents(cluster10) <- cluster10$predicted.id
pred17 <- WhichCells(cluster10, idents = '17')

int17 <- intersect(sub17, pred17) #1043

#cluster 18 
Idents(cluster10) <- cluster10$seurat_clusters
sub18 <- WhichCells(cluster10, idents = c('18'))

Idents(cluster10) <- cluster10$predicted.id
pred18 <- WhichCells(cluster10, idents = '18')

int18 <- intersect(sub18, pred18) #1007

#cluster 19 
Idents(cluster10) <- cluster10$seurat_clusters
sub19 <- WhichCells(cluster10, idents = c('19'))

Idents(cluster10) <- cluster10$predicted.id
pred19 <- WhichCells(cluster10, idents = '19')

int19 <- intersect(sub19, pred19) #710

#cluster 20 
Idents(cluster10) <- cluster10$seurat_clusters
sub20 <- WhichCells(cluster10, idents = c('8', '27'))

Idents(cluster10) <- cluster10$predicted.id
pred20 <- WhichCells(cluster10, idents = '20')

int20 <- intersect(sub20, pred20) #960

#cluster 21
Idents(cluster10) <- cluster10$seurat_clusters
sub21 <- WhichCells(cluster10, idents = c('23'))

Idents(cluster10) <- cluster10$predicted.id
pred21 <- WhichCells(cluster10, idents = '21')

int21 <- intersect(sub21, pred21) #627

#cluster 22
Idents(cluster10) <- cluster10$seurat_clusters
sub22 <- WhichCells(cluster10, idents = c('20'))

Idents(cluster10) <- cluster10$predicted.id
pred22 <- WhichCells(cluster10, idents = '22')

int22 <- intersect(sub22, pred22) #568

#cluster 23
Idents(cluster10) <- cluster10$seurat_clusters
sub23 <- WhichCells(cluster10, idents = c('26'))

Idents(cluster10) <- cluster10$predicted.id
pred23 <- WhichCells(cluster10, idents = '23')

int23 <- intersect(sub23, pred23) #456

#cluster 24
Idents(cluster10) <- cluster10$seurat_clusters
sub24 <- WhichCells(cluster10, idents = c('28'))

Idents(cluster10) <- cluster10$predicted.id
pred24 <- WhichCells(cluster10, idents = '24')

int24 <- intersect(sub24, pred24) #409

#cluster 25
Idents(cluster10) <- cluster10$seurat_clusters
sub25 <- WhichCells(cluster10, idents = c('11', '24'))

Idents(cluster10) <- cluster10$predicted.id
pred25 <- WhichCells(cluster10, idents = '25')

int25 <- intersect(sub25, pred25) #392

#cluster 26
Idents(cluster10) <- cluster10$seurat_clusters
sub26 <- WhichCells(cluster10, idents = c('26', '31'))

Idents(cluster10) <- cluster10$predicted.id
pred26 <- WhichCells(cluster10, idents = '26')

int26 <- intersect(sub26, pred26) #349

#cluster 27
Idents(cluster10) <- cluster10$seurat_clusters
sub27 <- WhichCells(cluster10, idents = c('30'))

Idents(cluster10) <- cluster10$predicted.id
pred27 <- WhichCells(cluster10, idents = '27')

int27 <- intersect(sub27, pred27) #323

#clulster 28 
Idents(cluster10) <- cluster10$seurat_clusters
sub28 <- WhichCells(cluster10, idents = c('24'))

Idents(cluster10) <- cluster10$predicted.id
pred28 <- WhichCells(cluster10, idents = '28')

int28 <- intersect(sub28, pred28) #231

#cluster 29 
Idents(cluster10) <- cluster10$seurat_clusters
sub29 <- WhichCells(cluster10, idents = c('33'))

Idents(cluster10) <- cluster10$predicted.id
pred29 <- WhichCells(cluster10, idents = '29')

int29 <- intersect(sub29, pred29) #177

#cluster 30 
Idents(cluster10) <- cluster10$seurat_clusters
sub30 <- WhichCells(cluster10, idents = c('32'))

Idents(cluster10) <- cluster10$predicted.id
pred30 <- WhichCells(cluster10, idents = '30')

int30 <- intersect(sub30, pred30) #205

#cluster 31
Idents(cluster10) <- cluster10$seurat_clusters
sub31<- WhichCells(cluster10, idents = c('34'))

Idents(cluster10) <- cluster10$predicted.id
pred31 <- WhichCells(cluster10, idents = '31')

int31 <- intersect(sub31, pred31) #101

#cluster 32
Idents(cluster10) <- cluster10$seurat_clusters
sub32 <- WhichCells(cluster10, idents = c('35'))

Idents(cluster10) <- cluster10$predicted.id
pred32 <- WhichCells(cluster10, idents = '32')

int32 <- intersect(sub32, pred32) #83

#cluster 33
Idents(cluster10) <- cluster10$seurat_clusters
sub33 <- WhichCells(cluster10, idents = c('36'))

Idents(cluster10) <- cluster10$predicted.id
pred33 <- WhichCells(cluster10, idents = '33')

int33 <- intersect(sub33, pred33) #79

#cluster 34
Idents(cluster10) <- cluster10$seurat_clusters
sub34 <- WhichCells(cluster10, idents = c('5'))

Idents(cluster10) <- cluster10$predicted.id
pred34 <- WhichCells(cluster10, idents = '34')

int34 <- intersect(sub34, pred34) #46

#cluster 35
Idents(cluster10) <- cluster10$seurat_clusters
sub35 <- WhichCells(cluster10, idents = c('37'))

Idents(cluster10) <- cluster10$predicted.id
pred35 <- WhichCells(cluster10, idents = '35')

int35 <- intersect(sub35, pred35) #44

#cluster 36
Idents(cluster10) <- cluster10$seurat_clusters
sub36 <- WhichCells(cluster10, idents = c('20'))

Idents(cluster10) <- cluster10$predicted.id
pred36 <- WhichCells(cluster10, idents = '36')

int36 <- intersect(sub36, pred36) #57

#cluster 37 
Idents(cluster10) <- cluster10$seurat_clusters
sub37 <- WhichCells(cluster10, idents = c('11'))

Idents(cluster10) <- cluster10$predicted.id
pred37 <- WhichCells(cluster10, idents = '37')

int37 <- intersect(sub37, pred37) #36

#cluster 38 
Idents(cluster10) <- cluster10$seurat_clusters
sub38 <- WhichCells(cluster10, idents = c('13'))

Idents(cluster10) <- cluster10$predicted.id
pred38 <- WhichCells(cluster10, idents = '38')

int38 <- intersect(sub38, pred38) #16

#total 
#pre:43893
#40467
total10 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15, int16, int17, int18, int19, int20, int21, int22, int23, int24, int25, int26, int27, int28, int29, int30, int31, int32, int33, int34, int35, int36, int37, int38))
write.csv(total10,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query10_high_res_barcodes_10_25_24.csv', row.names=F)

#create confusion matrix heatmap 
cluster10$query_label <- paste(cluster10$seurat_clusters, "query", sep = "_")
cluster10$ref_label <- paste(cluster10$predicted.id, "ref", sep = "_")
library(ArchR)
cm1 <- confusionMatrix(paste0(cluster10$ref_label), paste0(cluster10$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)
cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
)

#-----load in the stable barcodes----------#
#load in barcodes 
total1 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query1_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)
total2 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query2_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)
total3 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query3_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)
total4 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query4_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)
total5 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query5_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)
total6 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query6_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)
total7 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query7_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)
total8 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query8_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)
total9 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query9_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)
total10 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/query10_high_res_barcodes_10_25_24.csv', what = "", sep = ",", skip = 1)

#merge the barcodes into a list 
#400171 barcodes 
total_merge <- c(total1, total2, total3, total4, total5, total6, total7, total8, total9, total10)

#counts the number of times each barcode is present
value_counts <- table(total_merge)
table(value_counts)

#1: 1202
#2: 2727
#3: 5915
#4: 12425
#5: 65214

#define a threshold - more times a barcode is present, more stable it is
threshold <- 4

#filter values that meet your threshold 
filtered_values <- names(value_counts[value_counts >= threshold])  

#1: 87438 - keeps ~99.9% of data
#2: 86281 - keeps ~98.28% of data 
#3: 83554 - keeps ~95.18% of data
#4: 77639 - keeps ~88.44% of data ----select threshold of 4 (biggest jump after threshold of 4)
#5: 65214 - keeps ~74.28% of data

#save the filtered barcodes 
write.csv(filtered_values, file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/barcodes/filtered_high_res_barcodes_10_28_24.csv', row.names = F)

#subset the ileum seurat object on the server 

#load SO
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_rpca_res1_10_21_24.rds")

#load in stably assigned cells 
filtered_values <- scan('/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/filtered_high_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)

#subset ileum based on stable barcodes 
ileum_sub <- subset(ileum, cells = filtered_values)

#save 
saveRDS(ileum_sub, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_13_ileum_high_res_stable_cell_10_28_24.rds")

#----------highly stably assigned and unstably assigned cells high resolution-------------#

#10/28/2024

#load in ileum SO 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_rPCA_10_17_24.rds")

#ileum barcodes 
barcodes <- colnames(ileum)

unstable <- setdiff(barcodes, filtered_values)

DimPlot(ileum, reduction = "umap", cells.highlight = unstable, sizes.highlight = 0.1) + scale_color_manual(labels = c("Stably Assigned Cells", "Unstably Assigned Cells"), values = c("grey", "red"))

#--------Assess ileum high res ARI score------------#

#merged ARI scores
df <- read.csv('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/df_ri_helm_batch1_13_ileum_high_res_dune_merge_10_29_24.csv', row.names = NULL, sep = ",")

df <- df[,-1]
names(df) <- gsub("^paste0\\.", "", names(df))

columns <- colnames(df)

#plot the ARI - without paste0
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/figures/RI_helm_batch1_13_ileum_high_res_all_param_merge_10_29_24.pdf", width = 15, height = 15)
plot(plotARIs(df %>% select(columns)) + theme(text=element_text(size=8)) + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) + theme(axis.text.y = element_text(size = 10)))
dev.off()


#----mean and median----#
#ARI values 
plot <- plotARIs(clusMat = df)

ari_values <- data.frame(plot[["plot_env"]][["Mat"]])

#calculate mean ARI
ari_values$mean_ari <- rowMeans(ari_values)

#top 10 mean ari 
mean_ari <- as.data.frame(ari_values$mean_ari)
rownames(mean_ari) <- rownames(ari_values)
top_n(mean_ari, 10)

#overall mean ari
summarise(ari_values, overall_mean = mean(mean_ari)) #0.8145946

#remove mean_ari from ari_values 
ari_values <- ari_values[, -c(37)]

#median ARI

ari_values$row_median = apply(ari_values, 1, median, na.rm=TRUE)

overall_median <- median(ari_values$row_median)  #0.8330813


#take average of median 
summarise(ari_values, overall_mean = mean(row_median)) #0.8160038

#top 10 median ARI scores 
median_ari <- as.data.frame(ari_values$row_median)
rownames(median_ari) <- rownames(ari_values)
top_n(median_ari, 10)

#------------Stably Assigned Cells Clustering results-------------#

#10/31/2024

#-------cluster 1---------#

#parameters: res: 0.5, PC: 15, HVG: 2000

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_50_15_2000_rpca_10_29_24.rds')

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24'))

library(randomcoloR)
no_of_colors <- 25

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#run FM on cluster 24
cluster24 <- FindMarkers(ileum, ident.1 = "24")

#----------preliminary cluster annotation for stability assessment comparison-------#

#11/17/2024

ileum_low@meta.data$cell_typev1 <- ileum_low@meta.data$seurat_clusters
ileum_low$cell_typev1 <- plyr::mapvalues(
  x = ileum_low$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24'),
  to = c("Naive B Cell", 'Stem/Paneth Cell', 'Enterocyte 1', 'CD4 T Cell', 'Enterocyte 2', 'Macrophage', "CD8 T Cell", "Plasma Cell 1", "Enterocyte 3", "Plasma Cell 2", "Monocyte", "Goblet Cell 1", "Goblet Cell 2", "Ambiguous T Cell", "Mesenchymal Cell", "NK/T Cell", "Inflammatory Monocyte", "TA Cell", "Memory B Cell", 'Endothelial Cell', "Mast Cell", "Tuft Cell", "Doublet Cluster", "ILC3", "Plasmablast")
)


ileum_low@meta.data$cell_typev2 <- ileum_low@meta.data$seurat_clusters
ileum_low$cell_typev2 <- plyr::mapvalues(
  x = ileum_low$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24'),
  to = c("Naive B Cell", 'Stem/Paneth Cell', 'Enterocyte', 'CD4 T Cell', 'Enterocyte', 'Macrophage', "CD8 T Cell", "Plasma Cell", "Enterocyte", "Plasma Cell", "Monocyte", "Goblet Cell", "Goblet Cell", "Ambiguous T Cell", "Mesenchymal Cell", "NK/T Cell", "Inflammatory Monocyte", "TA Cell", "Memory B Cell", 'Endothelial Cell', "Mast Cell", "Tuft Cell", "Doublet Cluster", "ILC3", "Plasmablast")
)

#colors for UMAP and stacked bar plot 

#blue = t cells 
blue <- randomColor(8, hue = c("blue"))
blue <- c("#34bfe5", "#4a93ce", "#b0daf2", "#6db5ed", "#1c5682","#077b9e", "#70c8e0", "#1f86f4") 

#green = b cells
green <- randomColor(8, hue = c("green"))
green <- c("#05fc95", "#9de587", "#1da54f", "#61e291", "#2bd147","#80cc1e", "#c5ef99", "#048c45")


#orange = myeloid 
orange <- randomColor(5, hue = c("orange"))
orange <- c("#fcab85", "#d3b050", "#f75f13", "#cc5414", "#f98a52")

#pink = epithelial
pink <- randomColor(12, hue = c("pink"), luminosity = c("light"))
pink <- c("#f41f8d", "#ffb7de", "#ed6add", "#fd8cff", "#cc0076", "#ffbff9", "#f271d6", "#f279b5", "#f47cb8","#eec4fc", "#ff7cc6", "#ed6ab2")


#purple = mesenchymal and endothelial 
purple <- randomColor(3, hue = c("purple"))
purple <- c("#d3adea", "#7214ad", "#ab26ff")

#red - ribosomal and doublet cluster 
red <- randomColor(2, hue = c("red"))
red <- c("#e8480d", "#dd2a02")

#group cell types 
ileum_low$cell_typev1 <- factor(ileum_low$cell_typev1, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell 1", "Plasma Cell 2", "Plasmablast", "CD4 T Cell", "CD8 T Cell", "NK/T Cell", "Ambiguous T Cell", "ILC3", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte 1", "Enterocyte 2", "Enterocyte 3", "Stem/Paneth Cell", "Goblet Cell 1", "Goblet Cell 2", "TA Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell", "Doublet Cluster"))

ileum_low$cell_typev2 <- factor(ileum_low$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "CD8 T Cell", "NK/T Cell", "Ambiguous T Cell", "ILC3", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte", "Stem/Paneth Cell", "Goblet Cell", "TA Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell", "Doublet Cluster"))

DimPlot(ileum_low, reduction = "umap", group.by = "cell_typev1", cols = c(green[1], green[2], green[3], green[4], green[5], blue[1], blue[2], blue[3], blue[8], blue[5], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[4], pink[5], pink[6], pink[7], pink[10], purple[1], purple[2], red[1]))

DimPlot(ileum_low, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[3], green[5], blue[1], blue[2], blue[3], blue[8], blue[5], orange[1], orange[2], orange[3], orange[4], pink[1], pink[4], pink[5], pink[7], pink[10], purple[1], purple[2], red[1]))


#proportion of cells 
ileum_low$orig.ident <- factor(ileum_low$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum_low) <- ileum_low$cell_typev1
pt <- table(Idents(ileum_low), ileum_low$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell 1", "Plasma Cell 2", "Plasmablast", "CD4 T Cell", "CD8 T Cell", "NK/T Cell", "Ambiguous T Cell", "ILC3", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte 1", "Enterocyte 2", "Enterocyte 3", "Stem/Paneth Cell", "Goblet Cell 1", "Goblet Cell 2", "TA Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell", "Doublet Cluster"))


plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[4], green[5], blue[1], blue[2], blue[3], blue[8], blue[5], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[4], pink[5], pink[6], pink[7], pink[10], purple[1], purple[2], red[1])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#cell type v2
ileum_low$orig.ident <- factor(ileum_low$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum_low) <- ileum_low$cell_typev2
pt <- table(Idents(ileum_low), ileum_low$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "CD8 T Cell", "NK/T Cell", "Ambiguous T Cell", "ILC3", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte", "Stem/Paneth Cell", "Goblet Cell", "TA Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell", "Doublet Cluster"))


plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[5], blue[1], blue[2], blue[3], blue[8], blue[5], orange[1], orange[2], orange[3], orange[4], pink[1], pink[4], pink[5], pink[7], pink[10], purple[1], purple[2], red[1])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))


#--------cluster 2---------#
#parameters: res: 1, PC: 15, HVG: 2000

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_1_15_2000_rpca_10_29_24.rds')

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32'))

library(randomcoloR)
no_of_colors <- 33

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#--------cluster 3---------#
#parameters: res: 0.5, PC: 15, HVG: 5000

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_50_15_5000_rpca_10_29_24.rds')

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24'))

library(randomcoloR)
no_of_colors <- 25

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))


#-----use ileum cluster 1 and remove doublet cluster--------#

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_50_15_2000_rpca_10_29_24.rds')

DefaultAssay(ileum) <- "RNA"

#set the Ident
Idents(ileum) <- ileum$seurat_clusters

#confirm that cluster 22 is doublet 
cluster22 <- FindMarkers(ileum, ident.1 = "22")
#yes - doublet

#subset cluster 22
ileum_sub <- subset(ileum, idents = "22", invert = T)

#save the seurat obejct 
saveRDS(ileum_sub, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_50_15_2000_rpca_sub_10_30_24.rds")

#-----------re-assess clustering results of ileum cluster 1---------------#

#load seurat object 
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_50_15_2000_rpca_sub_cluster_10_31_24.rds')

#parameters: res: 0.5, PC: 15, HVG: 2000

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21'))

library(randomcoloR)
no_of_colors <- 25

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#---------Preliminary Annotation of Ileum Cluster 1 Subset--------------#

#ileum cluster 1 subset 
DefaultAssay(ileum) <- "RNA"
Idents(ileum) <- ileum$seurat_clusters


#-----T cells-----#

#naive/central memory CD4 T cell
DotPlot(ileum, features = c("CD3D", 'CD4', 'CCR7', 'SELL'))

#tissue resident memory CD4 T cell
DotPlot(ileum, features = c("CD3D", 'CD4', 'ITGAE', 'ITGA1', 'SPRY1'))

#tissue resident memory T helper 17 
DotPlot(ileum, features = c("CD3D", 'CD4', 'ITGAE', 'ITGA1', 'SPRY1', 'IL17A', 'CCR6', 'CCL20'))

#naive t follicular helper 
DotPlot(ileum, features = c("CD3D", 'CD4', 'SELL', 'CXCR5', 'BCL6', 'CD40LG'))

#t follicular helper cell 
DotPlot(ileum, features = c("CD3D", 'CD4', 'CXCR5', 'BCL6', 'PDCD1', 'CD40LG'))

#regulatory t cell 
DotPlot(ileum, features = c("CD3D", 'CD4', 'FOXP3', 'CTLA4', 'IL2RA'))

#regulatory t cell IL10+
DotPlot(ileum, features = c("CD3D", 'CD4', 'FOXP3', 'CTLA4', 'IL2RA', 'IL10'))

#naive/central memory CD8 t cell 
DotPlot(ileum, features = c("CD3D", 'CD8A', 'CD8B', 'CCR7', 'SELL'))

#tissue resident memory CD8 t cell 
DotPlot(ileum, features = c("CD3D", 'CD8A', 'CD8B', 'ITGAE', 'ITGA1', 'SPRY1'))

#tissue resident memory/effector memory Cd8 t cell 
DotPlot(ileum, features = c("CD3D", 'CD8A', 'CD8B', 'GZMK', 'CRTAM', 'EOMES'))

#naive gamma delta t cell 
DotPlot(ileum, features = c("CD3D", 'SELL', 'KLRC2', 'TRDC', 'TRGC1', 'KIR2DL4', 'CCR7'))

#gamma delta t cell 
DotPlot(ileum, features = c("CD3D", 'KLRC2', 'TRDC', 'TRGC1', 'KIR2DL4'))

#MAIT (mucosal associated invariant t cell)
DotPlot(ileum, features = c("CD3D", 'SLC4A10', 'TRAV1-2'))

#CD16+ NK cell 
DotPlot(ileum, features = c('KLRF1', 'NKG7', 'FCGR3A', 'GNLY', 'GZMB'))

#CD56 bright natural killer cell 
DotPlot(ileum, features = c('KLRF1', 'NKG7', 'XCL1', 'IL2RB', 'NCR1', 'FCER1G', 'NCAM1'))

#ILC3 innate lymphoid cell type 3 (cd3d negative)
DotPlot(ileum, features = c('IL7R', 'RORC', 'KIT', 'LST1', 'PCDH9', 'IL1R1', 'IL23R'))

#cycling t or NK cell
DotPlot(ileum, features = c('CD3D', 'MKI67', 'TOP2A'))

#--------myeloid------------#

#type 1 conventional DC 
DotPlot(ileum, features = c('HLA-DRA', 'HLA-DPA1', 'CLEC9A', 'XCR1', 'BATF3', 'CADM1', 'RAB7B'))

#type 2 conventional DC
DotPlot(ileum, features = c('HLA-DRA', 'HLA-DPA1', 'CLEC10A', 'FCER1A', 'CD1C'))

#migratory DC
DotPlot(ileum, features = c('HLA-DRA', 'HLA-DPA1', 'CCR7', 'LAMP3'))

#plasmacytoid DC
DotPlot(ileum, features = c('HLA-DRA', 'HLA-DPA1', 'IRF7', 'CLEC4C', 'JCHAIN', 'LILRA4', 'GZMB'))

#langerhans DC
DotPlot(ileum, features = c('HLA-DRA', 'HLA-DPA1', 'ITGAX', 'IL22RA2', 'CD207', 'RUNX3'))

#monocyte
DotPlot(ileum, features = c('FCN1', 'S100A8', 'S100A9', 'IL1B', 'EREG', 'NAMPT', 'PLAUR', 'VCAN', 'FPR1', 'CD3D00E'))

#LYVE1 macrophage
DotPlot(ileum, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'LYVE1', 'RNASE1', 'FOLR2'))

#MMP9 macrophage
DotPlot(ileum, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'MMP9', 'PLA2G2D', 'ADAMDEC1'))

#TREM2 macrophage 
DotPlot(ileum, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'TREM2', 'ACP5', 'CTSD', 'CSTB'))

#CD5L macrophage
DotPlot(ileum, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'CD5L', 'VCAM1', 'CXCL12', 'PDK4', 'RBP7'))

#macrophage
DotPlot(ileum, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'CD209'))

#mast cell 
DotPlot(ileum, features = c('CD69', 'KIT', 'TPSB2', 'TPSAB1'))

#eosinophil/basophil
DotPlot(ileum, features = c('GATA2', 'CNRIP1', 'PRG2', 'GIHCG', 'CLC'))

#erythrocytes
DotPlot(ileum, features = c('GATA1', 'HBZ', 'HBE1', 'HBG1'))

#monocyte/neutrophil progenitor
DotPlot(ileum, features = c('FCN1', 'S100A8', 'S100A9', 'MPO', 'RETN', 'RNASE2', 'PCLAF'))

#megakaryocyte/platelet
DotPlot(ileum, features = c('GATA1', 'TAL1', 'MMRN1', 'CMTM5', 'MPIG6B', 'ITGA2B', 'PF4'))

#----------b cells-------------#

#pre b cell 
DotPlot(ileum, features = c('CD19', 'HLA-DRA', 'CD79B', 'SPIB', 'TCL1A', 'CD3D7'))

#pro b cell 
DotPlot(ileum, features = c('CD19', 'HLA-DRA', 'CD79B', 'IGLL1', 'RAG1', 'DNTT', 'VPREB3'))

#naive b cell 
DotPlot(ileum, features = c('CD19', 'HLA-DRA', 'CD79A', 'MS4A1', 'SELL', 'TCL1A', 'IGHD'))

#memory b cell 
DotPlot(ileum, features = c('CD19', 'HLA-DRA', 'CD79A', 'MS4A1', 'CD27', 'TNFSF13B'))

#germinal center b cell 1
DotPlot(ileum, features = c('CD19', 'HLA-DRA', 'CD79A', 'MS4A1', 'MKI67', 'AICDA', 'BCL6', 'SUGCT'))

#plasmablast
DotPlot(ileum, features = c('MZB1', 'JCHAIN', 'XBP1'))

#IgM plasma cell 
DotPlot(ileum, features = c('MZB1', 'JCHAIN', 'IGHM'))

#IgA plasma cell 
DotPlot(ileum, features = c('MZB1', 'JCHAIN', 'IGHA1', 'IGHA2'))

#IgG plasma cell 
DotPlot(ileum, features = c('MZB1', 'JCHAIN', 'IGHG3', 'IGHG1', 'IGHG2', 'IGHG4'))

#---------endothelial cells-----------#

#capillary endothelial cell 
DotPlot(ileum, features = c('PECAM1', 'CD3D6', 'RGCC', 'COL4A1', 'COL4A2', 'IL32', 'MCAM', 'MYO1B'))

#arterial endothelial cell 
DotPlot(ileum, features = c('PECAM1', 'CD3D6', 'GJA4', 'HEY1', 'CXCL12', 'SEMA3G', 'IGFBP3', 'FBLN2', 'FBLN5', 'ELN', 'BTNL9', 'ALPL'))

#venous endothelial cell
DotPlot(ileum, features = c('PECAM1', 'CD3D6', 'ACKR1', 'CCL14', 'SELE', 'TNFRSF6B'))

#lymphatic endothelial cell 
DotPlot(ileum, features = c('PECAM1', 'CD3D6', 'CCL21', 'TFF3', 'PROX2', 'NTS'))

#cycling endothelial cell 
DotPlot(ileum, features = c('PECAM1', 'CD3D6', 'MKI67', 'TOP2A'))

#-------------mesenchymal cells-------------#

#vascular smooth muscle cell 
DotPlot(ileum, features = c('VIM', 'DCN', 'TAGLN', 'ACTA2', 'TPM2', 'MYH11', 'RERGL', 'MUSTN1', 'LBH', 'NET1', 'MAP3K20'))

#pericyte
DotPlot(ileum, features = c('VIM', 'DCN', 'TAGLN', 'ACTA2', 'TPM2', 'MYH11', 'COX4I2', 'HIGD1B', 'RGS5', 'NDUFA4L2'))

#immune recruiting pericyte
DotPlot(ileum, features = c('VIM', 'DCN', 'TAGLN', 'ACTA2', 'TPM2', 'MYH11', 'GPC3', 'COL14A1', 'ECRG4', 'ID4', 'FHL2', 'CXCL12'))

#myofibroblast
DotPlot(ileum, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'ACTG2', 'HHIP', 'SOSTDC1', 'NPNT'))

#follicular DC
DotPlot(ileum, features = c('VIM', 'FDCSP', 'SRGN', 'CR2', 'CLU', 'CSTA'))

#reticular fibroblast 
DotPlot(ileum, features = c('VIM', 'CCL21', 'CCL19', 'TNFSF13B', 'TDO2'))

#oral mucosa fibroblast
DotPlot(ileum, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'CTHRC1', 'COL12A1', 'COL1A1', 'CTSK', 'COL5A2'))

#oesophagus fibroblast
DotPlot(ileum, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'APOD', 'PLPP1', 'MFAP4', 'IFITM1', 'RASD1'))

#crypt fibroblast
DotPlot(ileum, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'PI16', 'RSPO3', 'SFRP1', 'TM2A'))

#villus fibroblast
DotPlot(ileum, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'F3', 'PLAT', 'HSD17B2', 'SOX6'))

#lamina propria fibroblast
DotPlot(ileum, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'ADAMDEC1', 'ADAM28', 'CCL11', 'CCL8', 'CCL13', 'CFD'))

#rectum fibroblast
DotPlot(ileum, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'KCNN3', 'THBS4', 'FNDC1', 'PPFIBP1'))

#mesothelium 
DotPlot(ileum, features = c('VIM', 'DCN', 'UPK3B', 'MSLN', 'SLPI', 'PLAT', 'KRT19'))

#-----------epithelium-------------#

#---------small intestine-----------#

#entereocyte
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'FABP1', 'APOA4', 'PRAP1', 'PCK1', 'RBP2', 'SI'))

#BEST4 enterocyte
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'BEST4', 'CA7', 'OTOP2', 'CFTR'))

#stem cell 
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'LGR5', 'RGMB', 'ASCL2', 'OLFM4'))

#enteroendocrine cell 
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'CHGA', 'PCSK1N', 'SCT', 'SCGN', 'NEUROD1'))

#microfold cell 
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'IL2RG', 'ICAM2', 'CCL20', 'CCL23'))

#tuft cell 
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'SH2D6', 'LRMP', 'MATK', 'FYB1', 'HPGDS', 'POU2F3', 'TRPM5'))

#goblet cell progenitor
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'GAU1', 'MUC2', 'TFF3', 'FCGBP', 'ZG16', 'OLFM4'))

#goblet cell 
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'MUC2', 'TFF3', 'FCGBP', 'ZG16'))

#proliferating goblet cell 
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'MUC2', 'TFF3', 'FCGBP', 'ZG16', 'MKI67', 'TOP2A'))

#transit amplifying cell 
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'MKI67', 'TOP2A', 'PCLAF', 'PCNA'))

#paneth cell 
DotPlot(ileum, features = c('CDH1', 'KRT19', 'EPCAM', 'DEFA6', 'DEFA5', 'REG3A', 'PLA2G2A'))

#-------preliminary cluster annotations-------#

#11/4/2024

ileum@meta.data$cell_typev1 <- ileum@meta.data$seurat_clusters
ileum$cell_typev1 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21'),
  to = c("Naive B Cell", 'Stem/Paneth Cell', 'Enterocyte 1', 'CD4 T Cell', 'NK/CD8 T Cell', 'Enterocyte 2', 'Macrophage', "Plasma Cell 1", "Plasma Cell 2", 'Monocyte', 'Enterocyte 3', 'Goblet Cell 1', 'Goblet Cell 2', 'DNT Cell', 'Mesenchymal Cell', "Pro-Inflammatory Monocyte", "Cycling Cell", "Memory B Cell", "Endothelial Cell", 'Mast Cell', "Tuft Cell", "Cycling Plasma Cell")
)


ileum@meta.data$cell_typev2 <- ileum@meta.data$seurat_clusters
ileum$cell_typev2 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21'),
  to = c("Naive B Cell", 'Stem/Paneth Cell', 'Enterocyte', 'CD4 T Cell', 'NK/CD8 T Cell', 'Enterocyte', 'Macrophage', "Plasma Cell", "Plasma Cell", 'Monocyte', 'Enterocyte', 'Goblet Cell', 'Goblet Cell', 'DNT Cell', 'Mesenchymal Cell', "Pro-Inflammatory Monocyte", "Cycling Cell", "Memory B Cell", "Endothelial Cell", 'Mast Cell', "Tuft Cell", "Cycling Plasma Cell")
)

#colors for UMAP and stacked bar plot 

#blue = t cells 
blue <- randomColor(3, hue = c("blue"))
blue <- c("#34bfe5", "#4a93ce", "#b0daf2") 


#green = b cells
green <- randomColor(5, hue = c("green"))
green <- c("#05fc95", "#9de587", "#1da54f", "#61e291", "#2bd147")

#orange = myeloid 
orange <- randomColor(4, hue = c("orange"))
orange <- c("#fcab85", "#d3b050", "#f75f13", "#cc5414")

#pink = epithelial
pink <- randomColor(8, hue = c("pink"), luminosity = c("light"))
pink <- c("#f41f8d", "#ffb7de", "#ed6add", "#fd8cff", "#cc0076", "#ffbff9", "#f271d6", "#f279b5")

#purple = mesenchymal and endothelial 
purple <- randomColor(2, hue = c("purple"))
purple <- c("#cd7af9", "#996cc9")

#group cell types 
ileum$cell_typev1 <- factor(ileum$cell_typev1, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell 1", "Plasma Cell 2", "Cycling Plasma Cell", "CD4 T Cell", "NK/CD8 T Cell", "DNT Cell", "Macrophage", "Monocyte", "Pro-Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Enterocyte 1", "Enterocyte 2", "Enterocyte 3", "Goblet Cell 1", "Goblet Cell 2", "Cycling Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

ileum$cell_typev2 <- factor(ileum$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Cycling Plasma Cell", "CD4 T Cell", "NK/CD8 T Cell", "DNT Cell", "Macrophage", "Monocyte", "Pro-Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Enterocyte", "Goblet Cell", "Cycling Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

DimPlot(ileum, reduction = "umap", group.by = "cell_typev1", cols = c(green[1], green[2], green[3], green[4], green[5], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[4], pink[5], pink[6], pink[7], pink[8], purple[1], purple[2]))


#proportion of cells 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$cell_typev1
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell 1", "Plasma Cell 2", "Cycling Plasma Cell", "CD4 T Cell", "NK/CD8 T Cell", "DNT Cell", "Macrophage", "Monocyte", "Pro-Inflammatory Monocyte", "Mast Cell", "Paneth/Stem Cell", "Enterocyte 1", "Enterocyte 2", "Enterocyte 3", "Goblet Cell 1", "Goblet Cell 2", "Cycling Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))


plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[4], green[5], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[4], pink[5], pink[6], pink[7], pink[8], purple[1], purple[2])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))



DimPlot(ileum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[3], green[5], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[3], pink[6], pink[7], pink[8], purple[1], purple[2]))


#proportion of cells 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$cell_typev2
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Cycling Plasma Cell", "CD4 T Cell", "NK/CD8 T Cell", "DNT Cell", "Macrophage", "Monocyte", "Pro-Inflammatory Monocyte", "Mast Cell", "Paneth/Stem Cell", "Enterocyte", "Goblet Cell", "Cycling Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))


plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[5], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[3], pink[6], pink[7], pink[8], purple[1], purple[2])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#do FM comparisons between clusters 

DefaultAssay(ileum) <- "RNA"
Idents(ileum) <- ileum$seurat_clusters

#compare cluster 2 vs. cluster 5 vs. cluster 10
c2_vs_5_10 <- FindMarkers(ileum, ident.1 = "2", ident.2 = c("5", "10"))

write.csv(c2_vs_5_10, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/fm_ileum_cluster2_vs_cluster5_10_param1_11_5_25.csv")

c5_vs_2_10 <- FindMarkers(ileum, ident.1 = "5", ident.2 = c("2", "10"))

write.csv(c5_vs_2_10, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/fm_ileum_cluster5_vs_cluster2_10_param1_11_5_25.csv")

c10_vs_2_5 <- FindMarkers(ileum, ident.1 = "10", ident.2 = c("2", "5"))

write.csv(c10_vs_2_5, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/fm_ileum_cluster10_vs_cluster2_5_param1_11_5_25.csv")

#compare cluster 7 vs. 8 
c7_vs_8 <- FindMarkers(ileum, ident.1 = "7", ident.2 = c("8"))

write.csv(c7_vs_8, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/fm_ileum_cluster7_vs_cluster8_param1_11_5_25.csv")

#compare cluster 11 vs. 12
c11_vs_12 <- FindMarkers(ileum, ident.1 = "11", ident.2 = c("12"))

write.csv(c11_vs_12, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/fm_ileum_cluster11_vs_cluster12_param1_11_5_25.csv")

#compare cluster 9 vs. 15
c9_vs_15 <- FindMarkers(ileum, ident.1 = "9", ident.2 = c("15"))

write.csv(c9_vs_15, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/fm_ileum_cluster9_vs_cluster15_param1_11_5_25.csv")



#include broader cell labels 
ileum@meta.data$ctypes <- ileum@meta.data$cell_typev2
ileum$ctypes <- plyr::mapvalues(
  x = ileum$cell_typev2,
  from = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Cycling Plasma Cell", "CD4 T Cell", "NK/CD8 T Cell", "DNT Cell", "Macrophage", "Monocyte", "Pro-Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Enterocyte", "Goblet Cell", "Cycling Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"),
  to = c("B Cell", "B Cell", "Plasma Cell", "Plasma Cell", "NK/T Cell", "NK/T Cell", "Ambig. T Cell", "Myeloid Cell", "Myeloid Cell", "Myeloid Cell", "Myeloid Cell", "Stem/Paneth Cell", "Enterocyte", "Goblet Cell", 'TA Cell', "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell")
)

#save the seurat object 
saveRDS(ileum, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_11_19_24.rds")

#--------add group metadata information from scITD analysis---------#
#11/20/2024

#subset donor 13, 29, 31, and 33 - was not kept in scITD analysis 
Idents(ileum) <- ileum$donor
ileum_sub <- subset(ileum, idents = c("donor13", "donor29", "donor31", "donor33"), invert = T)

#group 1 = positive score
#group 2 = negative score

ileum_sub@meta.data$group <- ileum_sub@meta.data$orig.ident
ileum_sub$group <- plyr::mapvalues(
  x = ileum_sub$orig.ident,
  from = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam147', 'helm_sam171'),
  to = c('group1', 'group1', 'group1', 'group1', 'group1', 'group2', 'group2', 'group1', 'group1', 'group2', 'group1', 'group2', 'group1', 'group2', 'group1', 'group1', 'group2', 'group1', 'group2', 'group2', 'group1', 'group2', 'group2')
)

#save the ileum seurat object with group information 
saveRDS(ileum_sub, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_scITD_groupInfo_11_20_24.rds")

#prep for dreamlet - have to convert SO to single cell experiment 

#code below 
DefaultAssay(ileum) <- "RNA"
ileum.sce <- as.SingleCellExperiment(ileum)
saveRDS(ileum.sce, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_scITD_group_sce_11_20_24.rds")

#------cell type proportion analysis with speckle---------#

library(speckle)
library(limma)
library(ggplot2)

# Run propeller testing for cell type proportion differences between the two 
# groups
macro_if <- propeller(clusters = ileum_all$cell_typev1, sample = ileum_all$orig.ident, group = ileum_all$macro_IF)

#not significant

macro_if <- propeller(clusters = ileum_all$cell_typev2, sample = ileum_all$orig.ident, group = ileum_all$macro_IF)

#not signficant

micro_if <- propeller(clusters = ileum_all$cell_typev1, sample = ileum_all$orig.ident, group = ileum_all$micro_IF)

#monocyte and cycling plasma cell are marginally significant (fdr p value = 0.051)

micro_if <- propeller(clusters = ileum_all$cell_typev2, sample = ileum_all$orig.ident, group = ileum_all$micro_IF)

#monocyte and cycling plasma cell associated with micro if status (fdr p value = 0.0359 for both)

batch <- propeller(clusters = ileum_all$cell_typev1, sample = ileum_all$orig.ident, group = ileum_all$batch)

#not significant 

batch <- propeller(clusters = ileum_all$cell_typev2, sample = ileum_all$orig.ident, group = ileum_all$batch)

#not significant 

sex <- propeller(clusters = ileum_all$cell_typev1, sample = ileum_all$orig.ident, group = ileum_all$sex)

#not significant

sex <- propeller(clusters = ileum_all$cell_typev2, sample = ileum_all$orig.ident, group = ileum_all$sex)

#not signicant 

ancestry <- propeller(clusters = ileum_all$cell_typev1, sample = ileum_all$orig.ident, group = ileum_all$ancestry)

#not significant 

ancestry <- propeller(clusters = ileum_all$cell_typev2, sample = ileum_all$orig.ident, group = ileum_all$ancestry)

#not significant 

#--------------High Resolution Stable Cells Top Clustering Parameters---------------#

#11/11/2024 

#-------cluster 1---------#

#parameters: res: 0.75, PC: 15, HVG: 2000

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_13_ileum_75_15_2000_rpca_11_11_24.rds')

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36'))

library(randomcoloR)
no_of_colors <- 37

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#AD6F3B", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#prelim annotation 
DefaultAssay(ileum) <- "RNA"
Idents(ileum) <- ileum$seurat_clusters

FeaturePlot(ileum, features = c("EPCAM", "PTPRC"))

#cluster 35?
#cluster 36?
cluster36 <- FindMarkers(ileum, ident.1 = "36")


#-------cluster 2---------#

#parameters: res: 0.50, PC: 15, HVG: 2000

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_13_ileum_50_15_2000_rpca_11_11_24.rds')

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30'))

library(randomcoloR)
no_of_colors <- 31

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#AD6F3B", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#preliminary annotation 
DefaultAssay(ileum) <- "RNA"
Idents(ileum) <- ileum$seurat_clusters

FeaturePlot(ileum, features = c("EPCAM", "PTPRC"))

#-------cluster 3---------#

#parameters: res: 1, PC: 15, HVG: 2000

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_13_ileum_1_15_2000_rpca_11_11_24.rds')

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37'))

library(randomcoloR)
no_of_colors <- 38

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#AD6F3B", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))


#preliminary annotation 
DefaultAssay(ileum) <- "RNA"
Idents(ileum) <- ileum$seurat_clusters

FeaturePlot(ileum, features = c("EPCAM", "PTPRC"))

#----------cell annotation of param 1 high res clustering results--------------#

#11/17/2024

ileum_high@meta.data$cell_typev1 <- ileum_high@meta.data$seurat_clusters
ileum_high$cell_typev1 <- plyr::mapvalues(
  x = ileum_high$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36'),
  to = c("Naive B Cell 1", "Memory B Cell", "Enterocyte 1", "Stem/Paneth Cell 1", "CD4 T Cell", "Enterocyte 2", "Stem/Paneth Cell 2", 'CD8 T Cell', "Plasma Cell 1", "Enterocyte 3", "Monocyte", "Plasma Cell 2", "Enterocyte 4", "Goblet Cell 1", "Stem Cell", "Regulatory T Cell", "Ribosomal Cluster", "Goblet Cell 2", "Inflammatory Macrophage", "Macrophage", "Ambiguous T Cell", "Inflammatory Monocyte", "TA Cell", "NK/T Cell", "Mesenchymal Cell 1", "Germinal Center B Cell", "Endothelial Cell", "Mesenchymal Cell 2", "Mast Cell", "Tuft Cell", "Doublet Cluster", "Plasma Cell 3", "Naive B Cell 2", "ILC3", "Enteroendocrine", "Plasmablast", "Cycling NK/T Cell")
)

ileum_high@meta.data$cell_typev2 <- ileum_high@meta.data$seurat_clusters
ileum_high$cell_typev2 <- plyr::mapvalues(
  x = ileum_high$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36'),
  to = c("Naive B Cell", "Memory B Cell", "Enterocyte", "Stem/Paneth Cell", "CD4 T Cell", "Enterocyte", "Stem/Paneth Cell", 'CD8 T Cell', "Plasma Cell", "Enterocyte", "Monocyte", "Plasma Cell", "Enterocyte", "Goblet Cell", "Stem Cell", "Regulatory T Cell", "Ribosomal Cluster", "Goblet Cell", "Inflammatory Macrophage", "Macrophage", "Ambiguous T Cell", "Inflammatory Monocyte", "TA Cell", "NK/T Cell", "Mesenchymal Cell", "Germinal Center B Cell", "Endothelial Cell", "Mesenchymal Cell", "Mast Cell", "Tuft Cell", "Doublet Cluster", "Plasma Cell", "Naive B Cell", "ILC3", "Enteroendocrine", "Plasmablast", "Cycling NK/T Cell")
)

ileum_high$cell_typev1 <- factor(ileum_high$cell_typev1, levels = c("Naive B Cell 1", "Naive B Cell 2", "Memory B Cell", "Germinal Center B Cell", "Plasma Cell 1", "Plasma Cell 2", "Plasma Cell 3", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/T Cell", "Cycling NK/T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Inflammatory Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte 1", "Enterocyte 2", "Enterocyte 3", "Enterocyte 4", "Stem/Paneth Cell 1", "Stem/Paneth Cell 2", "Stem Cell", "TA Cell", "Goblet Cell 1", "Goblet Cell 2", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell 1", "Mesenchymal Cell 2", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))

ileum_high$cell_typev2 <- factor(ileum_high$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Germinal Center B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/T Cell", "Cycling NK/T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Inflammatory Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte", "Stem/Paneth Cell", "Stem Cell", "TA Cell", "Goblet Cell", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))

#blue = t cells 
blue <- randomColor(8, hue = c("blue"))
blue <- c("#34bfe5", "#4a93ce", "#b0daf2", "#6db5ed", "#1c5682","#077b9e", "#70c8e0", "#1f86f4") 

#green = b cells
green <- randomColor(8, hue = c("green"))
green <- c("#05fc95", "#9de587", "#1da54f", "#61e291", "#2bd147","#80cc1e", "#c5ef99", "#048c45")


#orange = myeloid 
orange <- randomColor(5, hue = c("orange"))
orange <- c("#fcab85", "#d3b050", "#f75f13", "#cc5414", "#f98a52")

#pink = epithelial
pink <- randomColor(12, hue = c("pink"), luminosity = c("light"))
pink <- c("#f41f8d", "#ffb7de", "#ed6add", "#fd8cff", "#cc0076", "#ffbff9", "#f271d6", "#f279b5", "#f47cb8","#eec4fc", "#ff7cc6", "#ed6ab2")


#purple = mesenchymal and endothelial 
purple <- randomColor(3, hue = c("purple"))
purple <- c("#d3adea", "#7214ad", "#ab26ff")

#red - ribosomal and doublet cluster 
red <- randomColor(2, hue = c("red"))

red <- c("#e8480d", "#dd2a02")

DimPlot(ileum_high, reduction = "umap", group.by = "cell_typev1", cols = c(green[1], green[2], green[3], green[4], green[5], green[6], green[7], green[8], blue[1], blue[2], blue[3], blue[4], blue[5], blue[6], blue[8], orange[1], orange[2], orange[3], orange[4], orange[5], pink[1], pink[2], pink[3], pink[4], pink[5], pink[6], pink[7], pink[8], pink[9], pink[10], pink[11], pink[12], purple[1], purple[2], purple[3], red[2], red[1]))

DimPlot(ileum_high, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[3], green[4], green[5], green[8], blue[1], blue[2], blue[3], blue[4], blue[5], blue[6], blue[8], orange[1], orange[2], orange[3], orange[4], orange[5], pink[1], pink[5], pink[7], pink[8], pink[10], pink[11], pink[12], purple[1], purple[3], red[2], red[1]))

#stacked bar plots 

#proportion of cells - cell type v1
ileum_high$orig.ident <- factor(ileum_high$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum_high) <- ileum_high$cell_typev1
pt <- table(Idents(ileum_high), ileum_high$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell 1", "Naive B Cell 2", "Memory B Cell", "Germinal Center B Cell", "Plasma Cell 1", "Plasma Cell 2", "Plasma Cell 3", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/T Cell", "Cycling NK/T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Inflammatory Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte 1", "Enterocyte 2", "Enterocyte 3", "Enterocyte 4", "Stem/Paneth Cell 1", "Stem/Paneth Cell 2", "Stem Cell", "TA Cell", "Goblet Cell 1", "Goblet Cell 2", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell 1", "Mesenchymal Cell 2", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))


plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[4], green[5], green[6], green[7], green[8], blue[1], blue[2], blue[3], blue[4], blue[5], blue[6], blue[8], orange[1], orange[2], orange[3], orange[4], orange[5], pink[1], pink[2], pink[3], pink[4], pink[5], pink[6], pink[7], pink[8], pink[9], pink[10], pink[11], pink[12], purple[1], purple[2], purple[3], red[2], red[1])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells v2
ileum_high$orig.ident <- factor(ileum_high$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum_high) <- ileum_high$cell_typev2
pt <- table(Idents(ileum_high), ileum_high$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Germinal Center B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/T Cell", "Cycling NK/T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Inflammatory Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte", "Stem/Paneth Cell", "Stem Cell", "TA Cell", "Goblet Cell", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))


plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[3], green[4], green[5], green[8], blue[1], blue[2], blue[3], blue[4], blue[5], blue[6], blue[8], orange[1], orange[2], orange[3], orange[4], orange[5], pink[1], pink[5], pink[7], pink[8], pink[10], pink[11], pink[12], purple[1], purple[3], red[2], red[1])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#--------------Reference Top Clustering Parameters---------------#

#11/11/2024 

#-------cluster 1---------#

#parameters: res: 0.75, PC: 30, HVG: 2000

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/seurat_object/helm_batch1_13_ileum_75_30_2000_ref_rpca_11_11_24.rds')

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35'))

library(randomcoloR)
no_of_colors <- 36

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#AD6F3B", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#preliminary cell type assignment 
DefaultAssay(ileum) <- "RNA"
Idents(ileum) <- ileum$seurat_clusters

FeaturePlot(ileum, features = c("EPCAM", "PTPRC"))

#cluster 15? cluster 31?
cluster31 <- FindMarkers(ileum, ident.1 = "31")
cluster15 <- FindMarkers(ileum, ident.1 = "15")


#-------cluster 2---------#

#parameters: res: 0.75, PC: 15, HVG: 500

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/seurat_object/helm_batch1_13_ileum_75_15_500_ref_rpca_11_11_24.rds')

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30'))

library(randomcoloR)
no_of_colors <- 31

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#AD6F3B", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#preliminary annotation 
DefaultAssay(ileum) <- "RNA"
Idents(ileum) <- ileum$seurat_clusters

FeaturePlot(ileum, features = c("EPCAM", "PTPRC"))

#cluster 12?
cluster12 <- FindMarkers(ileum, ident.1 = "12")

#cluster 18?
cluster18 <- FindMarkers(ileum, ident.1 = "18")


#-------cluster 3---------#

#parameters: res: 1, PC: 30, HVG: 5000

#load in seurat object
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/seurat_object/helm_batch1_13_ileum_1_30_5000_ref_rpca_11_11_24.rds')

ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
ileum$batch <- factor(ileum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize QC parameters of clusters 
VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0)
VlnPlot(ileum, features = "nCount_RNA", pt.size = 0)
VlnPlot(ileum, features = "percent.mt", pt.size = 0)

#visualize PCA results
ElbowPlot(ileum)
DimPlot(ileum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(ileum, reduction = 'pca', group.by = "orig.ident")
DimPlot(ileum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(ileum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(ileum, reduction = 'umap', group.by = "orig.ident")
DimPlot(ileum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
ileum$orig.ident <- factor(ileum$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum) <- ileum$seurat_clusters
pt <- table(Idents(ileum), ileum$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37'))

library(randomcoloR)
no_of_colors <- 38

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = palette
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))

#proportion of cells per batch
Idents(ileum) <- ileum$batch
pt2 <- table(Idents(ileum), ileum$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c("#DA5724", "#CE50CA", "#74D944", "#3F4921", "#89C5DA", "#C0717C", "#AD6F3B", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#00CCCC", "#0033FF", "#FF3300","#599861", "#CC33CC", "#FFFF00", "#990033", "#3399FF"
       )) +
       theme_bw(base_size=10) +
       geom_col(position = "fill", width = 0.5) +
       # xlab("Sample") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       # ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=10)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
       theme(axis.text.y = element_text(size=10)))

#preliminary cell type assignment check 
DefaultAssay(ileum) <- "RNA"
Idents(ileum) <- ileum$seurat_clusters

#unknowns - cluster 14 and cluster 33

cluster14 <- FindMarkers(ileum, ident.1 = "14")

FeaturePlot(ileum, features = c("EPCAM", "PTPRC"))

#----------annotation for robustness check----------#

#used FM results and checked with marker genes 
ileum_ref@meta.data$cell_typev1 <- ileum_ref@meta.data$seurat_clusters
ileum_ref$cell_typev1 <- plyr::mapvalues(
  x = ileum_ref$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35'),
  to = c("Naive B Cell 1", "Enterocyte 1", 'Memory B Cell', 'Enterocyte 2', 'Stem/Paneth Cell 1', "Stem/Paneth Cell 2", 'Naive/CM CD4 T Cell', 'Macrophage', "Plasma Cell 1", 'Monocyte', 'Stem Cell', 'Plasma Cell 2', 'Goblet Cell 1', 'Enterocyte 3', 'Goblet Cell 2', 'Ribosomal Cluster', 'Regulatory T Cell', 'NK/CD8 T Cell', 'Tissue Resident Memory CD8 T Cell', 'Ambiguous T Cell', 'Inflammatory Monocyte', 'TA Cell', 'Cycling B Cell', 'Tissue Resident Memory CD4 T Cell', 'Mesenchymal Cell 1', 'Naive/CM CD8 T Cell', 'Mesenchymal Cell 2', 'Endothelial Cell', 'Tuft Cell', 'Mast Cell', 'Enterocyte 4', 'Naive B Cell 2', 'Doublet Cluster', 'Enteroendocrine', 'ILC3', 'Plasmablast')
)

ileum_ref@meta.data$cell_typev2 <- ileum_ref@meta.data$seurat_clusters
ileum_ref$cell_typev2 <- plyr::mapvalues(
  x = ileum_ref$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35'),
  to = c("Naive B Cell", "Enterocyte", 'Memory B Cell', 'Enterocyte', 'Stem/Paneth Cell', "Stem/Paneth Cell", 'CD4 T Cell', 'Macrophage', "Plasma Cell", 'Monocyte', 'Stem Cell', 'Plasma Cell', 'Goblet Cell', 'Enterocyte', 'Goblet Cell', 'Ribosomal Cluster', 'Regulatory T Cell', 'NK/CD8 T Cell', 'CD8 T Cell', 'Ambiguous T Cell', 'Inflammatory Monocyte', 'TA Cell', 'Cycling B Cell', 'CD4 T Cell', 'Mesenchymal Cell', 'CD8 T Cell', 'Mesenchymal Cell', 'Endothelial Cell', 'Tuft Cell', 'Mast Cell', 'Enterocyte', 'Naive B Cell', 'Doublet Cluster', 'Enteroendocrine', 'ILC3', 'Plasmablast')
)

ileum_ref$cell_typev1 <- factor(ileum_ref$cell_typev1, levels = c("Naive B Cell 1", "Naive B Cell 2", "Memory B Cell", "Cycling B Cell", "Plasma Cell 1", "Plasma Cell 2", "Plasmablast", "Naive/CM CD4 T Cell", "Tissue Resident Memory CD4 T Cell", "Regulatory T Cell", "Naive/CM CD8 T Cell", "Tissue Resident Memory CD8 T Cell", "NK/CD8 T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte 1", "Enterocyte 2", "Enterocyte 3", "Enterocyte 4", "Stem/Paneth Cell 1", "Stem/Paneth Cell 2", "Stem Cell", "TA Cell", "Goblet Cell 1", "Goblet Cell 2", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell 1", "Mesenchymal Cell 2", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))

ileum_ref$cell_typev2 <- factor(ileum_ref$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/CD8 T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte", "Stem/Paneth Cell", "Stem Cell", "TA Cell", "Goblet Cell", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))

#create UMAPs

#blue = t cells 
blue <- randomColor(8, hue = c("blue"))
blue <- c("#34bfe5", "#4a93ce", "#b0daf2", "#6db5ed", "#1c5682","#077b9e", "#70c8e0", "#1f86f4") 

#green = b cells
green <- randomColor(7, hue = c("green"))
green <- c("#05fc95", "#9de587", "#1da54f", "#61e291", "#2bd147","#80cc1e", "#c5ef99")

#orange = myeloid 
orange <- randomColor(4, hue = c("orange"))
orange <- c("#fcab85", "#d3b050", "#f75f13", "#cc5414")

#pink = epithelial
pink <- randomColor(12, hue = c("pink"), luminosity = c("light"))
pink <- c("#f41f8d", "#ffb7de", "#ed6add", "#fd8cff", "#cc0076", "#ffbff9", "#f271d6", "#f279b5", "#f47cb8","#eec4fc", "#ff7cc6", "#ed6ab2")


#purple = mesenchymal and endothelial 
purple <- randomColor(3, hue = c("purple"))
purple <- c("#d3adea", "#7214ad", "#ab26ff")

#red - ribosomal and doublet cluster 
red <- randomColor(2, hue = c("red"))

red <- c("#e8480d", "#dd2a02")

DimPlot(ileum_ref, reduction = "umap", group.by = "cell_typev1", cols = c(green[1], green[2], green[3], green[4], green[5], green[6], green[7], blue[1], blue[2], blue[3], blue[4], blue[5], blue[6], blue[7], blue[8], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[4], pink[5], pink[6], pink[7], pink[8], pink[9], pink[10], pink[11], pink[12], purple[1], purple[2], purple[3], red[2], red[1]))

DimPlot(ileum_ref, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[3], green[4], green[5], green[7], blue[1], blue[3], blue[4], blue[6], blue[7], blue[8], orange[1], orange[2], orange[3], orange[4], pink[1], pink[5], pink[7], pink[8], pink[9], pink[11], pink[12], purple[1], purple[3], red[2], red[1]))

#stacked bar plots 

#proportion of cells - cell type v1
ileum_ref$orig.ident <- factor(ileum_ref$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum_ref) <- ileum_ref$cell_typev1
pt <- table(Idents(ileum_ref), ileum_ref$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell 1", "Naive B Cell 2", "Memory B Cell", "Cycling B Cell", "Plasma Cell 1", "Plasma Cell 2", "Plasmablast", "Naive/CM CD4 T Cell", "Tissue Resident Memory CD4 T Cell", "Regulatory T Cell", "Naive/CM CD8 T Cell", "Tissue Resident Memory CD8 T Cell", "NK/CD8 T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte 1", "Enterocyte 2", "Enterocyte 3", "Enterocyte 4", "Stem/Paneth Cell 1", "Stem/Paneth Cell 2", "Stem Cell", "TA Cell", "Goblet Cell 1", "Goblet Cell 2", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell 1", "Mesenchymal Cell 2", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))


plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[4], green[5], green[6], green[7], blue[1], blue[2], blue[3], blue[4], blue[5], blue[6], blue[7], blue[8], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[4], pink[5], pink[6], pink[7], pink[8], pink[9], pink[10], pink[11], pink[12], purple[1], purple[2], purple[3], red[2], red[1])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))


#proportion of cells - cell type v2
ileum_ref$orig.ident <- factor(ileum_ref$orig.ident, levels = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'))
Idents(ileum_ref) <- ileum_ref$cell_typev2
pt <- table(Idents(ileum_ref), ileum_ref$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/CD8 T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Enterocyte", "Stem/Paneth Cell", "Stem Cell", "TA Cell", "Goblet Cell", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))


plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[3], green[4], green[5], green[7], blue[1], blue[3], blue[4], blue[6], blue[7], blue[8], orange[1], orange[2], orange[3], orange[4], pink[1], pink[5], pink[7], pink[8], pink[9], pink[11], pink[12], purple[1], purple[3], red[2], red[1])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  # xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  # ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(size=10))



#---------------Robustness check---------------#

#use cell annotations 

#11/18/2024

#load in the reference cell annotations

ambig_t_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/ambiguous_t_cell.csv")
cd4_t_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/cd4_t_cell.csv")
cd8_t_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/cd8_t_cell.csv")
cycling_b_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/cycling_b_cell.csv")
doublet_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/doublet_cluster.csv")
endo_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/endothelial_cell.csv")
ent_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/enterocyte.csv")
entero_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/enteroendocrine.csv")
gob_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/goblet_cell.csv")
ilc3_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/ilc3.csv")
inflam_mono_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/inflammatory_monocyte.csv")
macro_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/macrophage.csv")
mast_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/mast_cell.csv")
mem_b_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/memory_b_cell.csv")
mes_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/mesenchymal_cell.csv")
mono_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/monocyte.csv")
naive_b_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/naive_b_cell.csv")
nk_t_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/nk_t_cell.csv")
plasma_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/plasma_cell.csv")
plasmablast_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/plasmablast.csv")
treg_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/regulatory_t_cell.csv")
ribo_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/ribosomal_cluster.csv")
stem_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/stem_cell.csv")
stem_pan_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/stem_paneth_cell.csv")
ta_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/ta_cell.csv")
tuft_ref <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/cell_annotation/param1_cell_typev2/tuft_cell.csv")

#reference 
deg_list_celltype_ref <- list(ambig_t_ref = ambig_t_ref, cd4_t_ref = cd4_t_ref, cd8_t_ref = cd8_t_ref, cycling_b_ref = cycling_b_ref, doublet_ref = doublet_ref, endo_ref = endo_ref, ent_ref = ent_ref, entero_ref = entero_ref, gob_ref = gob_ref, ilc3_ref = ilc3_ref, inflam_mono_ref = inflam_mono_ref, macro_ref = macro_ref, mast_ref = mast_ref, mem_b_ref = mem_b_ref, mes_ref = mes_ref, mono_ref = mono_ref, naive_b_ref = naive_b_ref, nk_t_ref = nk_t_ref, plasma_ref = plasma_ref, plasmablast_ref = plasmablast_ref, treg_ref = treg_ref, ribo_ref = ribo_ref, stem_ref = stem_ref, stem_pan_ref = stem_pan_ref, ta_ref = ta_ref, tuft_ref = tuft_ref)



#load in low resolution cell annotation 

ambig_t_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/ambiguous_t_cell.csv")
cd4_t_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/cd4_t_cell.csv")
cd8_t_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/cd8_t_cell.csv")
doublet_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/doublet_cluster.csv")
endo_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/endothelial_cell.csv")
ent_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/enterocyte.csv")
gob_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/goblet_cell.csv")
ilc3_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/ILC3.csv")
inflam_mono_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/inflammatory_monocyte.csv")
macro_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/macrophage.csv")
mast_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/mast_cell.csv")
mem_b_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/memory_b_cell.csv")
mes_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/mesenchymal_cell.csv")
mono_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/monocyte.csv")
naive_b_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/naive_b_cell.csv")
nk_t_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/nk_t_cell.csv")
plasma_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/plasma_cell.csv")
plasmablast_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/plasmablast.csv")
stem_pan_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/stem_paneth_cell.csv")
ta_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/ta_cell.csv")
tuft_low <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/cell_annotation/param1_cell_typev2/tuft_cell.csv")


#low res 
deg_list_celltype_low <- list(ambig_t_low = ambig_t_low, cd4_t_low = cd4_t_low, cd8_t_low = cd8_t_low, doublet_low = doublet_low, endo_low = endo_low, ent_low = ent_low, gob_low = gob_low, ilc3_low = ilc3_low, inflam_mono_low = inflam_mono_low, macro_low = macro_low, mast_low = mast_low, mem_b_low = mem_b_low, mes_low = mes_low, mono_low = mono_low, naive_b_low = naive_b_low, nk_t_low = nk_t_low, plasma_low = plasma_low, plasmablast_low = plasmablast_low, stem_pan_low = stem_pan_low, ta_low = ta_low, tuft_low = tuft_low)


#load in high resolution clustering results 

ambig_t_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/ambiguous_t_cell.csv")
cd4_t_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/cd4_t_cell.csv")
cd8_t_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/cd8_t_cell.csv")
cyc_nk_t_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/cycling_nk_t_cell.csv")
doublet_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/doublet_cluster.csv")
endo_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/endothelial_cell.csv")
ent_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/enterocyte.csv")
entero_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/enteroendocrine.csv")
gc_b_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/gc_b_cell.csv")
gob_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/goblet_cell.csv")
ilc3_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/ilc3.csv")
inflam_macro_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/inflammatory_macrophage.csv")
inflam_mono_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/inflammatory_monocyte.csv")
macro_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/macrophage.csv")
mast_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/mast_cell.csv")
mem_b_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/memory_b_cell.csv")
mes_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/mesenchymal_cell.csv")
mono_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/monocyte.csv")
naive_b_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/naive_b_cell.csv")
nk_t_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/nk_t_cell.csv")
plasma_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/plasma_cell.csv")
plasmablast_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/plasmablast.csv")
treg_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/regulatory_t_cell.csv")
ribo_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/ribosomal_cluster.csv")
stem_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/stem_cell.csv")
stem_pan_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/stem_paneth_cell.csv")
ta_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/ta_cell.csv")
tuft_high <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/cell_annotation/param1_cell_typev2/tuft_cell.csv")


#high res 
deg_list_celltype_high <- list(ambig_t_high = ambig_t_high, cd4_t_high = cd4_t_high, cd8_t_high = cd8_t_high, cyc_nk_t_high = cyc_nk_t_high, doublet_high = doublet_high, endo_high = endo_high, ent_high = ent_high, entero_high = entero_high, gc_b_high = gc_b_high, gob_high = gob_high, ilc3_high = ilc3_high, inflam_macro_high = inflam_macro_high, inflam_mono_high = inflam_mono_high, macro_high = macro_high, mast_high = mast_high, mem_b_high = mem_b_high, mes_high = mes_high, mono_high = mono_high, naive_b_high = naive_b_high, nk_t_high = nk_t_high, plasma_high = plasma_high, plasmablast_high = plasmablast_high, treg_high = treg_high, ribo_high = ribo_high, stem_high = stem_high, stem_pan_high = stem_pan_high, tuft_high = tuft_high, ta_high = ta_high)

Idents(ileum_ref) <- ileum_ref$cell_typev2
Idents(ileum_low) <- ileum_low$cell_typev2
Idents(ileum_high) <- ileum_high$cell_typev2


#----------reference vs low res results-----------#


# vectors A and B 
jaccard_fun <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


clusters_ref <- unique(Idents(ileum_ref))
clusters_low <- unique(Idents(ileum_low))
clusters_high <- unique(Idents(ileum_high))

# Initialize a matrix to store the results
jaccard_matrix <- matrix(NA, nrow = length(clusters_ref), ncol = length(clusters_low),
                         dimnames = list(paste0("Ref_", clusters_ref), paste0("Low_", clusters_low)))

# Loop over each cluster in the reference and query sets
for (ref_cluster in clusters_ref) {
  for (low_cluster in clusters_low) {
    # Get cells in each cluster
    ref_cells <- WhichCells(ileum_ref, idents = as.character(ref_cluster))
    low_cells <- WhichCells(ileum_low, idents = as.character(low_cluster))
    
    # Calculate Jaccard index for the cluster pair
    jaccard_matrix[paste0("Ref_", ref_cluster), paste0("Low_", low_cluster)] <- jaccard_fun(ref_cells, low_cells)
  }
}

pheatmap(jaccard_matrix,
         main = "Jaccard Index Heatmap Cell Annotations: Reference vs. Low Resolution",
         color = colorRampPalette(c("white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE)  # Optionally display the actual Jaccard index values


# Initialize a matrix to store the results high res vs. reference 
jaccard_matrix_barcode_high_ref <- matrix(NA, nrow = length(clusters_ref), ncol = length(clusters_high),
                                          dimnames = list(paste0("Ref_", clusters_ref), paste0("High_", clusters_high)))

# Loop over each cluster in the reference and query sets
for (ref_cluster in clusters_ref) {
  for (high_cluster in clusters_high) {
    # Get cells in each cluster
    ref_cells <- WhichCells(ileum_ref, idents = as.character(ref_cluster))
    high_cells <- WhichCells(ileum_high, idents = as.character(high_cluster))
    
    # Calculate Jaccard index for the cluster pair
    jaccard_matrix_barcode_high_ref[paste0("Ref_", ref_cluster), paste0("High_", high_cluster)] <- jaccard_fun(ref_cells, high_cells)
  }
}

pheatmap(jaccard_matrix_barcode_high_ref,
         main = "Jaccard Index Heatmap Cell Annotations: Reference vs. High Resolution",
         color = colorRampPalette(c("white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE)  # Optionally display the actual Jaccard index values


# Initialize a matrix to store the results high res vs. low res 
jaccard_matrix_barcode_high_low <- matrix(NA, nrow = length(clusters_low), ncol = length(clusters_high),
                                          dimnames = list(paste0("Low_", clusters_low), paste0("High_", clusters_high)))

# Loop over each cluster in the reference and query sets
for (low_cluster in clusters_low) {
  for (high_cluster in clusters_high) {
    # Get cells in each cluster
    low_cells <- WhichCells(ileum_low, idents = as.character(low_cluster))
    high_cells <- WhichCells(ileum_high, idents = as.character(high_cluster))
    
    # Calculate Jaccard index for the cluster pair
    jaccard_matrix_barcode_high_low[paste0("Low_", low_cluster), paste0("High_", high_cluster)] <- jaccard_fun(low_cells, high_cells)
  }
}

pheatmap(jaccard_matrix_barcode_high_low,
         main = "Jaccard Index Heatmap Cell Annotations: High Resolution vs. Low Resolution",
         color = colorRampPalette(c("white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE)  # Optionally display the actual Jaccard index values

#jaccard index based on genes present in each FM result

clusters_ref <- c("ambig_t_ref", "cd4_t_ref", 'cd8_t_ref', 'cycling_b_ref', 'doublet_ref', 'endo_ref', 'ent_ref', 'entero_ref', 'gob_ref', 'ilc3_ref', 'inflam_mono_ref', 'macro_ref', 'mast_ref', 'mem_b_ref', 'mes_ref', 'mono_ref', 'naive_b_ref', 'nk_t_ref', 'plasma_ref', 'plasmablast_ref', 'treg_ref', 'ribo_ref', 'stem_ref', 'stem_pan_ref', 'ta_ref', 'tuft_ref')
clusters_low <- c('ambig_t_low', 'cd4_t_low', 'cd8_t_low', 'doublet_low', 'endo_low', 'ent_low', 'gob_low', 'ilc3_low', 'inflam_mono_low', 'macro_low', 'mast_low', 'mem_b_low', 'mes_low', 'mono_low', 'naive_b_low', 'nk_t_low', 'plasma_low', 'plasmablast_low', 'stem_pan_low', 'ta_low', 'tuft_low')
clusters_high <- c('ambig_t_high', 'cd4_t_high', 'cd8_t_high', 'cyc_nk_t_high', 'doublet_high', 'endo_high', 'ent_high', 'entero_high', 'gc_b_high', 'gob_high', 'ilc3_high', 'inflam_macro_high', 'inflam_mono_high', 'macro_high', 'mast_high', 'mem_b_high', 'mes_high', 'mono_high', 'naive_b_high', 'nk_t_high', 'plasma_high', 'plasmablast_high', 'treg_high', 'ribo_high', 'stem_high', 'stem_pan_high', 'tuft_high', 'ta_high')

#low res vs. high res jaccard index with genes 

# Initialize a matrix to store the results
jaccard_matrix_low_ref_gene <- matrix(NA, nrow = length(clusters_ref), ncol = length(clusters_low),
                                      dimnames = list(paste0("Ref_", clusters_ref), paste0("Low_", clusters_low)))


# Loop through each pair of clusters
for (ref_cluster in clusters_ref) {
  for (low_cluster in clusters_low) {
    
    # Extract DEG results for each cluster
    ref_data <- deg_list_celltype_ref[[ref_cluster]]
    low_data <- deg_list_celltype_low[[low_cluster]]
    
    ref_data_gene <- ref_data$gene
    low_data_gene <- low_data$gene
    
    # Calculate Jaccard index for the cluster pair
    jaccard_matrix_low_ref_gene[paste0("Ref_", ref_cluster), paste0("Low_", low_cluster)] <- jaccard_fun(ref_data_gene, low_data_gene)
    
  }
}

pheatmap(jaccard_matrix_low_ref_gene,
         main = "Jaccard Index Heatmap Cell Annotations FM Gene: Reference vs. Low Resolution",
         color = colorRampPalette(c("white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE)  # Optionally display the actual Jaccard index values

#high res vs. ref jaccard index with genes 

# Initialize a matrix to store the results
jaccard_matrix_high_ref_gene <- matrix(NA, nrow = length(clusters_ref), ncol = length(clusters_high),
                                       dimnames = list(paste0("Ref_", clusters_ref), paste0("High_", clusters_high)))


# Loop through each pair of clusters
for (ref_cluster in clusters_ref) {
  for (high_cluster in clusters_high) {
    
    # Extract DEG results for each cluster
    ref_data <- deg_list_celltype_ref[[ref_cluster]]
    high_data <- deg_list_celltype_high[[high_cluster]]
    
    ref_data_gene <- ref_data$gene
    high_data_gene <- high_data$gene
    
    # Calculate Jaccard index for the cluster pair
    jaccard_matrix_high_ref_gene[paste0("Ref_", ref_cluster), paste0("High_", high_cluster)] <- jaccard_fun(ref_data_gene, high_data_gene)
    
  }
}

pheatmap(jaccard_matrix_high_ref_gene,
         main = "Jaccard Index Heatmap Cell Annotations FM Gene: Reference vs. High Resolution",
         color = colorRampPalette(c("white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE)  # Optionally display the actual Jaccard index values

#high res vs. low res jaccard index with genes 

# Initialize a matrix to store the results
jaccard_matrix_high_low_gene <- matrix(NA, nrow = length(clusters_low), ncol = length(clusters_high),
                                       dimnames = list(paste0("Low_", clusters_low), paste0("High_", clusters_high)))


# Loop through each pair of clusters
for (low_cluster in clusters_low) {
  for (high_cluster in clusters_high) {
    
    # Extract DEG results for each cluster
    low_data <- deg_list_celltype_low[[low_cluster]]
    high_data <- deg_list_celltype_high[[high_cluster]]
    
    low_data_gene <- low_data$gene
    high_data_gene <- high_data$gene
    
    # Calculate Jaccard index for the cluster pair
    jaccard_matrix_high_low_gene[paste0("Low_", low_cluster), paste0("High_", high_cluster)] <- jaccard_fun(low_data_gene, high_data_gene)
    
  }
}

pheatmap(jaccard_matrix_high_low_gene,
         main = "Jaccard Index Heatmap Cell Annotations FM Gene: Low Resolution vs. High Resolution",
         color = colorRampPalette(c("white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE)  # Optionally display the actual Jaccard index values


#-----reference and low res-------#
clusters_ref <- c("ambig_t_ref", "cd4_t_ref", 'cd8_t_ref', 'cycling_b_ref', 'doublet_ref', 'endo_ref', 'ent_ref', 'entero_ref', 'gob_ref', 'ilc3_ref', 'inflam_mono_ref', 'macro_ref', 'mast_ref', 'mem_b_ref', 'mes_ref', 'mono_ref', 'naive_b_ref', 'nk_t_ref', 'plasma_ref', 'plasmablast_ref', 'treg_ref', 'ribo_ref', 'stem_ref', 'stem_pan_ref', 'ta_ref', 'tuft_ref')
clusters_low <- c('ambig_t_low', 'cd4_t_low', 'cd8_t_low', 'doublet_low', 'endo_low', 'ent_low', 'gob_low', 'ilc3_low', 'inflam_mono_low', 'macro_low', 'mast_low', 'mem_b_low', 'mes_low', 'mono_low', 'naive_b_low', 'nk_t_low', 'plasma_low', 'plasmablast_low', 'stem_pan_low', 'ta_low', 'tuft_low')
clusters_high <- c('ambig_t_high', 'cd4_t_high', 'cd8_t_high', 'cyc_nk_t_high', 'doublet_high', 'endo_high', 'ent_high', 'entero_high', 'gc_b_high', 'gob_high', 'ilc3_high', 'inflam_macro_high', 'inflam_mono_high', 'macro_high', 'mast_high', 'mem_b_high', 'mes_high', 'mono_high', 'naive_b_high', 'nk_t_high', 'plasma_high', 'plasmablast_high', 'treg_high', 'ribo_high', 'stem_high', 'stem_pan_high', 'tuft_high', 'ta_high')


# Initialize a matrix to store the results
spearman_matrix_low_ref_log2fc <- matrix(NA, nrow = length(clusters_low), ncol = length(clusters_ref),
                                         dimnames = list(paste0("Low_", clusters_low), paste0("Ref_", clusters_ref)))


# Loop through each pair of clusters
for (ref_cluster in clusters_ref) {
  for (low_cluster in clusters_low) {
    
    # Extract DEG results for each cluster
    ref_data <- deg_list_celltype_ref[[ref_cluster]]
    low_data <- deg_list_celltype_low[[low_cluster]]
    
    
    # Ensure that the data frames contain only the genes present in both clusters
    common_genes <- intersect(ref_data$gene, low_data$gene)
    
    # Subset the dataframes to only include the common genes
    ref_subset <- ref_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    low_subset <- low_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    
    #check if genes are in the same order 
    print(all(ref_subset$gene == low_subset$gene))
    # Check if there are any common genes
    if (length(common_genes) > 0) {
      # Calculate Spearman correlation for log2FC values of the common genes
      spearman_corr <- cor(ref_subset$avg_log2FC, low_subset$avg_log2FC, method = "spearman")
      
      # Store the correlation in the matrix
      spearman_matrix_low_ref_log2fc[paste0("Low_", low_cluster), paste0("Ref_", ref_cluster)] <- spearman_corr
    } else {
      # If no common genes, set correlation to NA
      spearman_matrix_low_ref_log2fc[paste0("Low_", low_cluster), paste0("Ref_", ref_cluster)] <- NA
    }
  }
}

pheatmap(spearman_matrix_low_ref_log2fc,
         main = "Average Log2FC Spearman Correlation: Low Resolution Stabley Assigned Cells vs. Reference",
         color = colorRampPalette(c("purple", "white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE, 
         number_color = "black")  # Optionally display the actual Jaccard index values


# Initialize a matrix to store the results - high res vs. reference 
spearman_matrix_high_ref_log2fc <- matrix(NA, nrow = length(clusters_high), ncol = length(clusters_ref),
                                          dimnames = list(paste0("High_", clusters_high), paste0("Ref_", clusters_ref)))


# Loop through each pair of clusters
for (ref_cluster in clusters_ref) {
  for (high_cluster in clusters_high) {
    
    # Extract DEG results for each cluster
    ref_data <- deg_list_celltype_ref[[ref_cluster]]
    high_data <- deg_list_celltype_high[[high_cluster]]
    
    
    # Ensure that the data frames contain only the genes present in both clusters
    common_genes <- intersect(ref_data$gene, high_data$gene)
    
    # Subset the dataframes to only include the common genes
    ref_subset <- ref_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    high_subset <- high_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    
    #check if genes are in the same order 
    print(all(ref_subset$gene == high_subset$gene))
    # Check if there are any common genes
    if (length(common_genes) > 0) {
      # Calculate Spearman correlation for log2FC values of the common genes
      spearman_corr <- cor(ref_subset$avg_log2FC, high_subset$avg_log2FC, method = "spearman")
      
      # Store the correlation in the matrix
      spearman_matrix_high_ref_log2fc[paste0("High_", high_cluster), paste0("Ref_", ref_cluster)] <- spearman_corr
    } else {
      # If no common genes, set correlation to NA
      spearman_matrix_high_ref_log2fc[paste0("High_", high_cluster), paste0("Ref_", ref_cluster)] <- NA
    }
  }
}

pheatmap(spearman_matrix_high_ref_log2fc,
         main = "Average Log2FC Spearman Correlation: High Resolution Stabley Assigned Cells vs. Reference",
         color = colorRampPalette(c("purple", "white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE, 
         number_color = "black")  # Optionally display the actual Jaccard index values

# Initialize a matrix to store the results - high res vs. low res 
spearman_matrix_high_low_log2fc <- matrix(NA, nrow = length(clusters_high), ncol = length(clusters_low),
                                          dimnames = list(paste0("High_", clusters_high), paste0("Low_", clusters_low)))


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
      spearman_matrix_high_low_log2fc[paste0("High_", high_cluster), paste0("Low_", low_cluster)] <- spearman_corr
    } else {
      # If no common genes, set correlation to NA
      spearman_matrix_high_low_log2fc[paste0("High_", high_cluster), paste0("Low_", low_cluster)] <- NA
    }
  }
}

pheatmap(spearman_matrix_high_low_log2fc,
         main = "Average Log2FC Spearman Correlation: High Resolution vs. Low Resolution Stabley Assigned Cells",
         color = colorRampPalette(c("purple", "white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE, 
         number_color = "black")  # Optionally display the actual Jaccard index values

#compute the pairwise correlation between adj p values 

# Initialize a matrix to store the results
spearman_matrix_low_ref_pval <- matrix(NA, nrow = length(clusters_low), ncol = length(clusters_ref),
                                       dimnames = list(paste0("Low_", clusters_low), paste0("Ref_", clusters_ref)))


# Loop through each pair of clusters
for (ref_cluster in clusters_ref) {
  for (low_cluster in clusters_low) {
    
    # Extract DEG results for each cluster
    ref_data <- deg_list_celltype_ref[[ref_cluster]]
    low_data <- deg_list_celltype_low[[low_cluster]]
    
    
    # Ensure that the data frames contain only the genes present in both clusters
    common_genes <- intersect(ref_data$gene, low_data$gene)
    
    # Subset the dataframes to only include the common genes
    ref_subset <- ref_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    low_subset <- low_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    
    #check if genes are in the same order 
    print(all(ref_subset$gene == low_subset$gene))
    # Check if there are any common genes
    if (length(common_genes) > 0) {
      # Calculate Spearman correlation for log2FC values of the common genes
      spearman_corr <- cor(ref_subset$p_val_adj, low_subset$p_val_adj, method = "spearman")
      
      # Store the correlation in the matrix
      spearman_matrix_low_ref_pval[paste0("Low_", low_cluster), paste0("Ref_", ref_cluster)] <- spearman_corr
    } else {
      # If no common genes, set correlation to NA
      spearman_matrix_low_ref_pval[paste0("Low_", low_cluster), paste0("Ref_", ref_cluster)] <- NA
    }
  }
}

pheatmap(spearman_matrix_low_ref_pval,
         main = "Adj. P value Spearman Correlation: Low Resolution Stabley Assigned Cells vs. Reference",
         color = colorRampPalette(c("purple", "white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE, 
         number_color = "black")  # Optionally display the actual Jaccard index values

#high res vs reference 
# Initialize a matrix to store the results
spearman_matrix_high_ref_pval <- matrix(NA, nrow = length(clusters_high), ncol = length(clusters_ref),
                                        dimnames = list(paste0("High_", clusters_high), paste0("Ref_", clusters_ref)))


# Loop through each pair of clusters
for (ref_cluster in clusters_ref) {
  for (high_cluster in clusters_high) {
    
    # Extract DEG results for each cluster
    ref_data <- deg_list_celltype_ref[[ref_cluster]]
    high_data <- deg_list_celltype_high[[high_cluster]]
    
    
    # Ensure that the data frames contain only the genes present in both clusters
    common_genes <- intersect(ref_data$gene, high_data$gene)
    
    # Subset the dataframes to only include the common genes
    ref_subset <- ref_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    high_subset <- high_data %>% filter(gene %in% common_genes) %>% arrange(gene)
    
    #check if genes are in the same order 
    print(all(ref_subset$gene == high_subset$gene))
    # Check if there are any common genes
    if (length(common_genes) > 0) {
      # Calculate Spearman correlation for log2FC values of the common genes
      spearman_corr <- cor(ref_subset$p_val_adj, high_subset$p_val_adj, method = "spearman")
      
      # Store the correlation in the matrix
      spearman_matrix_high_ref_pval[paste0("High_", high_cluster), paste0("Ref_", ref_cluster)] <- spearman_corr
    } else {
      # If no common genes, set correlation to NA
      spearman_matrix_high_ref_pval[paste0("High_", high_cluster), paste0("Ref_", ref_cluster)] <- NA
    }
  }
}

pheatmap(spearman_matrix_high_ref_pval,
         main = "Adj. P value Spearman Correlation: High Resolution Stabley Assigned Cells vs. Reference",
         color = colorRampPalette(c("purple", "white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE, 
         number_color = "black")  # Optionally display the actual Jaccard index values

#high res vs low res  
# Initialize a matrix to store the results
spearman_matrix_high_low_pval <- matrix(NA, nrow = length(clusters_high), ncol = length(clusters_low),
                                        dimnames = list(paste0("High_", clusters_high), paste0("Low_", clusters_low)))


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
      spearman_corr <- cor(low_subset$p_val_adj, high_subset$p_val_adj, method = "spearman")
      
      # Store the correlation in the matrix
      spearman_matrix_high_low_pval[paste0("High_", high_cluster), paste0("Low_", low_cluster)] <- spearman_corr
    } else {
      # If no common genes, set correlation to NA
      spearman_matrix_high_low_pval[paste0("High_", high_cluster), paste0("Low_", low_cluster)] <- NA
    }
  }
}

pheatmap(spearman_matrix_high_low_pval,
         main = "Adj. P value Spearman Correlation: High Resolution vs. Low Resolution Stabley Assigned Cells",
         color = colorRampPalette(c("purple", "white", "red"))(50), # Color gradient from white to blue
         cluster_rows = TRUE,  # Disable clustering if you don't want it
         cluster_cols = TRUE,
         display_numbers = TRUE, 
         number_color = "black")  # Optionally display the actual Jaccard index values

#-------------create Urko Barplot-------------#

#11/20/2024

library(GO.db)
library(biomaRt)

library(clusterProfiler)
library(org.Hs.eg.db)  # Use appropriate annotation package for your organism

#load in seurat object 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_scITD_groupInfo_11_20_24.rds")

#do this on server 
#set idents to group 
Idents(ileum) <- ileum$group

ileum_counts <- AverageExpression(ileum, assays = "RNA", group.by = "group", layer = "data")

write.csv(ileum_counts, file = "/storage/home/swashburn30/Helmsley/Batch1_13/DEG_analysis/helm_batch1_13_ileum_group_avg_counts.csv")

#---read in counts----# 
ileum_counts <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_group_avg_counts.csv")

rownames(ileum_counts) <- ileum_counts$X

colnames(ileum_counts) <- c("genes", "group_1", "group_2")


# Add a small constant (e.g., 1) to handle zeros
ileum_counts$log2FC <- log2((ileum_counts$group_1 + 1) / (ileum_counts$group_2 + 1))

#positive is up in group 1
#negative is up in group 2

# Assuming your dataframe is called df and the column you are checking is called column_name
#ileum_counts$sig <- ifelse(ileum_counts$log2FC > 0, 1, -1)

ileum_counts$log2FC <- as.numeric(ileum_counts$log2FC)

#ileum_counts$sig_strict <- ifelse(ileum_counts$log2FC > 1, 1, ifelse(ileum_counts$log2FC < -1, -1, NA))

#use this threshold (1 was too high)
ileum_counts$sig_strict <- ifelse(ileum_counts$log2FC > 0.25, 1,
                                  ifelse(ileum_counts$log2FC < -0.25, -1, NA))


# Example GO term: "GO:0006955" (immune response)

#metallopeptidase activity 
go_term <- "GO:0008237"  # Replace with your GO term ID


# Convert GO term to gene symbols
genes_GO0008237 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# View the resulting gene symbols
#head(genes)
# Connect to the Ensembl biomart
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve genes associated with the GO term
#genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
#filters = 'go_id', values = go_term, mart = ensembl)

intersect(ileum_counts$genes, genes_GO0008237$SYMBOL) #170 

overlap <- ileum_counts[ileum_counts$gene %in% genes_GO0008237$SYMBOL,]
table(overlap$sig_strict) #4, 18

#response to TNF
go_term <- "GO:0034612"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0034612 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db)


intersect(ileum_counts$genes, genes_GO0034612$SYMBOL) #235

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0034612$SYMBOL,]
table(overlap$sig_strict) #18,17

#cytokine activity 

#GO:0005125

go_term <- "GO:0005125"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0005125 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db)


intersect(ileum_counts$genes, genes_GO0005125$SYMBOL) #203

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0005125$SYMBOL,]
table(overlap$sig_strict) #13,5

#response to LPS 

#GO:0032496

go_term <- "GO:0032496"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0032496 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db)


intersect(ileum_counts$genes, genes_GO0032496$SYMBOL) #315

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0032496$SYMBOL,]
table(overlap$sig_strict) #22,14 

#collagen binding 

#GO:0005518

go_term <- "GO:0005518"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0005518 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #339


intersect(ileum_counts$genes, genes_GO0005518$SYMBOL) #67

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0005518$SYMBOL,]
table(overlap$sig_strict) #3,2

#inflammatory response 

go_term <- "GO:0006954"

# Convert GO term to gene symbols
genes_GO0006954 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #339


intersect(ileum_counts$genes, genes_GO0006954$SYMBOL) #752

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0006954$SYMBOL,]
table(overlap$sig_strict) #50,34

#collagen catabolic process 

go_term <- "GO:0030574"

genes_GO0030574 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #339


intersect(ileum_counts$genes, genes_GO0030574$SYMBOL) #43

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0030574$SYMBOL,]
table(overlap$sig_strict) #1,3

#neutrophil chemotaxis 

#GO:0030593

go_term <- "GO:0030593"

genes_GO0030593 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #106


intersect(ileum_counts$genes, genes_GO0030593$SYMBOL) #102

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0030593$SYMBOL,]
table(overlap$sig_strict) #20,6

#neutrophil migration 

#GO:1990266

go_term <- "GO:1990266"

genes_GO1990266 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #130


intersect(ileum_counts$genes, genes_GO1990266$SYMBOL) #124

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO1990266$SYMBOL,]
table(overlap$sig_strict) #20,6

#chemokine receptor binding 

#GO:0042379

go_term <- "GO:0042379"

genes_GO0042379 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #71


intersect(ileum_counts$genes, genes_GO0042379$SYMBOL) #59

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0042379$SYMBOL,]
table(overlap$sig_strict) #7,1

#response to bacterium 

#GO:0009617

go_term <- "GO:0009617"

genes_GO0009617 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #786


intersect(ileum_counts$genes, genes_GO0009617$SYMBOL) #668

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0009617$SYMBOL,]
table(overlap$sig_strict) #56, 40

#response to wounding 

#GO:0009611

go_term <- "GO:0009611"

genes_GO0009611 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #577


intersect(ileum_counts$genes, genes_GO0009611$SYMBOL) #549

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0009611$SYMBOL,]
table(overlap$sig_strict) #22,42

#cytokine receptor binding 

#GO:0005126

go_term <- "GO:0005126"

genes_GO0005126 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #273


intersect(ileum_counts$genes, genes_GO0005126$SYMBOL) #235

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0005126$SYMBOL,]
table(overlap$sig_strict) #13,7

#ECM disassembly 

#GO:0022617

go_term <- "GO:0022617"

genes_GO0022617 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #62


intersect(ileum_counts$genes, genes_GO0022617$SYMBOL) #60

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0022617$SYMBOL,]
table(overlap$sig_strict) #3,3


#growth factor binding 

#GO:0019838

go_term <- "GO:0019838"

genes_GO0019838 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #132


intersect(ileum_counts$genes, genes_GO0019838$SYMBOL) #128

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0019838$SYMBOL,]
table(overlap$sig_strict) #7,3


#ecm orangization 

#GO:0030198

go_term <- "GO:0030198"

genes_GO0030198 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #318


intersect(ileum_counts$genes, genes_GO0030198$SYMBOL) #303

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0030198$SYMBOL,]
table(overlap$sig_strict) #3,14


#fibronectin binding 

#GO:0001968

go_term <- "GO:0001968"

genes_GO0001968 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #31


intersect(ileum_counts$genes, genes_GO0001968$SYMBOL) #29

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0001968$SYMBOL,]
table(overlap$sig_strict) #2,2

#GO:0005201

go_term <- "GO:0005201"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0005201 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db)


intersect(ileum_counts$genes, genes_GO0005201$SYMBOL) #159

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0005201$SYMBOL,]
table(overlap$sig) #2,7

#-------read in pathways for bar chart--------#

#load in file that has pathways and number of genes in each group 
bar_plot_pathways <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_group_pathway_enrichment_11_20_24.csv")


test <- bar_plot_pathways %>%
  group_by(pathways) %>%
  mutate(TotalGenes = sum(number)) %>%   # Calculate total number of genes in each pathway
  mutate(PercentGenes = (number / TotalGenes) * 100) %>%   # Calculate percentage of genes
  ungroup()

test <- test %>%
  mutate(PercentGenes = ifelse(group == "group 1",
                               PercentGenes,
                               -1*PercentGenes))

test_reorder <- test %>%
  group_by(group) %>%
  arrange(desc(abs(PercentGenes)), .by_group = TRUE) %>%  # Arrange by PercentGenes within each Pathway
  ungroup()  # Ungroup after sorting

## calculate breaks values
breaks_values <- pretty(test$PercentGenes)

plot <- test %>%
  ggplot(aes(x = factor(pathways, levels = c("Chemokine Receptor Binding (59)", "Neutrophil Chemotaxis (102)", "Neutrophil Migration (124)", "Cytokine Activity (203)", "Growth Factor Binding (128)", "Cytokine Receptor Binding (235)",  "Response to LPS (315)", "Collagen Binding (67)", "Inflammatory Response (752)", "Response to Bacterium (668)", "Response to TNF (235)", "ECM Disassembly (60)", "Fibronectin Binding (29)", "Response to Wounding (549)", "Collagen Catabolic Process (43)", "ECM Structural Constituent (159)", "Metallopeptidase activity (170)", "ECM Organization (303)")), y = PercentGenes, fill = group))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c('#F8766D', '#619CFF')) +
  coord_flip()+
  scale_y_continuous(breaks = breaks_values,
                     labels = abs(breaks_values)) + 
  geom_col(position = "dodge") +
  xlab("Pathways") + ylab("Percent of Genes in Pathway") +
  theme_minimal()+
  theme(axis.text.x = element_text(color="black", size=14),axis.text.y = element_text(color="black", size=14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  guides(fill=guide_legend(title="Groups")) 
#geom_text(aes(label=abs(number))) 
plot

#---------module score with Ucell----------#

library(UCell)

#11/22/2024

ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_scITD_groupInfo_11_20_24.rds")

ileum_counts <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_group_avg_counts.csv")

DefaultAssay(ileum) <- "RNA"

ileum$group <- factor(ileum$group, levels = c("group1", "group2"))

ileum@meta.data$cell_typev2 <- ileum@meta.data$seurat_clusters
ileum$cell_typev2 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21'),
  to = c("Naive B Cell", 'Stem/Paneth Cell', 'Enterocyte', 'CD4 T Cell', 'NK/CD8 T Cell', 'Enterocyte', 'Macrophage', "Plasma Cell", "Plasma Cell", 'Monocyte', 'Enterocyte', 'Goblet Cell', 'Goblet Cell', 'Ambig. T Cell', 'Mesenchymal Cell', "Pro-Inflammatory Monocyte", "Cycling Cell", "Memory B Cell", "Endothelial Cell", 'Mast Cell', "Tuft Cell", "Cycling Plasma Cell")
)


ileum$cell_typev2 <- factor(ileum$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Cycling Plasma Cell", "CD4 T Cell", "NK/CD8 T Cell", "Ambig. T Cell", "Macrophage", "Monocyte", "Pro-Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Enterocyte", "Goblet Cell", "Cycling Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))


#ECM organization 

ecm_genes <- intersect(ileum_counts$genes, genes_GO0030198$SYMBOL) #303

ecm_genes <- paste(shQuote(ecm_genes), collapse=", ")

ecm_genes <- c('VWA1', 'MMP23B', 'ANGPTL7', 'TNFRSF1B', 'PDPN', 'MATN1', 'COL16A1', 'COL8A2', 'COL9A2', 'TIE1', 'CCN1', 'COL24A1', 'COL11A1', 'NTNG1', 'ADAMTSL4', 'CTSS', 'CTSK', 'PBXIP1', 'ADAM15', 'HAPLN2', 'ADAMTS4', 'OLFML2B', 'DDR2', 'DPT', 'TNR', 'QSOX1', 'LAMC1', 'COLGALT2', 'HMCN1', 'ELF3', 'FMOD', 'LAMB3', 'TGFB2', 'WNT3A', 'AGT', 'EXOC8', 'NID1', 'PXDN', 'MATN3', 'EMILIN1', 'MPV17', 'VIT', 'CYP1B1', 'ANTXR1', 'LOXL3', 'CLASP1', 'DPP4', 'FAP', 'COL3A1', 'COL5A2', 'CFLAR', 'IHH', 'COL4A4', 'COL4A3', 'FBLN2', 'COLQ', 'CRTAP', 'CLASP2', 'KIF9', 'LAMB2', 'DAG1', 'ADAMTS9', 'COL8A1', 'ABI3BP', 'IMPG2', 'PHLDB2', 'CCDC80', 'COL6A5', 'COL6A6', 'NPHP3', 'PTX3', 'MELTF', 'PDGFRA', 'ADAMTS3', 'SLC39A8', 'NPNT', 'PRDM5', 'NDNF', 'SFRP2', 'RXFP1', 'TLL1', 'ADAMTS16', 'ADAMTS12', 'EGFLAM', 'ADAMTS6', 'LOX', 'ADAMTS19', 'TGFBI', 'SPINK5', 'SH3PXD2B', 'FGFR4', 'COL23A1', 'ADAMTS2', 'FOXF2', 'FOXC1', 'ADTRP', 'FLOT1', 'DDR1', 'TNF', 'TNXB', 'COL11A2', 'COL19A1', 'COL9A1', 'COL12A1', 'IMPG1', 'COL10A1', 'LAMA2', 'CCN2', 'SERAC1', 'PLG', 'SMOC2', 'FSCN1', 'COL28A1', 'IL6', 'AEBP1', 'ELN', 'COL1A2', 'PLOD3', 'LAMB1', 'LAMB4', 'CAV2', 'CAV1', 'ST7', 'PRSS1', 'PRSS2', 'DNAJB6', 'BMP1', 'LOXL2', 'SULF1', 'MMP16', 'MATN2', 'EXT1', 'TNFRSF11B', 'COL14A1', 'HAS2', 'COL22A1', 'SCX', 'WASHC1', 'RIC1', 'ADAMTSL1', 'B4GALT1', 'RECK', 'ECM2', 'CTSV', 'COL15A1', 'TGFBR1', 'TMEM38B', 'COL27A1', 'OLFML2A', 'ENG', 'ABL1', 'POMT1', 'NTNG2', 'ADAMTS13', 'ADAMTSL2', 'COL5A1', 'NOTCH1', 'ITGA8', 'ITGB1', 'COL13A1', 'ADAMTS14', 'SPOCK2', 'P4HA1', 'TLL2', 'LOXL4', 'HPSE2', 'KAZALD1', 'NFKB2', 'COL17A1', 'MMP21', 'ADAM8', 'RIC8A', 'GAS2', 'WT1', 'HSD17B12', 'CREB3L1', 'LTBP3', 'EFEMP2', 'SERPINH1', 'MMP7', 'MMP20', 'MMP27', 'MMP8', 'MMP10', 'MMP1', 'MMP3', 'MMP12', 'MMP13', 'MPZL3', 'PHLDB1', 'ADAMTS8', 'ADAMTS15', 'TNFRSF1A', 'COL2A1', 'MMP19', 'LRP1', 'LUM', 'NTN4', 'MMP17', 'POSTN', 'RGCC', 'LCP1', 'RB1', 'COL4A1', 'COL4A2', 'GAS6', 'MMP14', 'CMA1', 'CTSG', 'NID2', 'ERO1A', 'SMOC1', 'PAPLN', 'VIPAS39', 'FLRT2', 'FBLN5', 'GREM1', 'SPINT1', 'WDR72', 'ADAM10', 'MYO1E', 'ANXA2', 'SMAD3', 'THSD4', 'LOXL1', 'ADAMTS7', 'ADAMTSL3', 'FURIN', 'VPS33B', 'ADAMTS17', 'TPSAB1', 'NOXO1', 'MMP25', 'MYH11', 'MMP2', 'MMP15', 'CARMIL2', 'GFOD2', 'SMPD3', 'HAS3', 'ATXN1L', 'ADAMTS18', 'CRISPLD2', 'FOXF1', 'FOXC2', 'ZNF469', 'SERPINF2', 'MFAP4', 'VTN', 'NF1', 'MMP28', 'P3H4', 'FKBP10', 'RAMP2', 'GFAP', 'ITGB3', 'COL1A1', 'SOX9', 'LAMA1', 'IER3IP1', 'SERPINB5', 'ELANE', 'ADAMTS10', 'COL5A3', 'COLGALT1', 'COMP', 'HPN', 'NPHS1', 'APLP1', 'SPINT2', 'MIA', 'TGFB1', 'BCL3', 'ERCC2', 'KLK4', 'KLK7', 'HAS1', 'TCF15', 'FERMT1', 'BMP2', 'CST3', 'MMP24', 'MATN4', 'MMP9', 'SLC2A10', 'SULF2', 'COL9A3', 'APP', 'ADAMTS1', 'ADAMTS5', 'RUNX1', 'COL18A1', 'MMP11', 'TMPRSS6', 'CHADL', 'FBLN1', 'EGFL6', 'GPM6B', 'PRDX4', 'ATP7A', 'NOX1', 'COL4A6', 'COL4A5')

ileum_gene <- rownames(ileum)

ecm_ileum <- intersect(ecm_genes, ileum_gene) #303

#number of genes in ecm_genes adn ecm_ileum are the same

ileum <- AddModuleScore_UCell(ileum, features = list(ecm_genes), name = "ecm_sig")

VlnPlot(ileum, features = "signature_1ecm_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") 

plot1 <- VlnPlot(ileum, features = "signature_1ecm_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") 
png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/figures/VlnPlot_ileum_group_ecm_org_module_11_22_24.png", width = 2500, height = 1800, res = 300)
plot(plot1 + ggtitle("ECM organization") +
       xlab("Cell types") + ylab("Module score"))
dev.off()

#metalopeptidase activity

met_pep_genes <- intersect(ileum_counts$genes, genes_GO0008237$SYMBOL) #170

met_pep_genes <- paste(shQuote(met_pep_genes), collapse=", ")

met_pep_genes <- c('MMP23B', 'MMEL1', 'ECE1', 'ZMPSTE24', 'TRABD2B', 'AGBL4', 'NRDC', 'OMA1', 'MYSM1', 'CLCA2', 'CLCA1', 'CLCA4', 'ADAM15', 'ADAMTS4', 'PAPPA2', 'RNPEP', 'VASH2', 'SPRTN', 'ADAM17', 'AGBL5', 'STAMBP', 'TRABD2A', 'ASTL', 'PSMD14', 'METAP1D', 'ADAM23', 'CPO', 'DNPEP', 'ECEL1', 'RNPEPL1', 'ADAMTS9', 'NUDT16', 'CPB1', 'CPA3', 'MME', 'LMLN', 'CPZ', 'LAP3', 'ADAMTS3', 'METAP1', 'ENPEP', 'CPE', 'TLL1', 'ADAM29', 'ADAMTS16', 'ADAMTS12', 'ADAMTS6', 'NLN', 'ERAP1', 'ERAP2', 'LNPEP', 'LVRN', 'ADAMTS19', 'ADAM19', 'ADAMTS2', 'MEP1A', 'AMZ1', 'ADAM22', 'PMPCB', 'CPA2', 'CPA4', 'CPA5', 'AGBL3', 'PRSS2', 'KEL', 'BMP1', 'ADAM28', 'ADAMDEC1', 'ADAM9', 'ADAM32', 'COPS5', 'CPA6', 'MMP16', 'CPQ', 'EIF3H', 'ERMP1', 'AGTPBP1', 'AOPEP', 'PAPPA', 'ADAMTS13', 'PMPCA', 'PITRM1', 'YME1L1', 'ADAMTS14', 'STAMBPL1', 'IDE', 'TLL2', 'XPNPEP1', 'CPXM2', 'MMP21', 'ADAM12', 'ADAM8', 'EIF3F', 'AGBL2', 'FOLH1', 'NAALADL1', 'DPP3', 'RCE1', 'NAALAD2', 'MMP7', 'MMP20', 'MMP27', 'MMP8', 'MMP10', 'MMP1', 'MMP3', 'MMP12', 'MMP13', 'ADAMTS8', 'ADAMTS15', 'MMP19', 'ATP23', 'CPM', 'TRHDE', 'METAP2', 'LTA4H', 'MMP17', 'MIPEP', 'MMP14', 'ADAM21', 'ADAM20', 'PAPLN', 'VASH1', 'ADAM10', 'ADAMTS7', 'AGBL1', 'ANPEP', 'ADAMTS17', 'MMP25', 'UQCRC2', 'MMP2', 'MMP15', 'DPEP3', 'DPEP2', 'ADAMTS18', 'SPG7', 'DPEP1', 'CHMP1A', 'CPD', 'MMP28', 'ADAM11', 'NPEPPS', 'ACE', 'AMZ2', 'AFG3L2', 'MEP1B', 'CNDP2', 'CNDP1', 'THOP1', 'MPND', 'ADAMTS10', 'PEPD', 'KLK7', 'CPXM1', 'ADAM33', 'MMP24', 'MMP9', 'NPEPL1', 'ADAMTS1', 'ADAMTS5', 'YBEY', 'MMP11', 'TMPRSS6', 'XPNPEP3', 'ACE2', 'MBTPS2', 'PHEX', 'XPNPEP2', 'BRCC3', 'CPN1')

met_pep_ileum <- intersect(met_pep_genes, ileum_gene) #170

ileum <- AddModuleScore_UCell(ileum, features = list(met_pep_genes), name = "met_pep_sig")


VlnPlot(ileum, features = "signature_1met_pep_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") 

plot2 <- VlnPlot(ileum, features = "signature_1met_pep_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") 
png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/figures/VlnPlot_ileum_group_met_pep_module_11_22_24.png", width = 2500, height = 1800, res = 300)
plot(plot2 + ggtitle("Metallopeptidase activity module") +
       xlab("Cell types") + ylab("Module score"))
dev.off()

#chemokine receptor binding 

chem_genes <- intersect(ileum_counts$genes, genes_GO0042379$SYMBOL) #59

chem_genes <- paste(shQuote(chem_genes), collapse=", ")

chem_genes <- c('JAK1', 'S100A14', 'NES', 'XCL2', 'XCL1', 'CNIH4', 'STAT1', 'CCL20', 'CX3CR1', 'CCR2', 'CCRL2', 'CXCL8', 'CXCL6', 'PF4V1', 'CXCL1', 'PF4', 'PPBP', 'CXCL5', 'CXCL3', 'CXCL2', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13', 'CCL28', 'CXCL14', 'CCL26', 'CCL24', 'DEFB1', 'DEFB4A', 'CCL27', 'CCL19', 'CCL21', 'MSMP', 'C5', 'CXCL12', 'CCL22', 'CX3CL1', 'CCL17', 'CKLF', 'CXCL16', 'CCL2', 'CCL7', 'CCL11', 'CCL8', 'CCL13', 'CCL1', 'CCL5', 'CCL16', 'CCL14', 'CCL15', 'CCL23', 'CCL18', 'CCL3', 'CCL4', 'CCL3L1', 'CCL25', 'ITCH', 'TFF2')

chem_genes_ileum <- intersect(chem_genes, ileum_gene) #59

ileum <- AddModuleScore_UCell(ileum, features = list(chem_genes), name = "chem_sig")

VlnPlot(ileum, features = "signature_1chem_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") 

plot3 <- VlnPlot(ileum, features = "signature_1chem_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") 
png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/figures/VlnPlot_ileum_group_chemokine_recept_module_11_22_24.png", width = 2500, height = 1800, res = 300)
plot(plot3 + ggtitle("Chemokine receptor binding module") +
       xlab("Cell types") + ylab("Module score"))
dev.off()

#neutrophil chemotaxis 

neutrophil_chemotaxis_genes <- intersect(ileum_counts$genes, genes_GO0030593$SYMBOL) #102

neutrophil_chemotaxis_genes <- paste(shQuote(neutrophil_chemotaxis_genes), collapse=", ")

neutrophil_chemotaxis_genes <- c('PIK3CD', 'CSF3R', 'EDN2', 'PDE4B', 'VAV3', 'S100A9', 'S100A12', 'S100A8', 'FCER1G', 'XCL2', 'XCL1', 'TGFB2', 'IL1B', 'TNFAIP6', 'DPP4', 'PIKFYVE', 'CXCR2', 'CXCR1', 'CCL20', 'ITGA9', 'BST1', 'SLIT2', 'CXCL8', 'CXCL6', 'PF4V1', 'CXCL1', 'PF4', 'PPBP', 'CXCL5', 'CXCL3', 'CXCL2', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13', 'ITGA1', 'THBS4', 'CD74', 'EDN1', 'RIPOR2', 'TREM1', 'RAC1', 'PPIA', 'CCL26', 'CCL24', 'PIK3CG', 'CCL19', 'CCL21', 'SYK', 'CAMK1D', 'MCU', 'GBF1', 'SAA1', 'MDK', 'JAML', 'TIRAP', 'JAM3', 'C3AR1', 'DNM1L', 'NCKAP1L', 'IL23A', 'PLA2G1B', 'SRP54', 'LGALS3', 'DAPK2', 'PPIB', 'NOD2', 'CCL22', 'CX3CL1', 'CCL17', 'CKLF', 'C1QBP', 'CCL2', 'CCL7', 'CCL11', 'CCL8', 'CCL13', 'CCL1', 'CCL5', 'CCL16', 'CCL14', 'CCL15', 'CCL23', 'CCL18', 'CCL3', 'CCL4', 'CCL3L1', 'CCR7', 'BSG', 'PIP5K1C', 'VAV1', 'CCL25', 'C5AR1', 'C5AR2', 'LBP', 'PREX1', 'EDN3', 'CXADR', 'ITGB2', 'RAC2', 'MOSPD2', 'MPP1')

neut_chem_ileum <- intersect(neutrophil_chemotaxis_genes, ileum_gene) #102

ileum <- AddModuleScore_UCell(ileum, features = list(neutrophil_chemotaxis_genes), name = "neut_chemo_sig")

VlnPlot(ileum, features = "signature_1neut_chemo_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") 

plot4 <- VlnPlot(ileum, features = "signature_1neut_chemo_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") 
png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/figures/VlnPlot_ileum_group_neutrophil_chemo_module_11_22_24.png", width = 2500, height = 1800, res = 300)
plot(plot4 + ggtitle("Neutrophil chemotaxis module") +
       xlab("Cell types") + ylab("Module score"))
dev.off()

#------------subset group 1 and group 2------------#

Idents(ileum) <- ileum$group

ileum_g1 <- subset(ileum, idents = "group1")

#save the seurat object 
saveRDS(ileum_g1, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_group1_11_22_24.rds")

ileum_g2 <- subset(ileum, idents = "group2")

#save the seurat object 
saveRDS(ileum_g2, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_group2_11_22_24.rds")

#----------------Hierarchical clustering--------------#

#11/25/2024

#use cell proportions 

#load in ileum data 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_11_19_24.rds")


ileum@meta.data$cell_typev2 <- ileum@meta.data$seurat_clusters
ileum$cell_typev2 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21'),
  to = c("Naive B Cell", 'Stem/Paneth Cell', 'Enterocyte', 'CD4 T Cell', 'NK/CD8 T Cell', 'Enterocyte', 'Macrophage', "Plasma Cell", "Plasma Cell", 'Monocyte', 'Enterocyte', 'Goblet Cell', 'Goblet Cell', 'Ambig. T Cell', 'Mesenchymal Cell', "Pro-Inflammatory Monocyte", "TA Cell", "Memory B Cell", "Endothelial Cell", 'Mast Cell', "Tuft Cell", "Cycling Plasma Cell")
)

ileum$cell_typev2 <- factor(ileum$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "NK/CD8 T Cell", "Ambig. T Cell", "Macrophage", "Monocyte", "Pro-Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Enterocyte", "Goblet Cell", "TA Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

metadata <- ileum@meta.data

df_unique <- metadata %>% distinct(orig.ident, .keep_all = TRUE)

library(speckle)

props <- getTransformedProps(ileum$cell_typev2, ileum$donor, transform="logit")

props_transformed <- props$TransformedProps

library(pheatmap)
library(factoextra)
props_mat <- as.matrix(props_transformed)
distance_matrix <- dist(props_mat, method = "euclidean")
hclust_object <- hclust(distance_matrix, method = "ward.D2")  # Perform Ward's hierarchical clustering

props_mat <- as.matrix(props_transformed)
props_df <- as.data.frame(props_transformed)
# Create the total within-cluster sum of squares (WSS) plot
#fviz_nbclust(as.data.frame(props_transformed), FUNcluster = hcut, method = "wss", k.max = 10)  # Adjust k.max based on your data

elbow_value <- fviz_nbclust(as.data.frame(props_transformed), FUN = hcut, method = "wss", k.max = 15, print.summary = T)

elbow_value 

#optimal number of clusters = 3
#save the heatmap 

save_pheatmap_pdf <- function(x, filename, width=5, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}



hmap <- pheatmap(props_transformed, 
                 clustering_method = "ward.D2", # Using Ward's method for hierarchical clustering
                 scale = "row",
                 cluster_rows = FALSE,
                 cluster_cols = TRUE, 
                 #clustering_distance_rows = "euclidean", 
                 clustering_distance_columns = "euclidean",
                 cutree_cols = 3,
                 color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
hmap


save_pheatmap_pdf(hmap, "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/ileum/heatmap_batch1_13_ileum_hierarchical_cluster_cellprop_v2_11_25_24.pdf", width = 8, height = 5)

#save the ileum seurat object 
saveRDS(ileum, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_udpate_11_26_24.rds")
#---------prep data for IF status comparison dreamlet----------#

#11/26/2024

#do this on server 
#load in seurat object 
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_udpate_11_26_24.rds")
DefaultAssay(ileum) <- "RNA"
ileum.sce <- as.SingleCellExperiment(ileum)
saveRDS(ileum.sce, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_scITD_group_sce_11_20_24.rds")



#--------test making heatmap of genes in ecm/chemokine pathways--------#
#pseudobulk data 

#first get average count of each gene per donor - do this on the server 
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_scITD_groupInfo_scale_11_26_24.rds")

DefaultAssay(ileum) <- "RNA"

ecm_genes <- c('VWA1', 'MMP23B', 'ANGPTL7', 'TNFRSF1B', 'PDPN', 'MATN1', 'COL16A1', 'COL8A2', 'COL9A2', 'TIE1', 'CCN1', 'COL24A1', 'COL11A1', 'NTNG1', 'ADAMTSL4', 'CTSS', 'CTSK', 'PBXIP1', 'ADAM15', 'HAPLN2', 'ADAMTS4', 'OLFML2B', 'DDR2', 'DPT', 'TNR', 'QSOX1', 'LAMC1', 'COLGALT2', 'HMCN1', 'ELF3', 'FMOD', 'LAMB3', 'TGFB2', 'WNT3A', 'AGT', 'EXOC8', 'NID1', 'PXDN', 'MATN3', 'EMILIN1', 'MPV17', 'VIT', 'CYP1B1', 'ANTXR1', 'LOXL3', 'CLASP1', 'DPP4', 'FAP', 'COL3A1', 'COL5A2', 'CFLAR', 'IHH', 'COL4A4', 'COL4A3', 'FBLN2', 'COLQ', 'CRTAP', 'CLASP2', 'KIF9', 'LAMB2', 'DAG1', 'ADAMTS9', 'COL8A1', 'ABI3BP', 'IMPG2', 'PHLDB2', 'CCDC80', 'COL6A5', 'COL6A6', 'NPHP3', 'PTX3', 'MELTF', 'PDGFRA', 'ADAMTS3', 'SLC39A8', 'NPNT', 'PRDM5', 'NDNF', 'SFRP2', 'RXFP1', 'TLL1', 'ADAMTS16', 'ADAMTS12', 'EGFLAM', 'ADAMTS6', 'LOX', 'ADAMTS19', 'TGFBI', 'SPINK5', 'SH3PXD2B', 'FGFR4', 'COL23A1', 'ADAMTS2', 'FOXF2', 'FOXC1', 'ADTRP', 'FLOT1', 'DDR1', 'TNF', 'TNXB', 'COL11A2', 'COL19A1', 'COL9A1', 'COL12A1', 'IMPG1', 'COL10A1', 'LAMA2', 'CCN2', 'SERAC1', 'PLG', 'SMOC2', 'FSCN1', 'COL28A1', 'IL6', 'AEBP1', 'ELN', 'COL1A2', 'PLOD3', 'LAMB1', 'LAMB4', 'CAV2', 'CAV1', 'ST7', 'PRSS1', 'PRSS2', 'DNAJB6', 'BMP1', 'LOXL2', 'SULF1', 'MMP16', 'MATN2', 'EXT1', 'TNFRSF11B', 'COL14A1', 'HAS2', 'COL22A1', 'SCX', 'WASHC1', 'RIC1', 'ADAMTSL1', 'B4GALT1', 'RECK', 'ECM2', 'CTSV', 'COL15A1', 'TGFBR1', 'TMEM38B', 'COL27A1', 'OLFML2A', 'ENG', 'ABL1', 'POMT1', 'NTNG2', 'ADAMTS13', 'ADAMTSL2', 'COL5A1', 'NOTCH1', 'ITGA8', 'ITGB1', 'COL13A1', 'ADAMTS14', 'SPOCK2', 'P4HA1', 'TLL2', 'LOXL4', 'HPSE2', 'KAZALD1', 'NFKB2', 'COL17A1', 'MMP21', 'ADAM8', 'RIC8A', 'GAS2', 'WT1', 'HSD17B12', 'CREB3L1', 'LTBP3', 'EFEMP2', 'SERPINH1', 'MMP7', 'MMP20', 'MMP27', 'MMP8', 'MMP10', 'MMP1', 'MMP3', 'MMP12', 'MMP13', 'MPZL3', 'PHLDB1', 'ADAMTS8', 'ADAMTS15', 'TNFRSF1A', 'COL2A1', 'MMP19', 'LRP1', 'LUM', 'NTN4', 'MMP17', 'POSTN', 'RGCC', 'LCP1', 'RB1', 'COL4A1', 'COL4A2', 'GAS6', 'MMP14', 'CMA1', 'CTSG', 'NID2', 'ERO1A', 'SMOC1', 'PAPLN', 'VIPAS39', 'FLRT2', 'FBLN5', 'GREM1', 'SPINT1', 'WDR72', 'ADAM10', 'MYO1E', 'ANXA2', 'SMAD3', 'THSD4', 'LOXL1', 'ADAMTS7', 'ADAMTSL3', 'FURIN', 'VPS33B', 'ADAMTS17', 'TPSAB1', 'NOXO1', 'MMP25', 'MYH11', 'MMP2', 'MMP15', 'CARMIL2', 'GFOD2', 'SMPD3', 'HAS3', 'ATXN1L', 'ADAMTS18', 'CRISPLD2', 'FOXF1', 'FOXC2', 'ZNF469', 'SERPINF2', 'MFAP4', 'VTN', 'NF1', 'MMP28', 'P3H4', 'FKBP10', 'RAMP2', 'GFAP', 'ITGB3', 'COL1A1', 'SOX9', 'LAMA1', 'IER3IP1', 'SERPINB5', 'ELANE', 'ADAMTS10', 'COL5A3', 'COLGALT1', 'COMP', 'HPN', 'NPHS1', 'APLP1', 'SPINT2', 'MIA', 'TGFB1', 'BCL3', 'ERCC2', 'KLK4', 'KLK7', 'HAS1', 'TCF15', 'FERMT1', 'BMP2', 'CST3', 'MMP24', 'MATN4', 'MMP9', 'SLC2A10', 'SULF2', 'COL9A3', 'APP', 'ADAMTS1', 'ADAMTS5', 'RUNX1', 'COL18A1', 'MMP11', 'TMPRSS6', 'CHADL', 'FBLN1', 'EGFL6', 'GPM6B', 'PRDX4', 'ATP7A', 'NOX1', 'COL4A6', 'COL4A5')

ileum_counts <- AverageExpression(ileum, assays = "RNA", features = ecm_genes, group.by = "donor", layer = "data")

write.csv(ileum_counts, file = "/storage/home/swashburn30/Helmsley/Batch1_13/DEG_analysis/helm_batch1_13_ileum_scITD_group_donor_avg_counts_11_26_24.csv")

#-------read in average counts per donor-------#

ileum_donor <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_scITD_group_donor_avg_counts_11_26_24.csv")

rownames(ileum_donor) <- ileum_donor$X

colnames(ileum_donor) <- c("genes", "donor1", "donor2", "donor3", "donor4", "donor5", "donor6", "donor7", "donor8", "donor9", "donor10", "donor11", "donor12", "donor14", "donor17", "donor18", "donor20", "donor21", "donor23", "donor24", "donor26", "donor27", "donor30", "donor34")

ileum_donor <- subset(ileum_donor, select = -genes)

# Example donor metadata
donor_metadata <- data.frame(
  donor = colnames(ileum_donor),
  group = c("group1", "group1", "group1", "group1", "group1", "group2", "group2", "group1", "group1", "group2", "group1", "group2", "group1", "group2", "group1", "group1", "group2", "group1", "group2", "group2", "group1", "group2", "group2") # Example groups
)

# Annotate donors
donor_annot <- HeatmapAnnotation(Group = donor_metadata$group)

#subset to 17 genes and try again; code above worked but hard to see each gene 
#scaled average gene expression 

subset_ileum_donor <- ileum_donor[rownames(ileum_donor) %in% group_gene, ]

mat_scaled_sub <- t(scale(t(subset_ileum_donor)))

#this worked - use this 
Heatmap(mat_scaled_sub, top_annotation = donor_annot, show_heatmap_legend = TRUE, cluster_columns = TRUE, name = "Z-score")

# Create heatmap
pheatmap(
  as.matrix(subset_ileum_donor),
  cluster_rows = TRUE,           # Cluster genes
  cluster_cols = TRUE,           # Cluster cells
  scale = "row",                 # Scale expression per gene
  annotation_col = donor_annot, # Annotate by group
  color = colorRampPalette(c("blue", "white", "red"))(50) # Color scale
)

#---------Heatmap for neutrophil chemotaxis genes-------#

#do this on server 

#first get average count of each gene per donor - do this on the server 
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_scITD_groupInfo_scale_11_26_24.rds")

DefaultAssay(ileum) <- "RNA"

neutrophil_chemotaxis_genes <- c('PIK3CD', 'CSF3R', 'EDN2', 'PDE4B', 'VAV3', 'S100A9', 'S100A12', 'S100A8', 'FCER1G', 'XCL2', 'XCL1', 'TGFB2', 'IL1B', 'TNFAIP6', 'DPP4', 'PIKFYVE', 'CXCR2', 'CXCR1', 'CCL20', 'ITGA9', 'BST1', 'SLIT2', 'CXCL8', 'CXCL6', 'PF4V1', 'CXCL1', 'PF4', 'PPBP', 'CXCL5', 'CXCL3', 'CXCL2', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13', 'ITGA1', 'THBS4', 'CD74', 'EDN1', 'RIPOR2', 'TREM1', 'RAC1', 'PPIA', 'CCL26', 'CCL24', 'PIK3CG', 'CCL19', 'CCL21', 'SYK', 'CAMK1D', 'MCU', 'GBF1', 'SAA1', 'MDK', 'JAML', 'TIRAP', 'JAM3', 'C3AR1', 'DNM1L', 'NCKAP1L', 'IL23A', 'PLA2G1B', 'SRP54', 'LGALS3', 'DAPK2', 'PPIB', 'NOD2', 'CCL22', 'CX3CL1', 'CCL17', 'CKLF', 'C1QBP', 'CCL2', 'CCL7', 'CCL11', 'CCL8', 'CCL13', 'CCL1', 'CCL5', 'CCL16', 'CCL14', 'CCL15', 'CCL23', 'CCL18', 'CCL3', 'CCL4', 'CCL3L1', 'CCR7', 'BSG', 'PIP5K1C', 'VAV1', 'CCL25', 'C5AR1', 'C5AR2', 'LBP', 'PREX1', 'EDN3', 'CXADR', 'ITGB2', 'RAC2', 'MOSPD2', 'MPP1')

ileum_counts <- AverageExpression(ileum, assays = "RNA", features = neutrophil_chemotaxis_genes, group.by = "donor", layer = "data")

write.csv(ileum_counts, file = "/storage/home/swashburn30/Helmsley/Batch1_13/DEG_analysis/helm_batch1_13_ileum_scITD_group_donor_avg_counts_neut_chem_11_27_24.csv")

#------read in the average gene expression data-------#
ileum_neut_chem <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_scITD_group_donor_avg_counts_neut_chem_11_27_24.csv")

rownames(ileum_neut_chem) <- ileum_neut_chem$X

colnames(ileum_neut_chem) <- c("genes", "donor1", "donor2", "donor3", "donor4", "donor5", "donor6", "donor7", "donor8", "donor9", "donor10", "donor11", "donor12", "donor14", "donor17", "donor18", "donor20", "donor21", "donor23", "donor24", "donor26", "donor27", "donor30", "donor34")

ileum_neut_chem <- subset(ileum_neut_chem, select = -genes)

#top neutrophil chemokine genes 
ileum_sub_neut <- ileum_counts[ileum_counts$genes %in% neutrophil_chemotaxis_genes, ]

ileum_g1_neut <- subset(ileum_sub_neut, sig_strict == "1")

ileum_g2_neut <- subset(ileum_sub_neut, sig_strict == "-1")

g1_gene_neut <- rownames(ileum_g1_neut)
g2_gene_neut <- rownames(ileum_g2_neut)
group_gene_neut <- c(g1_gene_neut, g2_gene_neut)

# Example donor metadata
donor_metadata <- data.frame(
  donor = colnames(ileum_donor),
  group = c("group1", "group1", "group1", "group1", "group1", "group2", "group2", "group1", "group1", "group2", "group1", "group2", "group1", "group2", "group1", "group1", "group2", "group1", "group2", "group2", "group1", "group2", "group2") # Example groups
)

subset_ileum_neut <- ileum_neut_chem[rownames(ileum_neut_chem) %in% group_gene_neut, ]

mat_scaled_sub_neut <- t(scale(t(subset_ileum_neut)))

#this worked 
Heatmap(mat_scaled_sub_neut, top_annotation = donor_annot, show_heatmap_legend = TRUE, cluster_columns = TRUE, name = "Z-score")



#---------Heatmap for chemokine receptor binding module genes-------#

#do this on server 

#first get average count of each gene per donor - do this on the server 
ileum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_scITD_groupInfo_scale_11_26_24.rds")

DefaultAssay(ileum) <- "RNA"

chem_genes <- c('JAK1', 'S100A14', 'NES', 'XCL2', 'XCL1', 'CNIH4', 'STAT1', 'CCL20', 'CX3CR1', 'CCR2', 'CCRL2', 'CXCL8', 'CXCL6', 'PF4V1', 'CXCL1', 'PF4', 'PPBP', 'CXCL5', 'CXCL3', 'CXCL2', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13', 'CCL28', 'CXCL14', 'CCL26', 'CCL24', 'DEFB1', 'DEFB4A', 'CCL27', 'CCL19', 'CCL21', 'MSMP', 'C5', 'CXCL12', 'CCL22', 'CX3CL1', 'CCL17', 'CKLF', 'CXCL16', 'CCL2', 'CCL7', 'CCL11', 'CCL8', 'CCL13', 'CCL1', 'CCL5', 'CCL16', 'CCL14', 'CCL15', 'CCL23', 'CCL18', 'CCL3', 'CCL4', 'CCL3L1', 'CCL25', 'ITCH', 'TFF2')

ileum_counts <- AverageExpression(ileum, assays = "RNA", features = chem_genes, group.by = "donor", layer = "data")

write.csv(ileum_counts, file = "/storage/home/swashburn30/Helmsley/Batch1_13/DEG_analysis/helm_batch1_13_ileum_scITD_group_donor_avg_counts_chemokine_receptor_12_18_24.csv")

#------read in the average gene expression data-------#
ileum_chem <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_scITD_group_donor_avg_counts_chemokine_receptor_12_18_24.csv")

rownames(ileum_chem) <- ileum_chem$X

colnames(ileum_chem) <- c("genes", "donor1", "donor2", "donor3", "donor4", "donor5", "donor6", "donor7", "donor8", "donor9", "donor10", "donor11", "donor12", "donor14", "donor17", "donor18", "donor20", "donor21", "donor23", "donor24", "donor26", "donor27", "donor30", "donor34")

ileum_chem <- subset(ileum_chem, select = -genes)

#top neutrophil chemokine genes 
ileum_sub_chem <- ileum_counts[ileum_counts$genes %in% chem_genes, ]

ileum_g1_chem <- subset(ileum_sub_chem, sig_strict == "1")

ileum_g2_chem <- subset(ileum_sub_chem, sig_strict == "-1")

g1_gene_chem <- rownames(ileum_g1_chem)
g2_gene_chem <- rownames(ileum_g2_chem)
group_gene_chem <- c(g1_gene_chem, g2_gene_chem)

# Example donor metadata
donor_metadata <- data.frame(
  donor = colnames(ileum_chem),
  group = c("group1", "group1", "group1", "group1", "group1", "group2", "group2", "group1", "group1", "group2", "group1", "group2", "group1", "group2", "group1", "group1", "group2", "group1", "group2", "group2", "group1", "group2", "group2") # Example groups
)

subset_ileum_chem <- ileum_chem[rownames(ileum_chem) %in% group_gene_chem, ]

mat_scaled_sub_chem <- t(scale(t(subset_ileum_chem)))

Heatmap(mat_scaled_sub_chem, top_annotation = donor_annot, show_heatmap_legend = TRUE, cluster_columns = TRUE, name = "Z-score")


#------------1/3/2025------------#

library(speckle)

#test if monocyte proportion is statistically different between group 1 and group 2 donors 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_scITD_groupInfo_11_20_24.rds")

ileum$donor <- factor(ileum$donor, levels = c("donor1", "donor2", "donor3", "donor4", "donor5", "donor6", "donor7", "donor8", "donor9", "donor10", "donor11", "donor12", "donor14", "donor17", "donor18", "donor20", "donor21", "donor23", "donor24", "donor26", "donor27", "donor30", "donor34"))

ileum$group <- factor(ileum$group, levels = c("group1", "group2"))

group <- propeller(clusters = ileum$cell_typev2, sample = ileum$donor, group = ileum$group)
