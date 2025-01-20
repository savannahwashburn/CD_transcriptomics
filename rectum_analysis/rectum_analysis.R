#Helmsley Batch 1-13 RECTUM Analysis


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

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_clustered_10_17_24.rds")

rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
rectum$batch <- factor(rectum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(rectum, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(rectum, features = 'nCount_RNA', pt.size = 0)
VlnPlot(rectum, features = 'percent.mt', pt.size = 0)

#visualize distribution of data (sample)
VlnPlot(rectum, features = 'nFeature_RNA', group.by = "orig.ident", cols = c("helm_sam3" = "#BBCCEE", "helm_sam9" = "#BBCCEE", "helm_sam25" = "#99DDFF", "helm_sam30" = "#99DDFF", "helm_sam34" = "#99DDFF", "helm_sam37" = "#99DDFF", "helm_sam46" = "#44BB99", "helm_sam49" = "#44BB99", "helm_sam51" = "#EE8866", "helm_sam55" = "#EE8866", "helm_sam58" = "#EE8866", "helm_sam61" = "#EE8866", "helm_sam64" = "#BBCC33", "helm_sam67" = "#BBCC33", "helm_sam73" = "#BBCC33", "helm_sam79" = "#882255", "helm_sam85" = "#882255", "helm_sam87" = "#882255", "helm_sam90" = "#EECC66", "helm_sam93" = "#EECC66", "helm_sam97" = "#EECC66", "helm_sam103" = "#EECC66", "helm_sam111" = "#FE036A", "helm_sam114" = "#FE036A", "helm_sam115" = "#FE036A", "helm_sam121" = "#A8CED2", "helm_sam130" = "#A8CED2", "helm_sam136" = "#B2DE7C", "helm_sam149" = "#B2DE7C", "helm_sam155" = "#E19390", "helm_sam158" = "#E19390", "helm_sam164" = "#E19390", "helm_sam173" = "#BA6BD4"), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
VlnPlot(rectum, features = 'nCount_RNA', group.by = "orig.ident", cols = c("helm_sam3" = "#BBCCEE", "helm_sam9" = "#BBCCEE", "helm_sam25" = "#99DDFF", "helm_sam30" = "#99DDFF", "helm_sam34" = "#99DDFF", "helm_sam37" = "#99DDFF", "helm_sam46" = "#44BB99", "helm_sam49" = "#44BB99", "helm_sam51" = "#EE8866", "helm_sam55" = "#EE8866", "helm_sam58" = "#EE8866", "helm_sam61" = "#EE8866", "helm_sam64" = "#BBCC33", "helm_sam67" = "#BBCC33", "helm_sam73" = "#BBCC33", "helm_sam79" = "#882255", "helm_sam85" = "#882255", "helm_sam87" = "#882255", "helm_sam90" = "#EECC66", "helm_sam93" = "#EECC66", "helm_sam97" = "#EECC66", "helm_sam103" = "#EECC66", "helm_sam111" = "#FE036A", "helm_sam114" = "#FE036A", "helm_sam115" = "#FE036A", "helm_sam121" = "#A8CED2", "helm_sam130" = "#A8CED2", "helm_sam136" = "#B2DE7C", "helm_sam149" = "#B2DE7C", "helm_sam155" = "#E19390", "helm_sam158" = "#E19390", "helm_sam164" = "#E19390", "helm_sam173" = "#BA6BD4"), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
VlnPlot(rectum, features = 'percent.mt', group.by = "orig.ident", cols = c("helm_sam3" = "#BBCCEE", "helm_sam9" = "#BBCCEE", "helm_sam25" = "#99DDFF", "helm_sam30" = "#99DDFF", "helm_sam34" = "#99DDFF", "helm_sam37" = "#99DDFF", "helm_sam46" = "#44BB99", "helm_sam49" = "#44BB99", "helm_sam51" = "#EE8866", "helm_sam55" = "#EE8866", "helm_sam58" = "#EE8866", "helm_sam61" = "#EE8866", "helm_sam64" = "#BBCC33", "helm_sam67" = "#BBCC33", "helm_sam73" = "#BBCC33", "helm_sam79" = "#882255", "helm_sam85" = "#882255", "helm_sam87" = "#882255", "helm_sam90" = "#EECC66", "helm_sam93" = "#EECC66", "helm_sam97" = "#EECC66", "helm_sam103" = "#EECC66", "helm_sam111" = "#FE036A", "helm_sam114" = "#FE036A", "helm_sam115" = "#FE036A", "helm_sam121" = "#A8CED2", "helm_sam130" = "#A8CED2", "helm_sam136" = "#B2DE7C", "helm_sam149" = "#B2DE7C", "helm_sam155" = "#E19390", "helm_sam158" = "#E19390", "helm_sam164" = "#E19390", "helm_sam173" = "#BA6BD4"), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#visualize PCA results
ElbowPlot(rectum)
DimPlot(rectum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(rectum, reduction = 'pca', group.by = "orig.ident")
DimPlot(rectum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(rectum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(rectum, reduction = 'umap', group.by = "orig.ident")
DimPlot(rectum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
Idents(rectum) <- rectum$seurat_clusters
pt1 <- table(Idents(rectum), rectum$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", '14'))

library(randomcoloR)
no_of_colors <- 17

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
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
Idents(rectum) <- rectum$batch
pt2 <- table(Idents(rectum), rectum$seurat_clusters)
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

#------------preliminary rPCA clustering results-------------#

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_rpca_10_17_24.rds")

rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
rectum$batch <- factor(rectum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(rectum, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(rectum, features = 'nCount_RNA', pt.size = 0)
VlnPlot(rectum, features = 'percent.mt', pt.size = 0)

#visualize PCA results
ElbowPlot(rectum)
DimPlot(rectum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(rectum, reduction = 'pca', group.by = "orig.ident")
DimPlot(rectum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(rectum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(rectum, reduction = 'umap', group.by = "orig.ident")
DimPlot(rectum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
Idents(rectum) <- rectum$seurat_clusters
pt1 <- table(Idents(rectum), rectum$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", '14', '15'))

library(randomcoloR)
no_of_colors <- 17

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
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
Idents(rectum) <- rectum$batch
pt2 <- table(Idents(rectum), rectum$seurat_clusters)
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


#-----------preliminary harmony clustering results-----------#

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_harmony_10_17_24.rds")

rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
rectum$batch <- factor(rectum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(rectum, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(rectum, features = 'nCount_RNA', pt.size = 0)
VlnPlot(rectum, features = 'percent.mt', pt.size = 0)

#visualize PCA results
ElbowPlot(rectum)
DimPlot(rectum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(rectum, reduction = 'pca', group.by = "orig.ident")
DimPlot(rectum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(rectum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(rectum, reduction = 'umap', group.by = "orig.ident")
DimPlot(rectum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
Idents(rectum) <- rectum$seurat_clusters
pt1 <- table(Idents(rectum), rectum$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

library(randomcoloR)
no_of_colors <- 17

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
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
Idents(rectum) <- rectum$batch
pt2 <- table(Idents(rectum), rectum$seurat_clusters)
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

#-----------preliminary CCA clustering results----------#

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_cca_10_17_24.rds")

rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
rectum$batch <- factor(rectum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(rectum, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(rectum, features = 'nCount_RNA', pt.size = 0)
VlnPlot(rectum, features = 'percent.mt', pt.size = 0)

#visualize PCA results
ElbowPlot(rectum)
DimPlot(rectum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(rectum, reduction = 'pca', group.by = "orig.ident")
DimPlot(rectum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(rectum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(rectum, reduction = 'umap', group.by = "orig.ident")
DimPlot(rectum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
Idents(rectum) <- rectum$seurat_clusters
pt1 <- table(Idents(rectum), rectum$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", '14'))

library(randomcoloR)
no_of_colors <- 17

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
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
Idents(rectum) <- rectum$batch
pt2 <- table(Idents(rectum), rectum$seurat_clusters)
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


#------------Low res stability assessment-------------#


#-----Query 1------#

cluster1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub_res20_1_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #13031

#cluster 1
Idents(cluster1) <- cluster1$seurat_clusters
sub1 <- WhichCells(cluster1, idents = c('1'))

Idents(cluster1) <- cluster1$predicted.id
pred1 <- WhichCells(cluster1, idents = '1')

int1 <- intersect(sub1, pred1) #6548

#cluster 2
Idents(cluster1) <- cluster1$seurat_clusters
sub2 <- WhichCells(cluster1, idents = c('2', '3'))

Idents(cluster1) <- cluster1$predicted.id
pred2 <- WhichCells(cluster1, idents = '2')

int2 <- intersect(sub2, pred2) #6053

#cluster 3
Idents(cluster1) <- cluster1$seurat_clusters
sub3 <- WhichCells(cluster1, idents = c('2', '5'))

Idents(cluster1) <- cluster1$predicted.id
pred3 <- WhichCells(cluster1, idents = '3')

int3 <- intersect(sub3, pred3) #5760

#cluster 4
Idents(cluster1) <- cluster1$seurat_clusters
sub4 <- WhichCells(cluster1, idents = c('4'))

Idents(cluster1) <- cluster1$predicted.id
pred4 <- WhichCells(cluster1, idents = '4')

int4 <- intersect(sub4, pred4) #3257

#cluster 5
Idents(cluster1) <- cluster1$seurat_clusters
sub5 <- WhichCells(cluster1, idents = c('6'))

Idents(cluster1) <- cluster1$predicted.id
pred5 <- WhichCells(cluster1, idents = '5')

int5 <- intersect(sub5, pred5) #2234

#cluster 6
Idents(cluster1) <- cluster1$seurat_clusters
sub6 <- WhichCells(cluster1, idents = c('7'))

Idents(cluster1) <- cluster1$predicted.id
pred6 <- WhichCells(cluster1, idents = '6')

int6 <- intersect(sub6, pred6) #2050

#cluster 7
Idents(cluster1) <- cluster1$seurat_clusters
sub7 <- WhichCells(cluster1, idents = c('8'))

Idents(cluster1) <- cluster1$predicted.id
pred7 <- WhichCells(cluster1, idents = '7')

int7 <- intersect(sub7, pred7) #1075

#cluster 8
Idents(cluster1) <- cluster1$seurat_clusters
sub8 <- WhichCells(cluster1, idents = c('9'))

Idents(cluster1) <- cluster1$predicted.id
pred8 <- WhichCells(cluster1, idents = '8')

int8 <- intersect(sub8, pred8) #1010

#cluster 9 
Idents(cluster1) <- cluster1$seurat_clusters
sub9 <- WhichCells(cluster1, idents = c('13'))

Idents(cluster1) <- cluster1$predicted.id
pred9 <- WhichCells(cluster1, idents = '9')

int9 <- intersect(sub9, pred9) #278

#cluster 10 
Idents(cluster1) <- cluster1$seurat_clusters
sub10 <- WhichCells(cluster1, idents = c('11'))

Idents(cluster1) <- cluster1$predicted.id
pred10 <- WhichCells(cluster1, idents = '10')

int10 <- intersect(sub10, pred10) #364

#cluster 11
Idents(cluster1) <- cluster1$seurat_clusters
sub11 <- WhichCells(cluster1, idents = c('10'))

Idents(cluster1) <- cluster1$predicted.id
pred11 <- WhichCells(cluster1, idents = '11')

int11 <- intersect(sub11, pred11) #320

#cluster 12
Idents(cluster1) <- cluster1$seurat_clusters
sub12 <- WhichCells(cluster1, idents = c('12'))

Idents(cluster1) <- cluster1$predicted.id
pred12 <- WhichCells(cluster1, idents = '12')

int12 <- intersect(sub12, pred12) #291

#cluster 13
Idents(cluster1) <- cluster1$seurat_clusters
sub13 <- WhichCells(cluster1, idents = c('10'))

Idents(cluster1) <- cluster1$predicted.id
pred13 <- WhichCells(cluster1, idents = '13')

int13 <- intersect(sub13, pred13) #162

#cluster 14 
Idents(cluster1) <- cluster1$seurat_clusters
sub14 <- WhichCells(cluster1, idents = c('10'))

Idents(cluster1) <- cluster1$predicted.id
pred14 <- WhichCells(cluster1, idents = '14')

int14 <- intersect(sub14, pred14) #364

#cluster 15
Idents(cluster1) <- cluster1$seurat_clusters
sub15 <- WhichCells(cluster1, idents = c('4'))

Idents(cluster1) <- cluster1$predicted.id
pred15 <- WhichCells(cluster1, idents = '15')

int15 <- intersect(sub15, pred15) #98

#total 
#pre:44178
#42674
total1 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total1,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query1_low_res_barcodes_11_2_24.csv', row.names=F)

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


#-----Query 2------#

cluster2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub2_res20_1_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #13247

#cluster 1
Idents(cluster2) <- cluster2$seurat_clusters
sub1 <- WhichCells(cluster2, idents = c('1'))

Idents(cluster2) <- cluster2$predicted.id
pred1 <- WhichCells(cluster2, idents = '1')

int1 <- intersect(sub1, pred1) #6533

#cluster 2
Idents(cluster2) <- cluster2$seurat_clusters
sub2 <- WhichCells(cluster2, idents = c('3'))

Idents(cluster2) <- cluster2$predicted.id
pred2 <- WhichCells(cluster2, idents = '2')

int2 <- intersect(sub2, pred2) #5449

#cluster 3 
Idents(cluster2) <- cluster2$seurat_clusters
sub3 <- WhichCells(cluster2, idents = c('2'))

Idents(cluster2) <- cluster2$predicted.id
pred3 <- WhichCells(cluster2, idents = '3')

int3 <- intersect(sub3, pred3) #5078

#cluster 4
Idents(cluster2) <- cluster2$seurat_clusters
sub4 <- WhichCells(cluster2, idents = c('4'))

Idents(cluster2) <- cluster2$predicted.id
pred4 <- WhichCells(cluster2, idents = '4')

int4 <- intersect(sub4, pred4) #3280

#cluster 5
Idents(cluster2) <- cluster2$seurat_clusters
sub5 <- WhichCells(cluster2, idents = c('5'))

Idents(cluster2) <- cluster2$predicted.id
pred5 <- WhichCells(cluster2, idents = '5')

int5 <- intersect(sub5, pred5) #2467

#cluster 6
Idents(cluster2) <- cluster2$seurat_clusters
sub6 <- WhichCells(cluster2, idents = c('6'))

Idents(cluster2) <- cluster2$predicted.id
pred6 <- WhichCells(cluster2, idents = '6')

int6 <- intersect(sub6, pred6) #2067

#cluster 7 
Idents(cluster2) <- cluster2$seurat_clusters
sub7 <- WhichCells(cluster2, idents = c('7'))

Idents(cluster2) <- cluster2$predicted.id
pred7 <- WhichCells(cluster2, idents = '7')

int7 <- intersect(sub7, pred7) #1067

#cluster 8 
Idents(cluster2) <- cluster2$seurat_clusters
sub8 <- WhichCells(cluster2, idents = c('8'))

Idents(cluster2) <- cluster2$predicted.id
pred8 <- WhichCells(cluster2, idents = '8')

int8 <- intersect(sub8, pred8) #945

#cluster 9 
Idents(cluster2) <- cluster2$seurat_clusters
sub9 <- WhichCells(cluster2, idents = c('9'))

Idents(cluster2) <- cluster2$predicted.id
pred9 <- WhichCells(cluster2, idents = '9')

int9 <- intersect(sub9, pred9) #395

#cluster 10 
Idents(cluster2) <- cluster2$seurat_clusters
sub10 <- WhichCells(cluster2, idents = c('11'))

Idents(cluster2) <- cluster2$predicted.id
pred10 <- WhichCells(cluster2, idents = '10')

int10 <- intersect(sub10, pred10) #372

#cluster 11
Idents(cluster2) <- cluster2$seurat_clusters
sub11 <- WhichCells(cluster2, idents = c('10'))

Idents(cluster2) <- cluster2$predicted.id
pred11 <- WhichCells(cluster2, idents = '11')

int11 <- intersect(sub11, pred11) #301

#cluster 12
Idents(cluster2) <- cluster2$seurat_clusters
sub12 <- WhichCells(cluster2, idents = c('12'))

Idents(cluster2) <- cluster2$predicted.id
pred12 <- WhichCells(cluster2, idents = '12')

int12 <- intersect(sub12, pred12) #295

#cluster 13
Idents(cluster2) <- cluster2$seurat_clusters
sub13 <- WhichCells(cluster2, idents = c('7'))

Idents(cluster2) <- cluster2$predicted.id
pred13 <- WhichCells(cluster2, idents = '13')

int13 <- intersect(sub13, pred13) #145

#cluster 14
Idents(cluster2) <- cluster2$seurat_clusters
sub14 <- WhichCells(cluster2, idents = c('10'))

Idents(cluster2) <- cluster2$predicted.id
pred14 <- WhichCells(cluster2, idents = '14')

int14 <- intersect(sub14, pred14) #130

#cluster 15
Idents(cluster2) <- cluster2$seurat_clusters
sub15 <- WhichCells(cluster2, idents = c('4'))

Idents(cluster2) <- cluster2$predicted.id
pred15 <- WhichCells(cluster2, idents = '15')

int15 <- intersect(sub15, pred15) #82

#total 
#pre:44177
#41817
total2 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total2,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query2_low_res_barcodes_11_2_24.csv', row.names=F)

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


#-----Query 3------#

cluster3 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub_res20_2_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #12924

#cluster 1
Idents(cluster3) <- cluster3$seurat_clusters
sub1 <- WhichCells(cluster3, idents = c('2'))

Idents(cluster3) <- cluster3$predicted.id
pred1 <- WhichCells(cluster3, idents = '1')

int1 <- intersect(sub1, pred1) #5894

#cluster 2
Idents(cluster3) <- cluster3$seurat_clusters
sub2 <- WhichCells(cluster3, idents = c('3'))

Idents(cluster3) <- cluster3$predicted.id
pred2 <- WhichCells(cluster3, idents = '2')

int2 <- intersect(sub2, pred2) #5193

#cluster 3
Idents(cluster3) <- cluster3$seurat_clusters
sub3 <- WhichCells(cluster3, idents = c('1'))

Idents(cluster3) <- cluster3$predicted.id
pred3 <- WhichCells(cluster3, idents = '3')

int3 <- intersect(sub3, pred3) #5742

#cluster 4
Idents(cluster3) <- cluster3$seurat_clusters
sub4 <- WhichCells(cluster3, idents = c('4'))

Idents(cluster3) <- cluster3$predicted.id
pred4 <- WhichCells(cluster3, idents = '4')

int4 <- intersect(sub4, pred4) #3261

#cluster 5
Idents(cluster3) <- cluster3$seurat_clusters
sub5 <- WhichCells(cluster3, idents = c('5'))

Idents(cluster3) <- cluster3$predicted.id
pred5 <- WhichCells(cluster3, idents = '5')

int5 <- intersect(sub5, pred5) #2073

#cluster 6
Idents(cluster3) <- cluster3$seurat_clusters
sub6 <- WhichCells(cluster3, idents = c('6'))

Idents(cluster3) <- cluster3$predicted.id
pred6 <- WhichCells(cluster3, idents = '6')

int6 <- intersect(sub6, pred6) #2014

#cluster 7 
Idents(cluster3) <- cluster3$seurat_clusters
sub7 <- WhichCells(cluster3, idents = c('8'))

Idents(cluster3) <- cluster3$predicted.id
pred7 <- WhichCells(cluster3, idents = '7')

int7 <- intersect(sub7, pred7) #1094

#cluster 8 
Idents(cluster3) <- cluster3$seurat_clusters
sub8 <- WhichCells(cluster3, idents = c('7'))

Idents(cluster3) <- cluster3$predicted.id
pred8 <- WhichCells(cluster3, idents = '8')

int8 <- intersect(sub8, pred8) #1046

#cluster 9 
Idents(cluster3) <- cluster3$seurat_clusters
sub9 <- WhichCells(cluster3, idents = c('9'))

Idents(cluster3) <- cluster3$predicted.id
pred9 <- WhichCells(cluster3, idents = '9')

int9 <- intersect(sub9, pred9) #409

#cluster 10 
Idents(cluster3) <- cluster3$seurat_clusters
sub10 <- WhichCells(cluster3, idents = c('10'))

Idents(cluster3) <- cluster3$predicted.id
pred10 <- WhichCells(cluster3, idents = '10')

int10 <- intersect(sub10, pred10) #367

#cluster 11
Idents(cluster3) <- cluster3$seurat_clusters
sub11 <- WhichCells(cluster3, idents = c('13'))

Idents(cluster3) <- cluster3$predicted.id
pred11 <- WhichCells(cluster3, idents = '11')

int11 <- intersect(sub11, pred11) #306

#cluster 12
Idents(cluster3) <- cluster3$seurat_clusters
sub12 <- WhichCells(cluster3, idents = c('11'))

Idents(cluster3) <- cluster3$predicted.id
pred12 <- WhichCells(cluster3, idents = '12')

int12 <- intersect(sub12, pred12) #309

#cluster 13
Idents(cluster3) <- cluster3$seurat_clusters
sub13 <- WhichCells(cluster3, idents = c('12'))

Idents(cluster3) <- cluster3$predicted.id
pred13 <- WhichCells(cluster3, idents = '13')

int13 <- intersect(sub13, pred13) #144

#cluster 14
Idents(cluster3) <- cluster3$seurat_clusters
sub14 <- WhichCells(cluster3, idents = c('12'))

Idents(cluster3) <- cluster3$predicted.id
pred14 <- WhichCells(cluster3, idents = '14')

int14 <- intersect(sub14, pred14) #153

#cluster 15
Idents(cluster3) <- cluster3$seurat_clusters
sub15 <- WhichCells(cluster3, idents = c('4'))

Idents(cluster3) <- cluster3$predicted.id
pred15 <- WhichCells(cluster3, idents = '15')

int15 <- intersect(sub15, pred15) #87

#total 
#pre:44178
#41016
total3 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total3,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query3_low_res_barcodes_11_2_24.csv', row.names=F)

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


#-----Query 4------#

cluster4 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub2_res20_2_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #13029

#cluster 1
Idents(cluster4) <- cluster4$seurat_clusters
sub1 <- WhichCells(cluster4, idents = c('2'))

Idents(cluster4) <- cluster4$predicted.id
pred1 <- WhichCells(cluster4, idents = '1')

int1 <- intersect(sub1, pred1) #5868

#cluster 2
Idents(cluster4) <- cluster4$seurat_clusters
sub2 <- WhichCells(cluster4, idents = c('0', '3'))

Idents(cluster4) <- cluster4$predicted.id
pred2 <- WhichCells(cluster4, idents = '2')

int2 <- intersect(sub2, pred2) #5912

#cluster 3
Idents(cluster4) <- cluster4$seurat_clusters
sub3 <- WhichCells(cluster4, idents = c('1'))

Idents(cluster4) <- cluster4$predicted.id
pred3 <- WhichCells(cluster4, idents = '3')

int3 <- intersect(sub3, pred3) #5197

#cluster 4
Idents(cluster4) <- cluster4$seurat_clusters
sub4 <- WhichCells(cluster4, idents = c('4'))

Idents(cluster4) <- cluster4$predicted.id
pred4 <- WhichCells(cluster4, idents = '4')

int4 <- intersect(sub4, pred4) #3275

#cluster 5
Idents(cluster4) <- cluster4$seurat_clusters
sub5 <- WhichCells(cluster4, idents = c('5'))

Idents(cluster4) <- cluster4$predicted.id
pred5 <- WhichCells(cluster4, idents = '5')

int5 <- intersect(sub5, pred5) #2450

#cluster 6
Idents(cluster4) <- cluster4$seurat_clusters
sub6 <- WhichCells(cluster4, idents = c('6'))

Idents(cluster4) <- cluster4$predicted.id
pred6 <- WhichCells(cluster4, idents = '6')

int6 <- intersect(sub6, pred6) #2034

#cluster 7
Idents(cluster4) <- cluster4$seurat_clusters
sub7 <- WhichCells(cluster4, idents = c('8'))

Idents(cluster4) <- cluster4$predicted.id
pred7 <- WhichCells(cluster4, idents = '7')

int7 <- intersect(sub7, pred7) #1036

#cluster 8 
Idents(cluster4) <- cluster4$seurat_clusters
sub8 <- WhichCells(cluster4, idents = c('9'))

Idents(cluster4) <- cluster4$predicted.id
pred8 <- WhichCells(cluster4, idents = '8')

int8 <- intersect(sub8, pred8) #916

#cluster 9 
Idents(cluster4) <- cluster4$seurat_clusters
sub9 <- WhichCells(cluster4, idents = c('12'))

Idents(cluster4) <- cluster4$predicted.id
pred9 <- WhichCells(cluster4, idents = '9')

int9 <- intersect(sub9, pred9) #293

#cluster 10 
Idents(cluster4) <- cluster4$seurat_clusters
sub10 <- WhichCells(cluster4, idents = c('11'))

Idents(cluster4) <- cluster4$predicted.id
pred10 <- WhichCells(cluster4, idents = '10')

int10 <- intersect(sub10, pred10) #364

#cluster 11
Idents(cluster4) <- cluster4$seurat_clusters
sub11 <- WhichCells(cluster4, idents = c('10'))

Idents(cluster4) <- cluster4$predicted.id
pred11 <- WhichCells(cluster4, idents = '11')

int11 <- intersect(sub11, pred11) #312

#cluster 12
Idents(cluster4) <- cluster4$seurat_clusters
sub12 <- WhichCells(cluster4, idents = c('13'))

Idents(cluster4) <- cluster4$predicted.id
pred12 <- WhichCells(cluster4, idents = '12')

int12 <- intersect(sub12, pred12) #276

#cluster 13
Idents(cluster4) <- cluster4$seurat_clusters
sub13 <- WhichCells(cluster4, idents = c('14'))

Idents(cluster4) <- cluster4$predicted.id
pred13 <- WhichCells(cluster4, idents = '13')

int13 <- intersect(sub13, pred13) #159

#cluster 14
Idents(cluster4) <- cluster4$seurat_clusters
sub14 <- WhichCells(cluster4, idents = c('10'))

Idents(cluster4) <- cluster4$predicted.id
pred14 <- WhichCells(cluster4, idents = '14')

int14 <- intersect(sub14, pred14) #124

#cluster 15
Idents(cluster4) <- cluster4$seurat_clusters
sub15 <- WhichCells(cluster4, idents = c('15'))

Idents(cluster4) <- cluster4$predicted.id
pred15 <- WhichCells(cluster4, idents = '15')

int15 <- intersect(sub15, pred15) #85

#total 
#pre:44177
#41330
total4 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total4,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query4_low_res_barcodes_11_2_24.csv', row.names=F)

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


#-----Query 5------#

cluster5 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub_res20_3_low_res_cluster_pred.rds")

cluster5$predicted.id <- factor(cluster5$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster5, reduction = 'umap', label = TRUE) + DimPlot(cluster5, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

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

int0 <- intersect(sub0, pred0) #13029

#cluster 1
Idents(cluster5) <- cluster5$seurat_clusters
sub1 <- WhichCells(cluster5, idents = c('2'))

Idents(cluster5) <- cluster5$predicted.id
pred1 <- WhichCells(cluster5, idents = '1')

int1 <- intersect(sub1, pred1) #5883

#cluster 2
Idents(cluster5) <- cluster5$seurat_clusters
sub2 <- WhichCells(cluster5, idents = c('1'))

Idents(cluster5) <- cluster5$predicted.id
pred2 <- WhichCells(cluster5, idents = '2')

int2 <- intersect(sub2, pred2) #5703

#cluster 3
Idents(cluster5) <- cluster5$seurat_clusters
sub3 <- WhichCells(cluster5, idents = c('3', '5'))

Idents(cluster5) <- cluster5$predicted.id
pred3 <- WhichCells(cluster5, idents = '3')

int3 <- intersect(sub3, pred3) #5314

#cluster 4
Idents(cluster5) <- cluster5$seurat_clusters
sub4 <- WhichCells(cluster5, idents = c('4'))

Idents(cluster5) <- cluster5$predicted.id
pred4 <- WhichCells(cluster5, idents = '4')

int4 <- intersect(sub4, pred4) #3290

#cluster 5
Idents(cluster5) <- cluster5$seurat_clusters
sub5 <- WhichCells(cluster5, idents = c('6'))

Idents(cluster5) <- cluster5$predicted.id
pred5 <- WhichCells(cluster5, idents = '5')

int5 <- intersect(sub5, pred5) #2347

#cluster 6
Idents(cluster5) <- cluster5$seurat_clusters
sub6 <- WhichCells(cluster5, idents = c('7'))

Idents(cluster5) <- cluster5$predicted.id
pred6 <- WhichCells(cluster5, idents = '6')

int6 <- intersect(sub6, pred6) #2056

#cluster 7
Idents(cluster5) <- cluster5$seurat_clusters
sub7 <- WhichCells(cluster5, idents = c('9'))

Idents(cluster5) <- cluster5$predicted.id
pred7 <- WhichCells(cluster5, idents = '7')

int7 <- intersect(sub7, pred7) #1037

#cluster 8 
Idents(cluster5) <- cluster5$seurat_clusters
sub8 <- WhichCells(cluster5, idents = c('8'))

Idents(cluster5) <- cluster5$predicted.id
pred8 <- WhichCells(cluster5, idents = '8')

int8 <- intersect(sub8, pred8) #1034

#cluster 9 
Idents(cluster5) <- cluster5$seurat_clusters
sub9 <- WhichCells(cluster5, idents = c('13'))

Idents(cluster5) <- cluster5$predicted.id
pred9 <- WhichCells(cluster5, idents = '9')

int9 <- intersect(sub9, pred9) #295

#cluster 10 
Idents(cluster5) <- cluster5$seurat_clusters
sub10 <- WhichCells(cluster5, idents = c('11'))

Idents(cluster5) <- cluster5$predicted.id
pred10 <- WhichCells(cluster5, idents = '10')

int10 <- intersect(sub10, pred10) #373

#cluster 11
Idents(cluster5) <- cluster5$seurat_clusters
sub11 <- WhichCells(cluster5, idents = c('10'))

Idents(cluster5) <- cluster5$predicted.id
pred11 <- WhichCells(cluster5, idents = '11')

int11 <- intersect(sub11, pred11) #303

#cluster 12
Idents(cluster5) <- cluster5$seurat_clusters
sub12 <- WhichCells(cluster5, idents = c('12'))

Idents(cluster5) <- cluster5$predicted.id
pred12 <- WhichCells(cluster5, idents = '12')

int12 <- intersect(sub12, pred12) #297

#cluster 13
Idents(cluster5) <- cluster5$seurat_clusters
sub13 <- WhichCells(cluster5, idents = c('14'))

Idents(cluster5) <- cluster5$predicted.id
pred13 <- WhichCells(cluster5, idents = '13')

int13 <- intersect(sub13, pred13) #154

#cluster 14
Idents(cluster5) <- cluster5$seurat_clusters
sub14 <- WhichCells(cluster5, idents = c('10'))

Idents(cluster5) <- cluster5$predicted.id
pred14 <- WhichCells(cluster5, idents = '14')

int14 <- intersect(sub14, pred14) #145

#cluster 15
Idents(cluster5) <- cluster5$seurat_clusters
sub15 <- WhichCells(cluster5, idents = c('4'))

Idents(cluster5) <- cluster5$predicted.id
pred15 <- WhichCells(cluster5, idents = '15')

int15 <- intersect(sub15, pred15) #95

#total 
#pre:44178
#40609
total5 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total5,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query5_low_res_barcodes_11_2_24.csv', row.names=F)

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


#-----Query 6------#

cluster6 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub2_res20_3_low_res_cluster_pred.rds")

cluster6$predicted.id <- factor(cluster6$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster6, reduction = 'umap', label = TRUE) + DimPlot(cluster6, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

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

int0 <- intersect(sub0, pred0) #13243

#cluster 1
Idents(cluster6) <- cluster6$seurat_clusters
sub1 <- WhichCells(cluster6, idents = c('1'))

Idents(cluster6) <- cluster6$predicted.id
pred1 <- WhichCells(cluster6, idents = '1')

int1 <- intersect(sub1, pred1) #6572

#cluster 2
Idents(cluster6) <- cluster6$seurat_clusters
sub2 <- WhichCells(cluster6, idents = c('3'))

Idents(cluster6) <- cluster6$predicted.id
pred2 <- WhichCells(cluster6, idents = '2')

int2 <- intersect(sub2, pred2) #5283

#cluster 3
Idents(cluster6) <- cluster6$seurat_clusters
sub3 <- WhichCells(cluster6, idents = c('2'))

Idents(cluster6) <- cluster6$predicted.id
pred3 <- WhichCells(cluster6, idents = '3')

int3 <- intersect(sub3, pred3) #5223

#cluster 4
Idents(cluster6) <- cluster6$seurat_clusters
sub4 <- WhichCells(cluster6, idents = c('4'))

Idents(cluster6) <- cluster6$predicted.id
pred4 <- WhichCells(cluster6, idents = '4')

int4 <- intersect(sub4, pred4) #3256

#cluster 5
Idents(cluster6) <- cluster6$seurat_clusters
sub5 <- WhichCells(cluster6, idents = c('5'))

Idents(cluster6) <- cluster6$predicted.id
pred5 <- WhichCells(cluster6, idents = '5')

int5 <- intersect(sub5, pred5) #2293

#cluster 6
Idents(cluster6) <- cluster6$seurat_clusters
sub6 <- WhichCells(cluster6, idents = c('6'))

Idents(cluster6) <- cluster6$predicted.id
pred6 <- WhichCells(cluster6, idents = '6')

int6 <- intersect(sub6, pred6) #2027

#cluster 7
Idents(cluster6) <- cluster6$seurat_clusters
sub7 <- WhichCells(cluster6, idents = c('7'))

Idents(cluster6) <- cluster6$predicted.id
pred7 <- WhichCells(cluster6, idents = '7')

int7 <- intersect(sub7, pred7) #1091

#cluster 8 
Idents(cluster6) <- cluster6$seurat_clusters
sub8 <- WhichCells(cluster6, idents = c('8'))

Idents(cluster6) <- cluster6$predicted.id
pred8 <- WhichCells(cluster6, idents = '8')

int8 <- intersect(sub8, pred8) #1030

#cluster 9
Idents(cluster6) <- cluster6$seurat_clusters
sub9 <- WhichCells(cluster6, idents = c('11'))

Idents(cluster6) <- cluster6$predicted.id
pred9 <- WhichCells(cluster6, idents = '9')

int9 <- intersect(sub9, pred9) #352

#cluster 10
Idents(cluster6) <- cluster6$seurat_clusters
sub10 <- WhichCells(cluster6, idents = c('9'))

Idents(cluster6) <- cluster6$predicted.id
pred10 <- WhichCells(cluster6, idents = '10')

int10 <- intersect(sub10, pred10) #366

#cluster 11
Idents(cluster6) <- cluster6$seurat_clusters
sub11 <- WhichCells(cluster6, idents = c('10'))

Idents(cluster6) <- cluster6$predicted.id
pred11 <- WhichCells(cluster6, idents = '11')

int11 <- intersect(sub11, pred11) #321

#cluster 12
Idents(cluster6) <- cluster6$seurat_clusters
sub12 <- WhichCells(cluster6, idents = c('13'))

Idents(cluster6) <- cluster6$predicted.id
pred12 <- WhichCells(cluster6, idents = '12')

int12 <- intersect(sub12, pred12) #285

#cluster 13
Idents(cluster6) <- cluster6$seurat_clusters
sub13 <- WhichCells(cluster6, idents = c('12'))

Idents(cluster6) <- cluster6$predicted.id
pred13 <- WhichCells(cluster6, idents = '13')

int13 <- intersect(sub13, pred13) #162

#cluster 14
Idents(cluster6) <- cluster6$seurat_clusters
sub14 <- WhichCells(cluster6, idents = c('12'))

Idents(cluster6) <- cluster6$predicted.id
pred14 <- WhichCells(cluster6, idents = '14')

int14 <- intersect(sub14, pred14) #131

#cluster 15
Idents(cluster6) <- cluster6$seurat_clusters
sub15 <- WhichCells(cluster6, idents = c('4'))

Idents(cluster6) <- cluster6$predicted.id
pred15 <- WhichCells(cluster6, idents = '15')

int15 <- intersect(sub15, pred15) #82

#total 
#pre:44177
#42017
total6 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total6,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query6_low_res_barcodes_11_2_24.csv', row.names=F)

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

#-----Query 7------#

cluster7 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub_res20_4_low_res_cluster_pred.rds")

cluster7$predicted.id <- factor(cluster7$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster7, reduction = 'umap', label = TRUE) + DimPlot(cluster7, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

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

int0 <- intersect(sub0, pred0) #12031

#cluster 1
Idents(cluster7) <- cluster7$seurat_clusters
sub1 <- WhichCells(cluster7, idents = c('2'))

Idents(cluster7) <- cluster7$predicted.id
pred1 <- WhichCells(cluster7, idents = '1')

int1 <- intersect(sub1, pred1) #5912

#cluster 2
Idents(cluster7) <- cluster7$seurat_clusters
sub2 <- WhichCells(cluster7, idents = c('1'))

Idents(cluster7) <- cluster7$predicted.id
pred2 <- WhichCells(cluster7, idents = '2')

int2 <- intersect(sub2, pred2) #6201

#cluster 3
Idents(cluster7) <- cluster7$seurat_clusters
sub3 <- WhichCells(cluster7, idents = c('3', '4'))

Idents(cluster7) <- cluster7$predicted.id
pred3 <- WhichCells(cluster7, idents = '3')

int3 <- intersect(sub3, pred3) #5573

#cluster 4
Idents(cluster7) <- cluster7$seurat_clusters
sub4 <- WhichCells(cluster7, idents = c('7', '9'))

Idents(cluster7) <- cluster7$predicted.id
pred4 <- WhichCells(cluster7, idents = '4')

int4 <- intersect(sub4, pred4) #3257

#cluster 5
Idents(cluster7) <- cluster7$seurat_clusters
sub5 <- WhichCells(cluster7, idents = c('5'))

Idents(cluster7) <- cluster7$predicted.id
pred5 <- WhichCells(cluster7, idents = '5')

int5 <- intersect(sub5, pred5) #2331

#cluster 6
Idents(cluster7) <- cluster7$seurat_clusters
sub6 <- WhichCells(cluster7, idents = c('6'))

Idents(cluster7) <- cluster7$predicted.id
pred6 <- WhichCells(cluster7, idents = '6')

int6 <- intersect(sub6, pred6) #2067

#cluster 7 
Idents(cluster7) <- cluster7$seurat_clusters
sub7 <- WhichCells(cluster7, idents = c('10'))

Idents(cluster7) <- cluster7$predicted.id
pred7 <- WhichCells(cluster7, idents = '7')

int7 <- intersect(sub7, pred7) #1089

#cluster 8 
Idents(cluster7) <- cluster7$seurat_clusters
sub8 <- WhichCells(cluster7, idents = c('8'))

Idents(cluster7) <- cluster7$predicted.id
pred8 <- WhichCells(cluster7, idents = '8')

int8 <- intersect(sub8, pred8) #1029

#cluster 9 
Idents(cluster7) <- cluster7$seurat_clusters
sub9 <- WhichCells(cluster7, idents = c('13'))

Idents(cluster7) <- cluster7$predicted.id
pred9 <- WhichCells(cluster7, idents = '9')

int9 <- intersect(sub9, pred9) #311

#cluster 10 
Idents(cluster7) <- cluster7$seurat_clusters
sub10 <- WhichCells(cluster7, idents = c('12'))

Idents(cluster7) <- cluster7$predicted.id
pred10 <- WhichCells(cluster7, idents = '10')

int10 <- intersect(sub10, pred10) #388

#cluster 11
Idents(cluster7) <- cluster7$seurat_clusters
sub11 <- WhichCells(cluster7, idents = c('11'))

Idents(cluster7) <- cluster7$predicted.id
pred11 <- WhichCells(cluster7, idents = '11')

int11 <- intersect(sub11, pred11) #284

#cluster 12
Idents(cluster7) <- cluster7$seurat_clusters
sub12 <- WhichCells(cluster7, idents = c('14'))

Idents(cluster7) <- cluster7$predicted.id
pred12 <- WhichCells(cluster7, idents = '12')

int12 <- intersect(sub12, pred12) #292

#cluster 13
Idents(cluster7) <- cluster7$seurat_clusters
sub13 <- WhichCells(cluster7, idents = c('11'))

Idents(cluster7) <- cluster7$predicted.id
pred13 <- WhichCells(cluster7, idents = '13')

int13 <- intersect(sub13, pred13) #163

#cluster 14
Idents(cluster7) <- cluster7$seurat_clusters
sub14 <- WhichCells(cluster7, idents = c('11'))

Idents(cluster7) <- cluster7$predicted.id
pred14 <- WhichCells(cluster7, idents = '14')

int14 <- intersect(sub14, pred14) #123

#cluster 15
Idents(cluster7) <- cluster7$seurat_clusters
sub15 <- WhichCells(cluster7, idents = c('9'))

Idents(cluster7) <- cluster7$predicted.id
pred15 <- WhichCells(cluster7, idents = '15')

int15 <- intersect(sub15, pred15) #91

#total 
#pre:44178
#41142
total7 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total7,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query7_low_res_barcodes_11_2_24.csv', row.names=F)

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


#-----Query 8------#

cluster8 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub2_res20_4_low_res_cluster_pred.rds")

cluster8$predicted.id <- factor(cluster8$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster8, reduction = 'umap', label = TRUE) + DimPlot(cluster8, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

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

int0 <- intersect(sub0, pred0) #13206

#cluster 1
Idents(cluster8) <- cluster8$seurat_clusters
sub1 <- WhichCells(cluster8, idents = c('1'))

Idents(cluster8) <- cluster8$predicted.id
pred1 <- WhichCells(cluster8, idents = '1')

int1 <- intersect(sub1, pred1) #6531

#cluster 2
Idents(cluster8) <- cluster8$seurat_clusters
sub2 <- WhichCells(cluster8, idents = c('3'))

Idents(cluster8) <- cluster8$predicted.id
pred2 <- WhichCells(cluster8, idents = '2')

int2 <- intersect(sub2, pred2) #5592

#cluster 3
Idents(cluster8) <- cluster8$seurat_clusters
sub3 <- WhichCells(cluster8, idents = c('2'))

Idents(cluster8) <- cluster8$predicted.id
pred3 <- WhichCells(cluster8, idents = '3')

int3 <- intersect(sub3, pred3) #5666

#cluster 4
Idents(cluster8) <- cluster8$seurat_clusters
sub4 <- WhichCells(cluster8, idents = c('4'))

Idents(cluster8) <- cluster8$predicted.id
pred4 <- WhichCells(cluster8, idents = '4')

int4 <- intersect(sub4, pred4) #3278

#cluster 5
Idents(cluster8) <- cluster8$seurat_clusters
sub5 <- WhichCells(cluster8, idents = c('5'))

Idents(cluster8) <- cluster8$predicted.id
pred5 <- WhichCells(cluster8, idents = '5')

int5 <- intersect(sub5, pred5) #2373

#cluster 6
Idents(cluster8) <- cluster8$seurat_clusters
sub6 <- WhichCells(cluster8, idents = c('6'))

Idents(cluster8) <- cluster8$predicted.id
pred6 <- WhichCells(cluster8, idents = '6')

int6 <- intersect(sub6, pred6) #2014

#cluster 7
Idents(cluster8) <- cluster8$seurat_clusters
sub7 <- WhichCells(cluster8, idents = c('7'))

Idents(cluster8) <- cluster8$predicted.id
pred7 <- WhichCells(cluster8, idents = '7')

int7 <- intersect(sub7, pred7) #1050

#cluster 8 
Idents(cluster8) <- cluster8$seurat_clusters
sub8 <- WhichCells(cluster8, idents = c('8'))

Idents(cluster8) <- cluster8$predicted.id
pred8 <- WhichCells(cluster8, idents = '8')

int8 <- intersect(sub8, pred8) #1034

#cluster 9 
Idents(cluster8) <- cluster8$seurat_clusters
sub9 <- WhichCells(cluster8, idents = c('11'))

Idents(cluster8) <- cluster8$predicted.id
pred9 <- WhichCells(cluster8, idents = '9')

int9 <- intersect(sub9, pred9) #311

#cluster 10 
Idents(cluster8) <- cluster8$seurat_clusters
sub10 <- WhichCells(cluster8, idents = c('10'))

Idents(cluster8) <- cluster8$predicted.id
pred10 <- WhichCells(cluster8, idents = '10')

int10 <- intersect(sub10, pred10) #344

#cluster 11
Idents(cluster8) <- cluster8$seurat_clusters
sub11 <- WhichCells(cluster8, idents = c('9'))

Idents(cluster8) <- cluster8$predicted.id
pred11 <- WhichCells(cluster8, idents = '11')

int11 <- intersect(sub11, pred11) #338

#cluster 12
Idents(cluster8) <- cluster8$seurat_clusters
sub12 <- WhichCells(cluster8, idents = c('12'))

Idents(cluster8) <- cluster8$predicted.id
pred12 <- WhichCells(cluster8, idents = '12')

int12 <- intersect(sub12, pred12) #294

#cluster 13
Idents(cluster8) <- cluster8$seurat_clusters
sub13 <- WhichCells(cluster8, idents = c('7'))

Idents(cluster8) <- cluster8$predicted.id
pred13 <- WhichCells(cluster8, idents = '13')

int13 <- intersect(sub13, pred13) #143

#cluster 14
Idents(cluster8) <- cluster8$seurat_clusters
sub14 <- WhichCells(cluster8, idents = c('9'))

Idents(cluster8) <- cluster8$predicted.id
pred14 <- WhichCells(cluster8, idents = '14')

int14 <- intersect(sub14, pred14) #149

#cluster 15
Idents(cluster8) <- cluster8$seurat_clusters
sub15 <- WhichCells(cluster8, idents = c('4'))

Idents(cluster8) <- cluster8$predicted.id
pred15 <- WhichCells(cluster8, idents = '15')

int15 <- intersect(sub15, pred15) #87

#total 
#pre:44177
#42410
total8 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total8,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query8_low_res_barcodes_11_2_24.csv', row.names=F)

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


#-----Query 9------#

cluster9 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub_res20_5_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #12490

#cluster 1
Idents(cluster9) <- cluster9$seurat_clusters
sub1 <- WhichCells(cluster9, idents = c('2'))

Idents(cluster9) <- cluster9$predicted.id
pred1 <- WhichCells(cluster9, idents = '1')

int1 <- intersect(sub1, pred1) #6516

#cluster 2
Idents(cluster9) <- cluster9$seurat_clusters
sub2 <- WhichCells(cluster9, idents = c('1'))

Idents(cluster9) <- cluster9$predicted.id
pred2 <- WhichCells(cluster9, idents = '2')

int2 <- intersect(sub2, pred2) #6148

#cluster 3
Idents(cluster9) <- cluster9$seurat_clusters
sub3 <- WhichCells(cluster9, idents = c('3'))

Idents(cluster9) <- cluster9$predicted.id
pred3 <- WhichCells(cluster9, idents = '3')

int3 <- intersect(sub3, pred3) #5540

#cluster 4
Idents(cluster9) <- cluster9$seurat_clusters
sub4 <- WhichCells(cluster9, idents = c('4'))

Idents(cluster9) <- cluster9$predicted.id
pred4 <- WhichCells(cluster9, idents = '4')

int4 <- intersect(sub4, pred4) #3258

#cluster 5
Idents(cluster9) <- cluster9$seurat_clusters
sub5 <- WhichCells(cluster9, idents = c('5'))

Idents(cluster9) <- cluster9$predicted.id
pred5 <- WhichCells(cluster9, idents = '5')

int5 <- intersect(sub5, pred5) #2210

#cluster 6
Idents(cluster9) <- cluster9$seurat_clusters
sub6 <- WhichCells(cluster9, idents = c('6'))

Idents(cluster9) <- cluster9$predicted.id
pred6 <- WhichCells(cluster9, idents = '6')

int6 <- intersect(sub6, pred6) #2025

#cluster 7 
Idents(cluster9) <- cluster9$seurat_clusters
sub7 <- WhichCells(cluster9, idents = c('7'))

Idents(cluster9) <- cluster9$predicted.id
pred7 <- WhichCells(cluster9, idents = '7')

int7 <- intersect(sub7, pred7) #1042

#cluster 8 
Idents(cluster9) <- cluster9$seurat_clusters
sub8 <- WhichCells(cluster9, idents = c('8'))

Idents(cluster9) <- cluster9$predicted.id
pred8 <- WhichCells(cluster9, idents = '8')

int8 <- intersect(sub8, pred8) #886

#cluster 9 
Idents(cluster9) <- cluster9$seurat_clusters
sub9 <- WhichCells(cluster9, idents = c('10'))

Idents(cluster9) <- cluster9$predicted.id
pred9 <- WhichCells(cluster9, idents = '9')

int9 <- intersect(sub9, pred9) #422

#cluster 10
Idents(cluster9) <- cluster9$seurat_clusters
sub10 <- WhichCells(cluster9, idents = c('11'))

Idents(cluster9) <- cluster9$predicted.id
pred10 <- WhichCells(cluster9, idents = '10')

int10 <- intersect(sub10, pred10) #354

#cluster 11
Idents(cluster9) <- cluster9$seurat_clusters
sub11 <- WhichCells(cluster9, idents = c('9'))

Idents(cluster9) <- cluster9$predicted.id
pred11 <- WhichCells(cluster9, idents = '11')

int11 <- intersect(sub11, pred11) #319

#cluster 12
Idents(cluster9) <- cluster9$seurat_clusters
sub12 <- WhichCells(cluster9, idents = c('12'))

Idents(cluster9) <- cluster9$predicted.id
pred12 <- WhichCells(cluster9, idents = '12')

int12 <- intersect(sub12, pred12) #300

#cluster 13
Idents(cluster9) <- cluster9$seurat_clusters
sub13 <- WhichCells(cluster9, idents = c('13'))

Idents(cluster9) <- cluster9$predicted.id
pred13 <- WhichCells(cluster9, idents = '13')

int13 <- intersect(sub13, pred13) #146

#cluster 14
Idents(cluster9) <- cluster9$seurat_clusters
sub14 <- WhichCells(cluster9, idents = c('9'))

Idents(cluster9) <- cluster9$predicted.id
pred14 <- WhichCells(cluster9, idents = '14')

int14 <- intersect(sub14, pred14) #422

#cluster 15
Idents(cluster9) <- cluster9$seurat_clusters
sub15 <- WhichCells(cluster9, idents = c('4'))

Idents(cluster9) <- cluster9$predicted.id
pred15 <- WhichCells(cluster9, idents = '15')

int15 <- intersect(sub15, pred15) #96

#total 
#pre:44178
#41903
total9 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total9,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query9_low_res_barcodes_11_2_24.csv', row.names=F)

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


#-----Query 10------#

cluster10 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/rectum_sub2_res20_5_low_res_cluster_pred.rds")

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
sub0 <- WhichCells(cluster10, idents = c('0'))

Idents(cluster10) <- cluster10$predicted.id
pred0 <- WhichCells(cluster10, idents = '0')

int0 <- intersect(sub0, pred0) #12837

#cluster 1
Idents(cluster10) <- cluster10$seurat_clusters
sub1 <- WhichCells(cluster10, idents = c('1'))

Idents(cluster10) <- cluster10$predicted.id
pred1 <- WhichCells(cluster10, idents = '1')

int1 <- intersect(sub1, pred1) #6549

#cluster 2
Idents(cluster10) <- cluster10$seurat_clusters
sub2 <- WhichCells(cluster10, idents = c('2', '3'))

Idents(cluster10) <- cluster10$predicted.id
pred2 <- WhichCells(cluster10, idents = '2')

int2 <- intersect(sub2, pred2) #6541

#cluster 3
Idents(cluster10) <- cluster10$seurat_clusters
sub3 <- WhichCells(cluster10, idents = c('2', '5'))

Idents(cluster10) <- cluster10$predicted.id
pred3 <- WhichCells(cluster10, idents = '3')

int3 <- intersect(sub3, pred3) #5554

#cluster 4
Idents(cluster10) <- cluster10$seurat_clusters
sub4 <- WhichCells(cluster10, idents = c('4'))

Idents(cluster10) <- cluster10$predicted.id
pred4 <- WhichCells(cluster10, idents = '4')

int4 <- intersect(sub4, pred4) #3250

#cluster 5
Idents(cluster10) <- cluster10$seurat_clusters
sub5 <- WhichCells(cluster10, idents = c('6'))

Idents(cluster10) <- cluster10$predicted.id
pred5 <- WhichCells(cluster10, idents = '5')

int5 <- intersect(sub5, pred5) #2347

#cluster 6
Idents(cluster10) <- cluster10$seurat_clusters
sub6 <- WhichCells(cluster10, idents = c('7'))

Idents(cluster10) <- cluster10$predicted.id
pred6 <- WhichCells(cluster10, idents = '6')

int6 <- intersect(sub6, pred6) #2018

#cluster 7
Idents(cluster10) <- cluster10$seurat_clusters
sub7 <- WhichCells(cluster10, idents = c('8'))

Idents(cluster10) <- cluster10$predicted.id
pred7 <- WhichCells(cluster10, idents = '7')

int7 <- intersect(sub7, pred7) #1093

#cluster 8
Idents(cluster10) <- cluster10$seurat_clusters
sub8 <- WhichCells(cluster10, idents = c('9'))

Idents(cluster10) <- cluster10$predicted.id
pred8 <- WhichCells(cluster10, idents = '8')

int8 <- intersect(sub8, pred8) #1052

#cluster 9
Idents(cluster10) <- cluster10$seurat_clusters
sub9 <- WhichCells(cluster10, idents = c('13'))

Idents(cluster10) <- cluster10$predicted.id
pred9 <- WhichCells(cluster10, idents = '9')

int9 <- intersect(sub9, pred9) #261

#cluster 10 
Idents(cluster10) <- cluster10$seurat_clusters
sub10 <- WhichCells(cluster10, idents = c('11'))

Idents(cluster10) <- cluster10$predicted.id
pred10 <- WhichCells(cluster10, idents = '10')

int10 <- intersect(sub10, pred10) #384

#cluster 11
Idents(cluster10) <- cluster10$seurat_clusters
sub11 <- WhichCells(cluster10, idents = c('10'))

Idents(cluster10) <- cluster10$predicted.id
pred11 <- WhichCells(cluster10, idents = '11')

int11 <- intersect(sub11, pred11) #295

#cluster 12
Idents(cluster10) <- cluster10$seurat_clusters
sub12 <- WhichCells(cluster10, idents = c('12'))

Idents(cluster10) <- cluster10$predicted.id
pred12 <- WhichCells(cluster10, idents = '12')

int12 <- intersect(sub12, pred12) #284

#cluster 13
Idents(cluster10) <- cluster10$seurat_clusters
sub13 <- WhichCells(cluster10, idents = c('10'))

Idents(cluster10) <- cluster10$predicted.id
pred13 <- WhichCells(cluster10, idents = '13')

int13 <- intersect(sub13, pred13) #159

#cluster 14
Idents(cluster10) <- cluster10$seurat_clusters
sub14 <- WhichCells(cluster10, idents = c('14'))

Idents(cluster10) <- cluster10$predicted.id
pred14 <- WhichCells(cluster10, idents = '14')

int14 <- intersect(sub14, pred14) #118

#cluster 15
Idents(cluster10) <- cluster10$seurat_clusters
sub15 <- WhichCells(cluster10, idents = c('4'))

Idents(cluster10) <- cluster10$predicted.id
pred15 <- WhichCells(cluster10, idents = '15')

int15 <- intersect(sub15, pred15) #118

#total 
#pre:44177
#42824
total10 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total10,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query10_low_res_barcodes_11_2_24.csv', row.names=F)

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

#----------identify stably assigned cells--------#

#load in barcodes 
total1 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query1_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)
total2 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query2_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)
total3 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query3_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)
total4 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query4_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)
total5 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query5_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)
total6 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query6_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)
total7 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query7_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)
total8 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query8_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)
total9 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query9_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)
total10 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/query10_low_res_barcodes_11_2_24.csv', what = "", sep = ",", skip = 1)

#merge the barcodes into a list 
#417742 barcodes 
total_merge <- c(total1, total2, total3, total4, total5, total6, total7, total8, total9, total10)

#counts the number of times each barcode is present
value_counts <- table(total_merge)
table(value_counts)

#define a threshold - more times a barcode is present, more stable it is
threshold <- 4

#filter values that meet your threshold 
filtered_values <- names(value_counts[value_counts >= threshold])  

#1: 88264 - keeps ~99.89% of data
#2: 87627 - keeps ~99.17% of data 
#3: 85619 - keeps ~96.90% of data
#4: 81500 - keeps ~92.24% of data ----select threshold of 4 (biggest jump after threshold of 4)
#5: 74732 - keeps ~84.58% of data

#save the filtered barcodes 
write.csv(filtered_values, file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/filtered_low_res_barcodes_11_4_24.csv', row.names = F)

#subset the rectum seurat object on the server 

#load SO
rectum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_low_res_ref_10_21_24.rds")

#load in stably assigned cells 
filtered_values <- scan('/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/filtered_low_res_barcodes_11_4_24.csv', what = "", sep = ",", skip = 1)

#subset rectum based on stable barcodes 
rectum_sub <- subset(rectum, cells = filtered_values)

#save 
saveRDS(rectum_sub, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_low_res_stable_cell_11_4_24.rds")

#---------highlight stable and unstable cells-------#

#load SO
rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_rpca_10_17_24.rds")

#rectum barcodes 
barcodes <- colnames(rectum)

unstable <- setdiff(barcodes, filtered_values)

DimPlot(rectum, reduction = "umap", cells.highlight = unstable, sizes.highlight = 0.1) + scale_color_manual(labels = c("Stably Assigned Cells", "Unstably Assigned Cells"), values = c("grey", "red"))


#-----------Low resolution ARI results------------#
#11/5/2024

#merged ARI scores
df <- read.csv('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/df_ri_helm_batch1_13_rectum_low_res_dune_merge_11_4_24.csv', row.names = NULL, sep = ",")

df <- df[,-1]
names(df) <- gsub("^paste0\\.", "", names(df))

columns <- colnames(df)

#plot the ARI - without paste0
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/figures/RI_helm_batch1_13_rectum_low_res_all_param_merge_11_5_24.pdf", width = 15, height = 15)
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
summarise(ari_values, overall_mean = mean(mean_ari)) #0.8225531

#remove mean_ari from ari_values 
ari_values <- ari_values[, -c(37)]

#median ARI

ari_values$row_median = apply(ari_values, 1, median, na.rm=TRUE)

overall_median <- median(ari_values$row_median)  #0.8544736

#top 10 median ARI scores 
median_ari <- as.data.frame(ari_values$row_median)
rownames(median_ari) <- rownames(ari_values)
top_n(median_ari, 10)

#---------low resolution top 3 clustering results----------#

#-----parameter set 1------#

#parameters: resolution: 0.75, PCs: 15, HVG: 2000

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_75_15_2000_rpca_11_5_24.rds")

rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
rectum$batch <- factor(rectum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(rectum, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(rectum, features = 'nCount_RNA', pt.size = 0)
VlnPlot(rectum, features = 'percent.mt', pt.size = 0)

#visualize PCA results
ElbowPlot(rectum)
DimPlot(rectum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(rectum, reduction = 'pca', group.by = "orig.ident")
DimPlot(rectum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(rectum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(rectum, reduction = 'umap', group.by = "orig.ident")
DimPlot(rectum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
Idents(rectum) <- rectum$seurat_clusters
pt1 <- table(Idents(rectum), rectum$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28'))

library(randomcoloR)
no_of_colors <- 29

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
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
Idents(rectum) <- rectum$batch
pt2 <- table(Idents(rectum), rectum$seurat_clusters)
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


#-----parameter set 2------#

#parameters: resolution: 1, PCs: 15, HVG: 2000

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_1_15_2000_rpca_11_5_24.rds")

rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
rectum$batch <- factor(rectum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(rectum, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(rectum, features = 'nCount_RNA', pt.size = 0)
VlnPlot(rectum, features = 'percent.mt', pt.size = 0)

#visualize PCA results
ElbowPlot(rectum)
DimPlot(rectum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(rectum, reduction = 'pca', group.by = "orig.ident")
DimPlot(rectum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(rectum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(rectum, reduction = 'umap', group.by = "orig.ident")
DimPlot(rectum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
Idents(rectum) <- rectum$seurat_clusters
pt1 <- table(Idents(rectum), rectum$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33'))

library(randomcoloR)
no_of_colors <- 34

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
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
Idents(rectum) <- rectum$batch
pt2 <- table(Idents(rectum), rectum$seurat_clusters)
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

#-----parameter set 3------#

#parameters: resolution: 0.5, PCs: 15, HVG: 2000

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_50_15_2000_rpca_11_5_24.rds")

rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
rectum$batch <- factor(rectum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(rectum, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(rectum, features = 'nCount_RNA', pt.size = 0)
VlnPlot(rectum, features = 'percent.mt', pt.size = 0)

#visualize PCA results
ElbowPlot(rectum)
DimPlot(rectum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(rectum, reduction = 'pca', group.by = "orig.ident")
DimPlot(rectum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(rectum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(rectum, reduction = 'umap', group.by = "orig.ident")
DimPlot(rectum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
Idents(rectum) <- rectum$seurat_clusters
pt1 <- table(Idents(rectum), rectum$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", '14', '15', '16', '17', '18', '19', '20', '21'))

library(randomcoloR)
no_of_colors <- 22

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
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
Idents(rectum) <- rectum$batch
pt2 <- table(Idents(rectum), rectum$seurat_clusters)
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


#subset cluster 28 
rectum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_75_15_2000_rpca_11_5_24.rds")
Idents(rectum) <- rectum$seurat_clusters
rectum_sub <- subset(rectum, ident = "28", invert = T)
saveRDS(rectum_sub, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_75_15_2000_rpca_sub_11_5_24.rds")

#-----parameter set 1 subset doublet cluster------#

#parameters: resolution: 0.75, PCs: 15, HVG: 2000

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_75_15_2000_rpca_sub_cluster_11_6_24.rds")

rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
rectum$batch <- factor(rectum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(rectum, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(rectum, features = 'nCount_RNA', pt.size = 0)
VlnPlot(rectum, features = 'percent.mt', pt.size = 0)

#visualize PCA results
ElbowPlot(rectum)
DimPlot(rectum, reduction = 'pca', group.by = "seurat_clusters")
DimPlot(rectum, reduction = 'pca', group.by = "orig.ident")
DimPlot(rectum, reduction = 'pca', group.by = "batch")

#visualize UMAP
DimPlot(rectum, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(rectum, reduction = 'umap', group.by = "orig.ident")
DimPlot(rectum, reduction = 'umap', group.by = "batch")

#proportion of cells per cluster 
rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
Idents(rectum) <- rectum$seurat_clusters
pt1 <- table(Idents(rectum), rectum$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27'))

library(randomcoloR)
no_of_colors <- 28

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

plot(ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
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
Idents(rectum) <- rectum$batch
pt2 <- table(Idents(rectum), rectum$seurat_clusters)
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


#------preliminary annotation------#

#11/6/2024

#set assay to RNA and make sure Idents are set to cluster
DefaultAssay(rectum) <- "RNA"

Idents(rectum) <- rectum$seurat_clusters

#-----T cells-----#

#naive/central memory CD4 T cell
DotPlot(rectum, features = c("CD3D", 'CD4', 'CCR7', 'SELL'))

#tissue resident memory CD4 T cell
DotPlot(rectum, features = c("CD3D", 'CD4', 'ITGAE', 'ITGA1', 'SPRY1'))

#tissue resident memory T helper 17 
DotPlot(rectum, features = c("CD3D", 'CD4', 'ITGAE', 'ITGA1', 'SPRY1', 'IL17A', 'CCR6', 'CCL20'))

#naive t follicular helper 
DotPlot(rectum, features = c("CD3D", 'CD4', 'SELL', 'CXCR5', 'BCL6', 'CD40LG'))

#t follicular helper cell 
DotPlot(rectum, features = c("CD3D", 'CD4', 'CXCR5', 'BCL6', 'PDCD1', 'CD40LG'))

#regulatory t cell 
DotPlot(rectum, features = c("CD3D", 'CD4', 'FOXP3', 'CTLA4', 'IL2RA'))

#regulatory t cell IL10+
DotPlot(rectum, features = c("CD3D", 'CD4', 'FOXP3', 'CTLA4', 'IL2RA', 'IL10'))

#naive/central memory CD8 t cell 
DotPlot(rectum, features = c("CD3D", 'CD8A', 'CD8B', 'CCR7', 'SELL'))

#tissue resident memory CD8 t cell 
DotPlot(rectum, features = c("CD3D", 'CD8A', 'CD8B', 'ITGAE', 'ITGA1', 'SPRY1'))

#tissue resident memory/effector memory Cd8 t cell 
DotPlot(rectum, features = c("CD3D", 'CD8A', 'CD8B', 'GZMK', 'CRTAM', 'EOMES'))

#naive gamma delta t cell 
DotPlot(rectum, features = c("CD3D", 'SELL', 'KLRC2', 'TRDC', 'TRGC1', 'KIR2DL4', 'CCR7'))

#gamma delta t cell 
DotPlot(rectum, features = c("CD3D", 'KLRC2', 'TRDC', 'TRGC1', 'KIR2DL4'))

#MAIT (mucosal associated invariant t cell)
DotPlot(rectum, features = c("CD3D", 'SLC4A10', 'TRAV1-2'))

#CD16+ NK cell 
DotPlot(rectum, features = c('KLRF1', 'NKG7', 'FCGR3A', 'GNLY', 'GZMB'))

#CD56 bright natural killer cell 
DotPlot(rectum, features = c('KLRF1', 'NKG7', 'XCL1', 'IL2RB', 'NCR1', 'FCER1G', 'NCAM1'))

#ILC3 innate lymphoid cell type 3 (cd3d negative)
DotPlot(rectum, features = c('IL7R', 'RORC', 'KIT', 'LST1', 'PCDH9', 'IL1R1', 'IL23R'))

#cycling t or NK cell
DotPlot(rectum, features = c('CD3D', 'MKI67', 'TOP2A'))

#--------myeloid------------#

#type 1 conventional DC 
DotPlot(seur_obj, features = c('HLA-DRA', 'HLA-DPA1', 'CLEC9A', 'XCR1', 'BATF3', 'CADM1', 'RAB7B'))

#type 2 conventional DC
DotPlot(seur_obj, features = c('HLA-DRA', 'HLA-DPA1', 'CLEC10A', 'FCER1A', 'CD1C'))

#migratory DC
DotPlot(seur_obj, features = c('HLA-DRA', 'HLA-DPA1', 'CCR7', 'LAMP3'))

#plasmacytoid DC
DotPlot(seur_obj, features = c('HLA-DRA', 'HLA-DPA1', 'IRF7', 'CLEC4C', 'JCHAIN', 'LILRA4', 'GZMB'))

#langerhans DC
DotPlot(seur_obj, features = c('HLA-DRA', 'HLA-DPA1', 'ITGAX', 'IL22RA2', 'CD207', 'RUNX3'))

#monocyte
DotPlot(rectum, features = c('FCN1', 'S100A8', 'S100A9', 'IL1B', 'EREG', 'NAMPT', 'PLAUR', 'VCAN', 'FPR1', 'CD3D00E'))

#LYVE1 macrophage
DotPlot(seur_obj, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'LYVE1', 'RNASE1', 'FOLR2'))

#MMP9 macrophage
DotPlot(seur_obj, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'MMP9', 'PLA2G2D', 'ADAMDEC1'))

#TREM2 macrophage 
DotPlot(seur_obj, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'TREM2', 'ACP5', 'CTSD', 'CSTB'))

#CD5L macrophage
DotPlot(seur_obj, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'CD5L', 'VCAM1', 'CXCL12', 'PDK4', 'RBP7'))

#macrophage
DotPlot(rectum, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'CD209'))

#mast cell 
DotPlot(rectum, features = c('CD69', 'KIT', 'TPSB2', 'TPSAB1'))

#eosinophil/basophil
DotPlot(seur_obj, features = c('GATA2', 'CNRIP1', 'PRG2', 'GIHCG', 'CLC'))

#erythrocytes
DotPlot(rectum, features = c('GATA1', 'HBZ', 'HBE1', 'HBG1'))

#monocyte/neutrophil progenitor
DotPlot(seur_obj, features = c('FCN1', 'S100A8', 'S100A9', 'MPO', 'RETN', 'RNASE2', 'PCLAF'))

#megakaryocyte/platelet
DotPlot(seur_obj, features = c('GATA1', 'TAL1', 'MMRN1', 'CMTM5', 'MPIG6B', 'ITGA2B', 'PF4'))

#----------b cells-------------#

#pre b cell 
DotPlot(rectum, features = c('CD19', 'HLA-DRA', 'CD79B', 'SPIB', 'TCL1A', 'CD3D7'))

#pro b cell 
DotPlot(rectum, features = c('CD19', 'HLA-DRA', 'CD79B', 'IGLL1', 'RAG1', 'DNTT', 'VPREB3'))

#naive b cell 
DotPlot(rectum, features = c('CD19', 'HLA-DRA', 'CD79A', 'MS4A1', 'SELL', 'TCL1A', 'IGHD'))

#memory b cell 
DotPlot(rectum, features = c('CD19', 'HLA-DRA', 'CD79A', 'MS4A1', 'CD27', 'TNFSF13B'))

#germinal center b cell 1
DotPlot(rectum, features = c('CD19', 'HLA-DRA', 'CD79A', 'MS4A1', 'MKI67', 'AICDA', 'BCL6', 'SUGCT'))

#plasmablast
DotPlot(rectum, features = c('MZB1', 'JCHAIN', 'XBP1'))

#IgM plasma cell 
DotPlot(seur_obj, features = c('MZB1', 'JCHAIN', 'IGHM'))

#IgA plasma cell 
DotPlot(rectum, features = c('MZB1', 'JCHAIN', 'IGHA1', 'IGHA2'))

#IgG plasma cell 
DotPlot(rectum, features = c('MZB1', 'JCHAIN', 'IGHG3', 'IGHG1', 'IGHG2', 'IGHG4'))

#---------endothelial cells-----------#

#capillary endothelial cell 
DotPlot(rectum, features = c('PECAM1', 'CD3D6', 'RGCC', 'COL4A1', 'COL4A2', 'IL32', 'MCAM', 'MYO1B'))

#arterial endothelial cell 
DotPlot(seur_obj, features = c('PECAM1', 'CD3D6', 'GJA4', 'HEY1', 'CXCL12', 'SEMA3G', 'IGFBP3', 'FBLN2', 'FBLN5', 'ELN', 'BTNL9', 'ALPL'))

#venous endothelial cell
DotPlot(rectum, features = c('PECAM1', 'CD3D6', 'ACKR1', 'CCL14', 'SELE', 'TNFRSF6B'))

#lymphatic endothelial cell 
DotPlot(seur_obj, features = c('PECAM1', 'CD3D6', 'CCL21', 'TFF3', 'PROX2', 'NTS'))

#cycling endothelial cell 
DotPlot(seur_obj, features = c('PECAM1', 'CD3D6', 'MKI67', 'TOP2A'))

#-------------mesenchymal cells-------------#

#vascular smooth muscle cell 
DotPlot(rectum, features = c('VIM', 'DCN', 'TAGLN', 'ACTA2', 'TPM2', 'MYH11', 'RERGL', 'MUSTN1', 'LBH', 'NET1', 'MAP3K20'))

#pericyte
DotPlot(seur_obj, features = c('VIM', 'DCN', 'TAGLN', 'ACTA2', 'TPM2', 'MYH11', 'COX4I2', 'HIGD1B', 'RGS5', 'NDUFA4L2'))

#immune recruiting pericyte
DotPlot(seur_obj, features = c('VIM', 'DCN', 'TAGLN', 'ACTA2', 'TPM2', 'MYH11', 'GPC3', 'COL14A1', 'ECRG4', 'ID4', 'FHL2', 'CXCL12'))

#myofibroblast
DotPlot(seur_obj, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'ACTG2', 'HHIP', 'SOSTDC1', 'NPNT'))

#follicular DC
DotPlot(seur_obj, features = c('VIM', 'FDCSP', 'SRGN', 'CR2', 'CLU', 'CSTA'))

#reticular fibroblast 
DotPlot(seur_obj, features = c('VIM', 'CCL21', 'CCL19', 'TNFSF13B', 'TDO2'))

#oral mucosa fibroblast
DotPlot(seur_obj, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'CTHRC1', 'COL12A1', 'COL1A1', 'CTSK', 'COL5A2'))

#oesophagus fibroblast
DotPlot(seur_obj, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'APOD', 'PLPP1', 'MFAP4', 'IFITM1', 'RASD1'))

#crypt fibroblast
DotPlot(seur_obj, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'PI16', 'RSPO3', 'SFRP1', 'TM2A'))

#villus fibroblast
DotPlot(seur_obj, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'F3', 'PLAT', 'HSD17B2', 'SOX6'))

#lamina propria fibroblast
DotPlot(seur_obj, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'ADAMDEC1', 'ADAM28', 'CCL11', 'CCL8', 'CCL13', 'CFD'))

#rectum fibroblast
DotPlot(rectum, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'KCNN3', 'THBS4', 'FNDC1', 'PPFIBP1'))

#mesothelium 
DotPlot(seur_obj, features = c('VIM', 'DCN', 'UPK3B', 'MSLN', 'SLPI', 'PLAT', 'KRT19'))

#-----------epithelium-------------#

#---------large intestine----------#

#colonocyte
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'CA1', 'GPT'))

#BEST4 colonocyte
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'BEST4', 'CA7', 'OTOP2'))

#enteroendocrine cell 
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'CHGA', 'PCSK1N', 'SCT', 'SCGN', 'NEUROD1'))

#goblet cell 
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'MUC2', 'TFF3', 'FCGBP', 'ZG16'))

#mature colonocyte
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'SLC26A3', 'AQP8', 'CEACAM7'))

#deep crypt secretory cell 
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'MUC17', 'TFF1', 'CD55', 'TM4SF1', 'DUOX2', 'DUOXA2'))

#tuft cell 
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'SH2D6', 'LRMP', 'MATK', 'FYB1', 'HPGDS', 'POU2F3', 'TRPM5'))

#enterocyte
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'FABP2', 'APOA1', 'ALDOB'))

#transit amplifying cell 
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'MKI67', 'TOP2A', 'PCLAF', 'PCNA'))

#stem cell 
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'LGR5', 'RGMB', 'ASCL2', 'OLFM4'))

#enteroendocrine progenitor
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'SOX4', 'SCGN', 'NEUROD1'))

#proliferating deep crypt secretory cell
DotPlot(rectum, features = c('CDH1', 'KRT19', 'EPCAM', 'MUC17', 'TFF1', 'MKI67', 'PCLAF', 'PCNA', 'TOP2A'))

#---------comparisons--------#

#T cells 

#cluster 8 vs. cluster 14 
c8_14 <- FindMarkers(rectum, ident.1 = "8", ident.2 = "14")

write.csv(c8_14, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster8_vs_cluster14_param1_11_9_25.csv")

#cluster 8 vs. 19 
c8_19 <- FindMarkers(rectum, ident.1 = "8", ident.2 = "19")

write.csv(c8_19, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster8_vs_cluster19_param1_11_9_25.csv")

#cluster 14 vs 19 
c14_19 <- FindMarkers(rectum, ident.1 = "14", ident.2 = "19")

write.csv(c14_19, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster14_vs_cluster19_param1_11_9_25.csv")


#cluster 8 and 14 vs 18 
c18_8_14 <- FindMarkers(rectum, ident.1 = "18", ident.2 = c("8", '14'))

write.csv(c18_8_14, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster18_vs_cluster8_14_param1_11_9_25.csv")

#plasma cell 
c10_16 <- FindMarkers(rectum, ident.1 = "10", ident.2 = "16")

write.csv(c10_16, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster10_vs_cluster16_param1_11_9_25.csv")

#goblet cell 
c9_15 <- FindMarkers(rectum, ident.1 = "9", ident.2 = "15")

write.csv(c9_15, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster9_vs_cluster15_param1_11_9_25.csv")

#mature colonocytes
c12_17 <- FindMarkers(rectum, ident.1 = "12", ident.2 = "17")

write.csv(c12_17, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster12_vs_cluster17_param1_11_9_25.csv")

c6_17 <- FindMarkers(rectum, ident.1 = "6", ident.2 = "17")

write.csv(c6_17, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster6_vs_cluster17_param1_11_9_25.csv")

#stem and TA cells 
c3_4 <- FindMarkers(rectum, ident.1 = "3", ident.2 = "4")

write.csv(c3_4, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster3_vs_cluster4_param1_11_9_25.csv")

#colnocyte clusters 
c0_1 <- FindMarkers(rectum, ident.1 = "0", ident.2 = "1")

write.csv(c0_1, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster0_vs_cluster1_param1_11_9_25.csv")

c0_5 <- FindMarkers(rectum, ident.1 = "0", ident.2 = "5")

write.csv(c0_5, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster0_vs_cluster5_param1_11_9_25.csv")

#b cells 
c2_26 <- FindMarkers(rectum, ident.1 = "2", ident.2 = "26")

write.csv(c2_26, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster2_vs_cluster26_param1_11_9_25.csv")

c7_26 <- FindMarkers(rectum, ident.1 = "7", ident.2 = "26")

write.csv(c7_26, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/cell_annotation/fm_rectum_cluster7_vs_cluster26_param1_11_9_25.csv")

#--------preliminary annotation of rectum------------#

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/seurat_object/helm_batch1_13_rectum_75_15_2000_rpca_sub_cluster_11_6_24.rds")

rectum$orig.ident <- factor(rectum$orig.ident, levels = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"))
rectum$batch <- factor(rectum$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

rectum@meta.data$cell_typev1 <- rectum@meta.data$seurat_clusters
rectum$cell_typev1 <- plyr::mapvalues(
  x = rectum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27'),
  to = c('Colonocyte 1', 'Colonocyte 2', 'Naive B Cell', 'Stem Cell', 'TA Cell', 'Colonocyte 3', 'Mature Colonocyte 1', 'Memory B Cell', 'CD4 T Cell', 'Goblet Cell 1', 'Plasma Cell 1', 'Macrophage', 'Mature Colonocyte 2', 'Cycling B Cell', 'CD8 T Cell', 'Goblet Cell 2', 'Plasma Cell 2', 'Mature Colonocyte 3', 'Ambig. T Cell', 'NK/T Cell', 'Enteroendocrine Cell', 'Tuft Cell', 'Mesenchymal Cell', 'Monocyte', 'Mast Cell', 'Endothelial Cell', 'Ambig. B Cell', 'ILC3')
)

rectum@meta.data$cell_typev2 <- rectum@meta.data$seurat_clusters
rectum$cell_typev2 <- plyr::mapvalues(
  x = rectum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27'),
  to = c('Colonocyte', 'Colonocyte', 'Naive B Cell', 'Stem Cell', 'TA Cell', 'Colonocyte', 'Mature Colonocyte', 'Memory B Cell', 'CD4 T Cell', 'Goblet Cell', 'Plasma Cell', 'Macrophage', 'Mature Colonocyte', 'Cycling B Cell', 'CD8 T Cell', 'Goblet Cell', 'Plasma Cell', 'Mature Colonocyte', 'Ambig. T Cell', 'NK/CD8 T Cell', 'Enteroendocrine Cell', 'Tuft Cell', 'Mesenchymal Cell', 'Monocyte', 'Mast Cell', 'Endothelial Cell', 'Ambig. B Cell', 'ILC3')
)

rectum$cell_typev1 <- factor(rectum$cell_typev1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell 1", "Plasma Cell 2", "Ambig. B Cell", "CD4 T Cell", "CD8 T Cell", "NK/T Cell", "ILC3", "Ambig. T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem Cell", "TA Cell", "Colonocyte 1", "Colonocyte 2", "Colonocyte 3", "Mature Colonocyte 1", "Mature Colonocyte 2", "Mature Colonocyte 3", "Goblet Cell 1", "Goblet Cell 2", "Tuft Cell", "Enteroendocrine Cell", "Mesenchymal Cell", "Endothelial Cell"))

rectum$cell_typev2 <- factor(rectum$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell", "Ambig. B Cell", "CD4 T Cell", "CD8 T Cell", "NK/CD8 T Cell", "ILC3", "Ambig. T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem Cell", "TA Cell", "Colonocyte", "Mature Colonocyte", "Goblet Cell", "Tuft Cell", "Enteroendocrine Cell", "Mesenchymal Cell", "Endothelial Cell"))

#blue = t cells 
blue <- randomColor(5, hue = c("blue"))
blue <- c("#34bfe5", "#4a93ce", "#b0daf2", "#5376e8", "#1d0fdd") 

#green = b cells
green <- randomColor(6, hue = c("green"))
green <- c("#05fc95", "#9de587", "#1da54f", "#61e291", "#2bd147", "#14ba85")


#orange = myeloid 
orange <- randomColor(4, hue = c("orange"))
orange <- c("#fcab85", "#d3b050", "#f75f13", "#cc5414")

#pink = epithelial
pink <- randomColor(12, hue = c("pink"))
pink <- c("#f41f8d", "#ffb7de", "#ed6add", "#fd8cff", "#cc0076", "#ffbff9", "#f271d6", "#f279b5", "#f963d6", "#fcbde3", "#fc05e7", "#ea62bf")




#purple = mesenchymal and endothelial 
purple <- randomColor(2, hue = c("purple"))
purple <- c("#cd7af9", "#996cc9")

#UMAPs
DimPlot(rectum, reduction = "umap", group.by = "cell_typev1", cols = c(green[1], green[2], green[3], green[4], green[5], green[6], blue[1], blue[2], blue[3], blue[4], blue[5], orange[1], orange[2], orange[3], pink[1], pink[2], pink[3], pink[4], pink[11], pink[6], pink[7], pink[8], pink[9], pink[10], pink[5], pink[12], purple[1], purple[2]), raster = F)

DimPlot(rectum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[3], green[4], green[6], blue[1], blue[2], blue[3], blue[4], blue[5], orange[1], orange[2], orange[3], pink[1], pink[2], pink[3], pink[6], pink[9], pink[5], pink[12], purple[1], purple[2]), raster = F)


#stacked bar plots 

Idents(rectum) <- rectum$cell_typev1
pt2 <- table(Idents(rectum), rectum$orig.ident)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

pt2$Var1 <- factor(pt2$Var1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell 1", "Plasma Cell 2", "Ambig. B Cell", "CD4 T Cell", "CD8 T Cell", "NK/T Cell", "ILC3", "Ambig. T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem Cell", "TA Cell", "Colonocyte 1", "Colonocyte 2", "Colonocyte 3", "Mature Colonocyte 1", "Mature Colonocyte 2", "Mature Colonocyte 3", "Goblet Cell 1", "Goblet Cell 2", "Tuft Cell", "Enteroendocrine Cell", "Mesenchymal Cell", "Endothelial Cell"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[4], green[5], green[6], blue[1], blue[2], blue[3], blue[4], blue[5], orange[1], orange[2], orange[3], pink[1], pink[2], pink[3], pink[4], pink[11], pink[6], pink[7], pink[8], pink[9], pink[10], pink[5], pink[12], purple[1], purple[2])) +
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

Idents(rectum) <- rectum$cell_typev2
pt2 <- table(Idents(rectum), rectum$orig.ident)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

pt2$Var1 <- factor(pt2$Var1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell", "Ambig. B Cell", "CD4 T Cell", "CD8 T Cell", "NK/T Cell", "ILC3", "Ambig. T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem Cell", "TA Cell", "Colonocyte", "Mature Colonocyte", "Goblet Cell", "Tuft Cell", "Enteroendocrine Cell", "Mesenchymal Cell", "Endothelial Cell"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[4], green[6], blue[1], blue[2], blue[3], blue[4], blue[5], orange[1], orange[2], orange[3], pink[1], pink[2], pink[3], pink[6], pink[9], pink[5], pink[12], purple[1], purple[2])) +
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

#------cell proportion test-------#

library(speckle)
library(limma)
library(ggplot2)

# Run propeller testing for cell type proportion differences between the two 
# groups
macro_if <- propeller(clusters = rectum$cell_typev1, sample = rectum$orig.ident, group = rectum$macro_IF)

#not significant

macro_if <- propeller(clusters = rectum$cell_typev2, sample = rectum$orig.ident, group = rectum$macro_IF)

#not signficant

micro_if <- propeller(clusters = rectum$cell_typev1, sample = rectum$orig.ident, group = rectum$micro_IF)

#not significant

micro_if <- propeller(clusters = rectum$cell_typev2, sample = rectum$orig.ident, group = rectum$micro_IF)

#not significant 

batch <- propeller(clusters = rectum$cell_typev1, sample = rectum$orig.ident, group = rectum$batch)

#ambig t, CD8 t cell, NK/T cell - significantly associated with batch (adj. p values = 0.028, 0.028, 0.039)

batch <- propeller(clusters = rectum$cell_typev2, sample = rectum$orig.ident, group = rectum$batch)

#ambig t, CD8 t cell, NK/T cell - significantly associated with batch (adj. p values = 0.028, 0.028, 0.039)

sex <- propeller(clusters = rectum$cell_typev1, sample = rectum$orig.ident, group = rectum$sex)

#mature colonocyte 2 - marginally significant (adj p value = 0.054)

sex <- propeller(clusters = rectum$cell_typev2, sample = rectum$orig.ident, group = rectum$sex)

#not signicant 

ancestry <- propeller(clusters = rectum$cell_typev1, sample = rectum$orig.ident, group = rectum$ancestry)

#not significant 

ancestry <- propeller(clusters = rectum$cell_typev2, sample = rectum$orig.ident, group = rectum$ancestry)

#not significant 

#---------Hierarchical clustering-----------#

props <- getTransformedProps(rectum$cell_typev2, rectum$donor, transform="logit")

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


save_pheatmap_pdf(hmap, "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/rectum/heatmap_batch1_13_rectum_hierarchical_cluster_cellprop_v2_11_25_24.pdf", width = 8, height = 5)

#-------scITD group info--------#

#12/8/2024

#subset donors not in scITD 
Idents(rectum) <- rectum$donor

rectum_sub <- subset(rectum, idents = c("donor3", "donor8", "donor9", "donor24"), invert = T)

rectum_sub@meta.data$group <- rectum_sub@meta.data$donor
rectum_sub$group <- plyr::mapvalues(
  x = rectum_sub$donor,
  from = c("donor1", "donor4", "donor5", "donor6", "donor7", 'donor10', "donor11", "donor12", "donor13", "donor14", "donor15", "donor16", "donor17", "donor18", "donor19", "donor20", "donor21", "donor22", "donor25", "donor26", "donor27", "donor28", "donor29", "donor30", "donor31", "donor32", "donor33", "donor34"),
  to = c("group2", 'group1', 'group2', 'group2', 'group2', 'group1', 'group2', 'group2', 'group2', 'group1', 'group2', 'group1', 'group2', 'group1', 'group2', 'group2', 'group1', 'group1', 'group1', 'group2', 'group2', 'group1', 'group2', 'group1', 'group1', 'group2', 'group1', 'group1')
)

#save the seurat object with group info 
saveRDS(rectum_sub, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_12_8_24.rds")

#convert to sce object on server 

#load in seurat object 
rectum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_12_8_24.rds")

#set assay to RNA
DefaultAssay(rectum) <- "RNA"

#convert to sce object 
rectum.sce <- as.SingleCellExperiment(rectum)

#save 
saveRDS(rectum.sce, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_sce_12_8_24.rds")

#--------prep data for cellchat---------#

#load in scITD group SO
rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_12_8_24.rds")

Idents(rectum) <- rectum$group

rectum_g1 <- subset(rectum, idents = "group1")

saveRDS(rectum_g1, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_g1_12_11_24.rds")

rectum_g2 <- subset(rectum, idents = "group2")

saveRDS(rectum_g2, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_g2_12_11_24.rds")


#-----------gene module score------------#



#use UCell to add TNF and Interferon gamma score 
library(UCell)

ifn_genes <- read.csv(file = "/storage/home/swashburn30/Helmsley/Batch1_13/DEG_analysis/type_ii_interferon_genes.csv", header = FALSE)

rectum <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_12_8_24.rds")

DefaultAssay(rectum) <- "RNA"

rectum$group <- factor(rectum$group, levels = c("group1", "group2"))

rectum_gene <- rownames(rectum)

ifn_rectum <- intersect(ifn_genes$V1, rectum_gene)

rectum <- AddModuleScore_UCell(rectum, features = list(ifn_rectum), name = "ifn_sig")

#TNF response module 

tnf_genes <- read.csv(file = "/storage/home/swashburn30/Helmsley/Batch1_13/DEG_analysis/response_tnf_gobp_genes.csv")

rectum_gene <- rownames(rectum)

tnf_rectum <- intersect(tnf_genes$GENE_SYMBOLS, rectum_gene)

rectum <- AddModuleScore_UCell(rectum, features = list(tnf_rectum), name = "tnf_sig")

#save the rectum object with module score 
saveRDS(rectum, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_rectum_group_scITD_mod_score_12_11_24.rds")

#----------IFN vs. TNF response---------#

#load in rectum with IFN and TNF response 
rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_group_scITD_mod_score_12_11_24.rds")

#re-order donors by group 
rectum$donor <- factor(rectum$donor, levels = c("donor4", "donor10", "donor14", "donor16", "donor18", "donor21", "donor22", "donor25", "donor28", "donor30", "donor31", "donor33", "donor34", "donor1", "donor5", "donor6", "donor7", "donor11", "donor12", "donor13", "donor15", "donor17", "donor19", "donor20", "donor26", "donor27", "donor29", "donor32"))

#group by donor
VlnPlot(rectum, features = "signature_1ifn_sig", group.by = "donor", pt.size = 0) + theme(legend.position = "none") 

VlnPlot(rectum, features = "signature_1tnf_sig", group.by = "donor", pt.size = 0) + theme(legend.position = "none") 

#group by cell type 
VlnPlot(rectum, features = "signature_1ifn_sig", group.by = "cell_typev2", pt.size = 0, split.by = "group", cols = c('#F8766D', '#619CFF'), sort = "decreasing") 

VlnPlot(rectum, features = "signature_1tnf_sig", group.by = "cell_typev2", pt.size = 0,  split.by = "group", cols = c('#F8766D', '#619CFF'), sort = "decreasing") 

#---lollipop chart----#

library(fmsb)
library(ggplot2)

#calculate average score per donor 
rectum_metadata <- rectum@meta.data

# Calculate average module scores for each donor
average_scores <- rectum_metadata %>%
  group_by(donor) %>%
  summarize(across(starts_with("signature"), ~ mean(.x, na.rm = TRUE)))

average_scores$group <- NA

average_scores$group <- c("group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2")

average_scores <- as.data.frame(average_scores)

#-------create lollipop plot for cell types driving score--------#

# Calculate the average gene module score per cell type split by group
average_scores_celltype <- rectum_metadata %>%
  group_by(cell_typev2, group) %>%
  summarize(
    mean_ifn_score = mean(signature_1ifn_sig, na.rm = TRUE),
    mean_tnf_score = mean(signature_1tnf_sig, na.rm = TRUE),
    .groups = "drop"
  )

# Reshape data to long format
df_long_cell <- average_scores_celltype %>%
  pivot_longer(cols = c(mean_ifn_score, mean_tnf_score), names_to = "score_type", values_to = "score")

test <- ggplot(df_long_cell, aes(x = cell_typev2, y = score, color = group, group = cell_typev2)) +
  geom_line(size = 1.2) +  # Connect IFN and TNF scores with a line
  geom_point(size = 4) +   # Add points for scores
  facet_wrap(~score_type, scales = "free_y", labeller = labeller(mean_ifn_score = ifn_sig_name, mean_tnf_score = tnf_sig_name)) + # Optional: Separate by group
  theme_minimal(base_size = 20) +
  labs(x = "Donor", y = "Score", color = "Group") +
  scale_color_manual(labels = c("Group 1", "Group 2"), values = c("#69b3a2", "#dab1da")) +
  coord_flip()

test
#-----overall mean module score per group---------#
# Calculate the mean module scores for each group
overall_mean <- rectum_metadata %>%
  group_by(group) %>%
  summarize(
    mean_ifn = mean(signature_1ifn_sig, na.rm = TRUE),
    mean_tnf = mean(signature_1tnf_sig, na.rm = TRUE),
    .groups = "drop"
  )

#mean IFN: g1 = 0.111 and g2 = 0.102
#mean TNF: g1 = 0.124 and g2 = 0.122


#--------calculate the mean donor score per cell---------#

average_scores_celltype_donor <- rectum_metadata %>%
  group_by(cell_typev2, donor, group) %>%
  summarize(
    mean_ifn_score = mean(signature_1ifn_sig, na.rm = TRUE),
    mean_tnf_score = mean(signature_1tnf_sig, na.rm = TRUE),
    .groups = "drop"
  )

# Reshape into wide format for paired testing

# Reshape into wide format for paired testing
df_wide <- average_scores_celltype_donor %>%
  pivot_wider(names_from = group, values_from = c(mean_ifn_score, mean_tnf_score))

# Run paired t-test for each cell type's IFN score
paired_ttest_ifn <- df_wide %>%
  group_by(cell_typev2, group) %>%
  summarise(
    p_value_ifn = t.test(mean_ifn_score_group1, mean_ifn_score_group2, paired = TRUE)$p.value
  )

results <- average_scores_celltype_donor %>%
  group_by(cell_typev2) %>%
  summarize(
    ifn_p_value = wilcox.test(mean_ifn_score ~ group, data = cur_data())$p.value,
    tnf_p_value = wilcox.test(mean_tnf_score ~ group, data = cur_data())$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    ifn_p_adj = p.adjust(ifn_p_value, method = "fdr"),
    tnf_p_adj = p.adjust(tnf_p_value, method = "fdr")
  )

#strict filters
results_filt_sig_ifn <- filter(results, ifn_p_adj <= 0.05) #CD4 t cell, macrophage, stem cell, TA cell, mature colonocyte
results_filt_sig_tnf <- filter(results, tnf_p_adj <= 0.05) #macrophage

#--------proportion difference between groups--------# 

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_12_8_24.rds")

rectum$group <- factor(rectum$group, levels = c("group1", "group2"))

group <- propeller(clusters = rectum$cell_typev2, sample = rectum$orig.ident, group = rectum$group)

#cell proportion not statistically significant between groups
