#Helmsley Batch 1-13 COLON Analysis

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

#perform on server 

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_clustered_10_17_24.rds")

colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))

#visualize distribution of data (clusters)
pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/VlnPlot_batch1_13_colon_nfeature_rna.pdf", width = 15, height = 8)
plot(VlnPlot(colon, features = 'nFeature_RNA', group.by = "orig.ident", cols = c("helm_sam2" = "#BBCCEE", "helm_sam8" = "#BBCCEE", "helm_sam26" = "#99DDFF", "helm_sam31" = "#99DDFF", "helm_sam33" = "#99DDFF", "helm_sam36" = "#99DDFF", "helm_sam45" = "#44BB99", "helm_sam48" = "#44BB99", "helm_sam52" = "#EE8866", "helm_sam54" = "#EE8866", "helm_sam57" = "#EE8866", "helm_sam60" = "#EE8866", "helm_sam63" = "#BBCC33", "helm_sam65" = "#BBCC33", "helm_sam66" = "#BBCC33", "helm_sam71" = "#BBCC33", "helm_sam72" = "#882255", "helm_sam78" = "#882255", "helm_sam84" = "#882255", "helm_sam86" = "#EECC66", "helm_sam89" = "#EECC66", "helm_sam92" = "#EECC66", "helm_sam96" = "#EECC66", "helm_sam99" = "#EECC66", "helm_sam102" = "#EECC66", "helm_sam110" = "#FE036A", "helm_sam113" = "#FE036A", "helm_sam120" = "#A8CED2", "helm_sam128" = "#A8CED2", "helm_sam129" = "#A8CED2", "helm_sam135" = "#B2DE7C", "helm_sam148" = "#B2DE7C", "helm_sam154" = "#E19390", "helm_sam157" = "#E19390", "helm_sam163" = "#E19390", "helm_sam172" = "#BA6BD4"), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/VlnPlot_batch1_13_colon_ncount_rna.pdf", width = 15, height = 8)
VlnPlot(colon, features = 'nCount_RNA', group.by = "orig.ident", cols = c("helm_sam2" = "#BBCCEE", "helm_sam8" = "#BBCCEE", "helm_sam26" = "#99DDFF", "helm_sam31" = "#99DDFF", "helm_sam33" = "#99DDFF", "helm_sam36" = "#99DDFF", "helm_sam45" = "#44BB99", "helm_sam48" = "#44BB99", "helm_sam52" = "#EE8866", "helm_sam54" = "#EE8866", "helm_sam57" = "#EE8866", "helm_sam60" = "#EE8866", "helm_sam63" = "#BBCC33", "helm_sam65" = "#BBCC33", "helm_sam66" = "#BBCC33", "helm_sam71" = "#BBCC33", "helm_sam72" = "#882255", "helm_sam78" = "#882255", "helm_sam84" = "#882255", "helm_sam86" = "#EECC66", "helm_sam89" = "#EECC66", "helm_sam92" = "#EECC66", "helm_sam96" = "#EECC66", "helm_sam99" = "#EECC66", "helm_sam102" = "#EECC66", "helm_sam110" = "#FE036A", "helm_sam113" = "#FE036A", "helm_sam120" = "#A8CED2", "helm_sam128" = "#A8CED2", "helm_sam129" = "#A8CED2", "helm_sam135" = "#B2DE7C", "helm_sam148" = "#B2DE7C", "helm_sam154" = "#E19390", "helm_sam157" = "#E19390", "helm_sam163" = "#E19390", "helm_sam172" = "#BA6BD4"), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/VlnPlot_batch1_13_colon_percentmt.pdf", width = 15, height = 8)
VlnPlot(colon, features = 'percent.mt', group.by = "orig.ident", cols = c("helm_sam2" = "#BBCCEE", "helm_sam8" = "#BBCCEE", "helm_sam26" = "#99DDFF", "helm_sam31" = "#99DDFF", "helm_sam33" = "#99DDFF", "helm_sam36" = "#99DDFF", "helm_sam45" = "#44BB99", "helm_sam48" = "#44BB99", "helm_sam52" = "#EE8866", "helm_sam54" = "#EE8866", "helm_sam57" = "#EE8866", "helm_sam60" = "#EE8866", "helm_sam63" = "#BBCC33", "helm_sam65" = "#BBCC33", "helm_sam66" = "#BBCC33", "helm_sam71" = "#BBCC33", "helm_sam72" = "#882255", "helm_sam78" = "#882255", "helm_sam84" = "#882255", "helm_sam86" = "#EECC66", "helm_sam89" = "#EECC66", "helm_sam92" = "#EECC66", "helm_sam96" = "#EECC66", "helm_sam99" = "#EECC66", "helm_sam102" = "#EECC66", "helm_sam110" = "#FE036A", "helm_sam113" = "#FE036A", "helm_sam120" = "#A8CED2", "helm_sam128" = "#A8CED2", "helm_sam129" = "#A8CED2", "helm_sam135" = "#B2DE7C", "helm_sam148" = "#B2DE7C", "helm_sam154" = "#E19390", "helm_sam157" = "#E19390", "helm_sam163" = "#E19390", "helm_sam172" = "#BA6BD4"), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


#visualize distribution of data (clusters)
pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/VlnPlot_batch1_13_colon_nfeature_rna_cluster.pdf", width = 12, height = 8)
plot(VlnPlot(colon, features = 'nFeature_RNA', pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/VlnPlot_batch1_13_colon_ncount_rna_cluster.pdf", width = 12, height = 8)
VlnPlot(colon, features = 'nCount_RNA', pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/VlnPlot_batch1_13_colon_percentmt_cluster.pdf", width = 12, height = 8)
VlnPlot(colon, features = 'percent.mt', pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#visualize PCA results
pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/PCA_batch1_13_colon_cluster.pdf", width = 10, height = 8)
plot(DimPlot(colon, reduction = 'pca', group.by = "seurat_clusters"))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/PCA_batch1_13_colon_cluster_sample.pdf", width = 12, height = 8)
plot(DimPlot(colon, reduction = 'pca', group.by = "orig.ident"))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/PCA_batch1_13_colon_cluster_batch.pdf", width = 10, height = 8)
plot(DimPlot(colon, reduction = 'pca', group.by = "batch"))
dev.off()

#visualize UMAP
pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/UMAP_batch1_13_colon_cluster.pdf", width = 10, height = 8)
plot(DimPlot(colon, reduction = 'umap', group.by = "seurat_clusters", label = TRUE))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/UMAP_batch1_13_colon_cluster_sample.pdf", width = 12, height = 8)
plot(DimPlot(colon, reduction = 'umap', group.by = "orig.ident"))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/UMAP_batch1_13_colon_cluster_batch.pdf", width = 10, height = 8)
DimPlot(colon, reduction = 'umap', group.by = "batch")
dev.off()

#proportion of cells per cluster 
colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
Idents(colon) <- colon$seurat_clusters
pt1 <- table(Idents(colon), colon$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

library(randomcoloR)
no_of_colors <- 15

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/StackedBarPlot_batch1_13_colon_cluster_sample.pdf", width = 10, height = 8)
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
dev.off()

#proportion of cells per batch
Idents(colon) <- colon$batch
pt2 <- table(Idents(colon), colon$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))


pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/StackedBarPlot_batch1_13_colon_cluster_batch.pdf", width = 10, height = 8)
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
dev.off()

#----------preliminary rPCA clustering results-----------#

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_rpca_10_17_24.rds")

colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
colon$batch <- factor(colon$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(colon, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(colon, features = 'nCount_RNA', pt.size = 0)
VlnPlot(colon, features = 'percent.mt', pt.size = 0)

#visualize PCA results
DimPlot(colon, reduction = 'pca', group.by = "seurat_clusters", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "batch", raster = F)

#visualize UMAP
DimPlot(colon, reduction = 'umap', group.by = "seurat_clusters", label = TRUE, raster = F)
DimPlot(colon, reduction = 'umap', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'umap', group.by = "batch", raster = F)

#proportion of cells per cluster 
colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
Idents(colon) <- colon$seurat_clusters
pt1 <- table(Idents(colon), colon$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

library(randomcoloR)
no_of_colors <- 16

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
Idents(colon) <- colon$batch
pt2 <- table(Idents(colon), colon$seurat_clusters)
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

#--------------preliminary harmony clustering results------------#

#perform on server
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_harmony_10_17_24.rds")

colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
colon$batch <- factor(colon$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))



#visualize distribution of data (clusters)
pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/VlnPlot_batch1_13_colon_nfeature_rna_harmony.pdf", width = 12, height = 8)
plot(VlnPlot(colon, features = 'nFeature_RNA', pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/VlnPlot_batch1_13_colon_ncount_rna_harmony.pdf", width = 12, height = 8)
VlnPlot(colon, features = 'nCount_RNA', pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/VlnPlot_batch1_13_colon_percentmt_harmony.pdf", width = 12, height = 8)
VlnPlot(colon, features = 'percent.mt', pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#visualize PCA results
pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/PCA_batch1_13_colon_harmony.pdf", width = 10, height = 8)
plot(DimPlot(colon, reduction = 'pca', group.by = "seurat_clusters"))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/PCA_batch1_13_colon_harmony_sample.pdf", width = 12, height = 8)
plot(DimPlot(colon, reduction = 'pca', group.by = "orig.ident"))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/PCA_batch1_13_colon_harmony_batch.pdf", width = 10, height = 8)
plot(DimPlot(colon, reduction = 'pca', group.by = "batch"))
dev.off()

#visualize UMAP
pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/UMAP_batch1_13_colon_harmony.pdf", width = 10, height = 8)
plot(DimPlot(colon, reduction = 'umap', group.by = "seurat_clusters", label = TRUE))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/UMAP_batch1_13_colon_harmony_sample.pdf", width = 12, height = 8)
plot(DimPlot(colon, reduction = 'umap', group.by = "orig.ident"))
dev.off()

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/UMAP_batch1_13_colon_harmony_batch.pdf", width = 10, height = 8)
DimPlot(colon, reduction = 'umap', group.by = "batch")
dev.off()

#proportion of cells per cluster 
colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
Idents(colon) <- colon$seurat_clusters
pt1 <- table(Idents(colon), colon$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"))

library(randomcoloR)
no_of_colors <- 15

# sample colors 
palette <- distinctColorPalette(no_of_colors)      

# hex color codes 
palette 

pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/StackedBarPlot_batch1_13_colon_harmony_sample.pdf", width = 10, height = 8)
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
dev.off()

#proportion of cells per batch
Idents(colon) <- colon$batch
pt2 <- table(Idents(colon), colon$seurat_clusters)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

pt2$Var1 <- factor(pt2$Var1, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))


pdf(file = "/storage/home/swashburn30/Helmsley/Batch1_13/figures/StackedBarPlot_batch1_13_colon_harmony_batch.pdf", width = 10, height = 8)
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
dev.off()

#-------------assess CCA preliminary clustering results------------#

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_cca_10_17_24.rds")

colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
colon$batch <- factor(colon$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(colon, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(colon, features = 'nCount_RNA', pt.size = 0)
VlnPlot(colon, features = 'percent.mt', pt.size = 0)

#visualize PCA results
DimPlot(colon, reduction = 'pca', group.by = "seurat_clusters", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "batch", raster = F)

#visualize UMAP
DimPlot(colon, reduction = 'umap', group.by = "seurat_clusters", label = TRUE, raster = F)
DimPlot(colon, reduction = 'umap', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'umap', group.by = "batch", raster = F)

#proportion of cells per cluster 
colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
Idents(colon) <- colon$seurat_clusters
pt1 <- table(Idents(colon), colon$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

library(randomcoloR)
no_of_colors <- 16

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
Idents(colon) <- colon$batch
pt2 <- table(Idents(colon), colon$seurat_clusters)
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


#--------------Low resolution Cluster Stability Assessment---------------#


#-----Query 1------#

cluster1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub_res20_1_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #13622

#cluster 1
Idents(cluster1) <- cluster1$seurat_clusters
sub1 <- WhichCells(cluster1, idents = c('1'))

Idents(cluster1) <- cluster1$predicted.id
pred1 <- WhichCells(cluster1, idents = '1')

int1 <- intersect(sub1, pred1) #10557

#cluster 2
Idents(cluster1) <- cluster1$seurat_clusters
sub2 <- WhichCells(cluster1, idents = c('3'))

Idents(cluster1) <- cluster1$predicted.id
pred2 <- WhichCells(cluster1, idents = '2')

int2 <- intersect(sub2, pred2) #5422

#cluster 3
Idents(cluster1) <- cluster1$seurat_clusters
sub3 <- WhichCells(cluster1, idents = c('2'))

Idents(cluster1) <- cluster1$predicted.id
pred3 <- WhichCells(cluster1, idents = '3')

int3 <- intersect(sub3, pred3) #5697

#cluster 4
Idents(cluster1) <- cluster1$seurat_clusters
sub4 <- WhichCells(cluster1, idents = c('5'))

Idents(cluster1) <- cluster1$predicted.id
pred4 <- WhichCells(cluster1, idents = '4')

int4 <- intersect(sub4, pred4) #4290

#cluster 5
Idents(cluster1) <- cluster1$seurat_clusters
sub5 <- WhichCells(cluster1, idents = c('4'))

Idents(cluster1) <- cluster1$predicted.id
pred5 <- WhichCells(cluster1, idents = '5')

int5 <- intersect(sub5, pred5) #3655

#cluster 6 
Idents(cluster1) <- cluster1$seurat_clusters
sub6 <- WhichCells(cluster1, idents = c('6'))

Idents(cluster1) <- cluster1$predicted.id
pred6 <- WhichCells(cluster1, idents = '6')

int6 <- intersect(sub6, pred6) #2483

#cluster 7 
Idents(cluster1) <- cluster1$seurat_clusters
sub7 <- WhichCells(cluster1, idents = c('7'))

Idents(cluster1) <- cluster1$predicted.id
pred7 <- WhichCells(cluster1, idents = '7')

int7 <- intersect(sub7, pred7) #2317

#cluster 8 
Idents(cluster1) <- cluster1$seurat_clusters
sub8 <- WhichCells(cluster1, idents = c('8'))

Idents(cluster1) <- cluster1$predicted.id
pred8 <- WhichCells(cluster1, idents = '8')

int8 <- intersect(sub8, pred8) #1898

#cluster 9 
Idents(cluster1) <- cluster1$seurat_clusters
sub9 <- WhichCells(cluster1, idents = c('9'))

Idents(cluster1) <- cluster1$predicted.id
pred9 <- WhichCells(cluster1, idents = '9')

int9 <- intersect(sub9, pred9) #915

#cluster 10 
Idents(cluster1) <- cluster1$seurat_clusters
sub10 <- WhichCells(cluster1, idents = c('11'))

Idents(cluster1) <- cluster1$predicted.id
pred10 <- WhichCells(cluster1, idents = '10')

int10 <- intersect(sub10, pred10) #866

#cluster 11
Idents(cluster1) <- cluster1$seurat_clusters
sub11 <- WhichCells(cluster1, idents = c('10'))

Idents(cluster1) <- cluster1$predicted.id
pred11 <- WhichCells(cluster1, idents = '11')

int11 <- intersect(sub11, pred11) #836

#cluster 12 
Idents(cluster1) <- cluster1$seurat_clusters
sub12 <- WhichCells(cluster1, idents = c('12'))

Idents(cluster1) <- cluster1$predicted.id
pred12 <- WhichCells(cluster1, idents = '12')

int12 <- intersect(sub12, pred12) #434

#cluster 13
Idents(cluster1) <- cluster1$seurat_clusters
sub13 <- WhichCells(cluster1, idents = c('13'))

Idents(cluster1) <- cluster1$predicted.id
pred13 <- WhichCells(cluster1, idents = '13')

int13 <- intersect(sub13, pred13) #350

#cluster 14
Idents(cluster1) <- cluster1$seurat_clusters
sub14 <- WhichCells(cluster1, idents = c('14'))

Idents(cluster1) <- cluster1$predicted.id
pred14 <- WhichCells(cluster1, idents = '14')

int14 <- intersect(sub14, pred14) #251

#cluster 15
Idents(cluster1) <- cluster1$seurat_clusters
sub15 <- WhichCells(cluster1, idents = c('15'))

Idents(cluster1) <- cluster1$predicted.id
pred15 <- WhichCells(cluster1, idents = '15')

int15 <- intersect(sub15, pred15) #59

#total
#pre:55641
#53652
total1 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total1,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query1_low_res_barcodes_10_28_24.csv', row.names=F)

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

cluster2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub2_res20_1_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #14246

#cluster 1
Idents(cluster2) <- cluster2$seurat_clusters
sub1 <- WhichCells(cluster2, idents = c('1'))

Idents(cluster2) <- cluster2$predicted.id
pred1 <- WhichCells(cluster2, idents = '1')

int1 <- intersect(sub1, pred1) #10517

#cluster 2 
Idents(cluster2) <- cluster2$seurat_clusters
sub2 <- WhichCells(cluster2, idents = c('2'))

Idents(cluster2) <- cluster2$predicted.id
pred2 <- WhichCells(cluster2, idents = '2')

int2 <- intersect(sub2, pred2) #5588

#cluster 3
Idents(cluster2) <- cluster2$seurat_clusters
sub3 <- WhichCells(cluster2, idents = c('3'))

Idents(cluster2) <- cluster2$predicted.id
pred3 <- WhichCells(cluster2, idents = '3')

int3 <- intersect(sub3, pred3) #5648

#cluster 4
Idents(cluster2) <- cluster2$seurat_clusters
sub4 <- WhichCells(cluster2, idents = c('4'))

Idents(cluster2) <- cluster2$predicted.id
pred4 <- WhichCells(cluster2, idents = '4')

int4 <- intersect(sub4, pred4) #4296

#cluster 5
Idents(cluster2) <- cluster2$seurat_clusters
sub5 <- WhichCells(cluster2, idents = c('5'))

Idents(cluster2) <- cluster2$predicted.id
pred5 <- WhichCells(cluster2, idents = '5')

int5 <- intersect(sub5, pred5) #3598

#cluster 6
Idents(cluster2) <- cluster2$seurat_clusters
sub6 <- WhichCells(cluster2, idents = c('6'))

Idents(cluster2) <- cluster2$predicted.id
pred6 <- WhichCells(cluster2, idents = '6')

int6 <- intersect(sub6, pred6) #2438

#cluster 7 
Idents(cluster2) <- cluster2$seurat_clusters
sub7 <- WhichCells(cluster2, idents = c('7'))

Idents(cluster2) <- cluster2$predicted.id
pred7 <- WhichCells(cluster2, idents = '7')

int7 <- intersect(sub7, pred7) #2388

#cluster 8 
Idents(cluster2) <- cluster2$seurat_clusters
sub8 <- WhichCells(cluster2, idents = c('8'))

Idents(cluster2) <- cluster2$predicted.id
pred8 <- WhichCells(cluster2, idents = '8')

int8 <- intersect(sub8, pred8) #1888

#cluster 9 
Idents(cluster2) <- cluster2$seurat_clusters
sub9 <- WhichCells(cluster2, idents = c('9'))

Idents(cluster2) <- cluster2$predicted.id
pred9 <- WhichCells(cluster2, idents = '9')

int9 <- intersect(sub9, pred9) #936

#cluster 10 
Idents(cluster2) <- cluster2$seurat_clusters
sub10 <- WhichCells(cluster2, idents = c('10'))

Idents(cluster2) <- cluster2$predicted.id
pred10 <- WhichCells(cluster2, idents = '10')

int10 <- intersect(sub10, pred10) #852

#cluster 11
Idents(cluster2) <- cluster2$seurat_clusters
sub11 <- WhichCells(cluster2, idents = c('11'))

Idents(cluster2) <- cluster2$predicted.id
pred11 <- WhichCells(cluster2, idents = '11')

int11 <- intersect(sub11, pred11) #809

#cluster 12
Idents(cluster2) <- cluster2$seurat_clusters
sub12 <- WhichCells(cluster2, idents = c('12'))

Idents(cluster2) <- cluster2$predicted.id
pred12 <- WhichCells(cluster2, idents = '12')

int12 <- intersect(sub12, pred12) #407

#cluster 13
Idents(cluster2) <- cluster2$seurat_clusters
sub13 <- WhichCells(cluster2, idents = c('13'))

Idents(cluster2) <- cluster2$predicted.id
pred13 <- WhichCells(cluster2, idents = '13')

int13 <- intersect(sub13, pred13) #334

#cluster 14
Idents(cluster2) <- cluster2$seurat_clusters
sub14 <- WhichCells(cluster2, idents = c('14'))

Idents(cluster2) <- cluster2$predicted.id
pred14 <- WhichCells(cluster2, idents = '14')

int14 <- intersect(sub14, pred14) #206

#cluster 15
Idents(cluster2) <- cluster2$seurat_clusters
sub15 <- WhichCells(cluster2, idents = c('4'))

Idents(cluster2) <- cluster2$predicted.id
pred15 <- WhichCells(cluster2, idents = '15')

int15 <- intersect(sub15, pred15) #69

#total 
#pre:55641
#54220
total2 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total2,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query2_low_res_barcodes_10_28_24.csv', row.names=F)

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

cluster3 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub_res20_2_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #13764

#cluster 1
Idents(cluster3) <- cluster3$seurat_clusters
sub1 <- WhichCells(cluster3, idents = c('1'))

Idents(cluster3) <- cluster3$predicted.id
pred1 <- WhichCells(cluster3, idents = '1')

int1 <- intersect(sub1, pred1) #10457

#cluster 2
Idents(cluster3) <- cluster3$seurat_clusters
sub2 <- WhichCells(cluster3, idents = c('2'))

Idents(cluster3) <- cluster3$predicted.id
pred2 <- WhichCells(cluster3, idents = '2')

int2 <- intersect(sub2, pred2) #5433

#cluster 3
Idents(cluster3) <- cluster3$seurat_clusters
sub3 <- WhichCells(cluster3, idents = c('3'))

Idents(cluster3) <- cluster3$predicted.id
pred3 <- WhichCells(cluster3, idents = '3')

int3 <- intersect(sub3, pred3) #5742

#cluster 4
Idents(cluster3) <- cluster3$seurat_clusters
sub4 <- WhichCells(cluster3, idents = c('4'))

Idents(cluster3) <- cluster3$predicted.id
pred4 <- WhichCells(cluster3, idents = '4')

int4 <- intersect(sub4, pred4) #4367

#cluster 5
Idents(cluster3) <- cluster3$seurat_clusters
sub5 <- WhichCells(cluster3, idents = c('5'))

Idents(cluster3) <- cluster3$predicted.id
pred5 <- WhichCells(cluster3, idents = '5')

int5 <- intersect(sub5, pred5) #3248

#cluster 6 
Idents(cluster3) <- cluster3$seurat_clusters
sub6 <- WhichCells(cluster3, idents = c('6'))

Idents(cluster3) <- cluster3$predicted.id
pred6 <- WhichCells(cluster3, idents = '6')

int6 <- intersect(sub6, pred6) #2479

#cluster 7 
Idents(cluster3) <- cluster3$seurat_clusters
sub7 <- WhichCells(cluster3, idents = c('7'))

Idents(cluster3) <- cluster3$predicted.id
pred7 <- WhichCells(cluster3, idents = '7')

int7 <- intersect(sub7, pred7) #2297

#cluster 8 
Idents(cluster3) <- cluster3$seurat_clusters
sub8 <- WhichCells(cluster3, idents = c('8'))

Idents(cluster3) <- cluster3$predicted.id
pred8 <- WhichCells(cluster3, idents = '8')

int8 <- intersect(sub8, pred8) #1792

#cluster 9 
Idents(cluster3) <- cluster3$seurat_clusters
sub9 <- WhichCells(cluster3, idents = c('9'))

Idents(cluster3) <- cluster3$predicted.id
pred9 <- WhichCells(cluster3, idents = '9')

int9 <- intersect(sub9, pred9) #883

#cluster 10 
Idents(cluster3) <- cluster3$seurat_clusters
sub10 <- WhichCells(cluster3, idents = c('10'))

Idents(cluster3) <- cluster3$predicted.id
pred10 <- WhichCells(cluster3, idents = '10')

int10 <- intersect(sub10, pred10) #872

#cluster 11
Idents(cluster3) <- cluster3$seurat_clusters
sub11 <- WhichCells(cluster3, idents = c('11'))

Idents(cluster3) <- cluster3$predicted.id
pred11 <- WhichCells(cluster3, idents = '11')

int11 <- intersect(sub11, pred11) #671

#cluster 12
Idents(cluster3) <- cluster3$seurat_clusters
sub12 <- WhichCells(cluster3, idents = c('12'))

Idents(cluster3) <- cluster3$predicted.id
pred12 <- WhichCells(cluster3, idents = '12')

int12 <- intersect(sub12, pred12) #408

#cluster 13
Idents(cluster3) <- cluster3$seurat_clusters
sub13 <- WhichCells(cluster3, idents = c('13'))

Idents(cluster3) <- cluster3$predicted.id
pred13 <- WhichCells(cluster3, idents = '13')

int13 <- intersect(sub13, pred13) #329

#cluster 14
Idents(cluster3) <- cluster3$seurat_clusters
sub14 <- WhichCells(cluster3, idents = c('14'))

Idents(cluster3) <- cluster3$predicted.id
pred14 <- WhichCells(cluster3, idents = '14')

int14 <- intersect(sub14, pred14) #233

#cluster 15
Idents(cluster3) <- cluster3$seurat_clusters
sub15 <- WhichCells(cluster3, idents = c('4'))

Idents(cluster3) <- cluster3$predicted.id
pred15 <- WhichCells(cluster3, idents = '15')

int15 <- intersect(sub15, pred15) #68

#total 
#pre:55641
#53043
total3 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total3,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query3_low_res_barcodes_10_28_24.csv', row.names=F)

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

cluster4 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub2_res20_2_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #12287

#cluster 1
Idents(cluster4) <- cluster4$seurat_clusters
sub1 <- WhichCells(cluster4, idents = c('1'))

Idents(cluster4) <- cluster4$predicted.id
pred1 <- WhichCells(cluster4, idents = '1')

int1 <- intersect(sub1, pred1) #10608

#cluster 2
Idents(cluster4) <- cluster4$seurat_clusters
sub2 <- WhichCells(cluster4, idents = c('2'))

Idents(cluster4) <- cluster4$predicted.id
pred2 <- WhichCells(cluster4, idents = '2')

int2 <- intersect(sub2, pred2) #5661

#cluster 3
Idents(cluster4) <- cluster4$seurat_clusters
sub3 <- WhichCells(cluster4, idents = c('3'))

Idents(cluster4) <- cluster4$predicted.id
pred3 <- WhichCells(cluster4, idents = '3')

int3 <- intersect(sub3, pred3) #5608

#cluster 4
Idents(cluster4) <- cluster4$seurat_clusters
sub4 <- WhichCells(cluster4, idents = c('4'))

Idents(cluster4) <- cluster4$predicted.id
pred4 <- WhichCells(cluster4, idents = '4')

int4 <- intersect(sub4, pred4) #3788

#cluster 5
Idents(cluster4) <- cluster4$seurat_clusters
sub5 <- WhichCells(cluster4, idents = c('5'))

Idents(cluster4) <- cluster4$predicted.id
pred5 <- WhichCells(cluster4, idents = '5')

int5 <- intersect(sub5, pred5) #3471

#cluster 6
Idents(cluster4) <- cluster4$seurat_clusters
sub6 <- WhichCells(cluster4, idents = c('6'))

Idents(cluster4) <- cluster4$predicted.id
pred6 <- WhichCells(cluster4, idents = '6')

int6 <- intersect(sub6, pred6) #2476

#cluster 7 
Idents(cluster4) <- cluster4$seurat_clusters
sub7 <- WhichCells(cluster4, idents = c('4', '8'))

Idents(cluster4) <- cluster4$predicted.id
pred7 <- WhichCells(cluster4, idents = '7')

int7 <- intersect(sub7, pred7) #2451

#cluster 8 
Idents(cluster4) <- cluster4$seurat_clusters
sub8 <- WhichCells(cluster4, idents = c('7'))

Idents(cluster4) <- cluster4$predicted.id
pred8 <- WhichCells(cluster4, idents = '8')

int8 <- intersect(sub8, pred8) #1960

#cluster 9 
Idents(cluster4) <- cluster4$seurat_clusters
sub9 <- WhichCells(cluster4, idents = c('9'))

Idents(cluster4) <- cluster4$predicted.id
pred9 <- WhichCells(cluster4, idents = '9')

int9 <- intersect(sub9, pred9) #974

#cluster 10 
Idents(cluster4) <- cluster4$seurat_clusters
sub10 <- WhichCells(cluster4, idents = c('10'))

Idents(cluster4) <- cluster4$predicted.id
pred10 <- WhichCells(cluster4, idents = '10')

int10 <- intersect(sub10, pred10) #857

#cluster 11
Idents(cluster4) <- cluster4$seurat_clusters
sub11 <- WhichCells(cluster4, idents = c('11'))

Idents(cluster4) <- cluster4$predicted.id
pred11 <- WhichCells(cluster4, idents = '11')

int11 <- intersect(sub11, pred11) #810

#cluster 12
Idents(cluster4) <- cluster4$seurat_clusters
sub12 <- WhichCells(cluster4, idents = c('13'))

Idents(cluster4) <- cluster4$predicted.id
pred12 <- WhichCells(cluster4, idents = '12')

int12 <- intersect(sub12, pred12) #419

#cluster 13
Idents(cluster4) <- cluster4$seurat_clusters
sub13 <- WhichCells(cluster4, idents = c('14'))

Idents(cluster4) <- cluster4$predicted.id
pred13 <- WhichCells(cluster4, idents = '13')

int13 <- intersect(sub13, pred13) #357

#cluster 14
Idents(cluster4) <- cluster4$seurat_clusters
sub14 <- WhichCells(cluster4, idents = c('15'))

Idents(cluster4) <- cluster4$predicted.id
pred14 <- WhichCells(cluster4, idents = '14')

int14 <- intersect(sub14, pred14) #240

#cluster 15
Idents(cluster4) <- cluster4$seurat_clusters
sub15 <- WhichCells(cluster4, idents = c('16'))

Idents(cluster4) <- cluster4$predicted.id
pred15 <- WhichCells(cluster4, idents = '15')

int15 <- intersect(sub15, pred15) #55

#total 
#pre:55641
#51986
total4 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total4,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query4_low_res_barcodes_10_28_24.csv', row.names=F)

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

cluster5 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub_res20_3_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #13835

#cluster 1
Idents(cluster5) <- cluster5$seurat_clusters
sub1 <- WhichCells(cluster5, idents = c('1'))

Idents(cluster5) <- cluster5$predicted.id
pred1 <- WhichCells(cluster5, idents = '1')

int1 <- intersect(sub1, pred1) #10492

#cluster 2
Idents(cluster5) <- cluster5$seurat_clusters
sub2 <- WhichCells(cluster5, idents = c('3'))

Idents(cluster5) <- cluster5$predicted.id
pred2 <- WhichCells(cluster5, idents = '2')

int2 <- intersect(sub2, pred2) #5238

#cluster 3
Idents(cluster5) <- cluster5$seurat_clusters
sub3 <- WhichCells(cluster5, idents = c('2'))

Idents(cluster5) <- cluster5$predicted.id
pred3 <- WhichCells(cluster5, idents = '3')

int3 <- intersect(sub3, pred3) #5677

#cluster 4
Idents(cluster5) <- cluster5$seurat_clusters
sub4 <- WhichCells(cluster5, idents = c('5'))

Idents(cluster5) <- cluster5$predicted.id
pred4 <- WhichCells(cluster5, idents = '4')

int4 <- intersect(sub4, pred4) #4285

#cluster 5
Idents(cluster5) <- cluster5$seurat_clusters
sub5 <- WhichCells(cluster5, idents = c('4'))

Idents(cluster5) <- cluster5$predicted.id
pred5 <- WhichCells(cluster5, idents = '5')

int5 <- intersect(sub5, pred5) #3535

#cluster 6
Idents(cluster5) <- cluster5$seurat_clusters
sub6 <- WhichCells(cluster5, idents = c('7'))

Idents(cluster5) <- cluster5$predicted.id
pred6 <- WhichCells(cluster5, idents = '6')

int6 <- intersect(sub6, pred6) #2448

#cluster 7 
Idents(cluster5) <- cluster5$seurat_clusters
sub7 <- WhichCells(cluster5, idents = c('6'))

Idents(cluster5) <- cluster5$predicted.id
pred7 <- WhichCells(cluster5, idents = '7')

int7 <- intersect(sub7, pred7) #2411

#cluster 8 
Idents(cluster5) <- cluster5$seurat_clusters
sub8 <- WhichCells(cluster5, idents = c('8'))

Idents(cluster5) <- cluster5$predicted.id
pred8 <- WhichCells(cluster5, idents = '8')

int8 <- intersect(sub8, pred8) #1921

#cluster 9 
Idents(cluster5) <- cluster5$seurat_clusters
sub9 <- WhichCells(cluster5, idents = c('9'))

Idents(cluster5) <- cluster5$predicted.id
pred9 <- WhichCells(cluster5, idents = '9')

int9 <- intersect(sub9, pred9) #940

#cluster 10 
Idents(cluster5) <- cluster5$seurat_clusters
sub10 <- WhichCells(cluster5, idents = c('10'))

Idents(cluster5) <- cluster5$predicted.id
pred10 <- WhichCells(cluster5, idents = '10')

int10 <- intersect(sub10, pred10) #871

#cluster 11
Idents(cluster5) <- cluster5$seurat_clusters
sub11 <- WhichCells(cluster5, idents = c('11'))

Idents(cluster5) <- cluster5$predicted.id
pred11 <- WhichCells(cluster5, idents = '11')

int11 <- intersect(sub11, pred11) #783

#cluster 12
Idents(cluster5) <- cluster5$seurat_clusters
sub12 <- WhichCells(cluster5, idents = c('12'))

Idents(cluster5) <- cluster5$predicted.id
pred12 <- WhichCells(cluster5, idents = '12')

int12 <- intersect(sub12, pred12) #422

#cluster 13
Idents(cluster5) <- cluster5$seurat_clusters
sub13 <- WhichCells(cluster5, idents = c('13'))

Idents(cluster5) <- cluster5$predicted.id
pred13 <- WhichCells(cluster5, idents = '13')

int13 <- intersect(sub13, pred13) #343

#cluster 14
Idents(cluster5) <- cluster5$seurat_clusters
sub14 <- WhichCells(cluster5, idents = c('14'))

Idents(cluster5) <- cluster5$predicted.id
pred14 <- WhichCells(cluster5, idents = '14')

int14 <- intersect(sub14, pred14) #243

#cluster 15
Idents(cluster5) <- cluster5$seurat_clusters
sub15 <- WhichCells(cluster5, idents = c('6'))

Idents(cluster5) <- cluster5$predicted.id
pred15 <- WhichCells(cluster5, idents = '15')

int15 <- intersect(sub15, pred15) #60

#total 
#pre:55641
#53504
total5 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total5,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query5_low_res_barcodes_10_28_24.csv', row.names=F)

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

cluster6 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub2_res20_3_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #13439

#cluster 1
Idents(cluster6) <- cluster6$seurat_clusters
sub1 <- WhichCells(cluster6, idents = c('1'))

Idents(cluster6) <- cluster6$predicted.id
pred1 <- WhichCells(cluster6, idents = '1')

int1 <- intersect(sub1, pred1) #10594

#cluster 2
Idents(cluster6) <- cluster6$seurat_clusters
sub2 <- WhichCells(cluster6, idents = c('2'))

Idents(cluster6) <- cluster6$predicted.id
pred2 <- WhichCells(cluster6, idents = '2')

int2 <- intersect(sub2, pred2) #5704

#cluster 3
Idents(cluster6) <- cluster6$seurat_clusters
sub3 <- WhichCells(cluster6, idents = c('3'))

Idents(cluster6) <- cluster6$predicted.id
pred3 <- WhichCells(cluster6, idents = '3')

int3 <- intersect(sub3, pred3) #5674

#cluster 4
Idents(cluster6) <- cluster6$seurat_clusters
sub4 <- WhichCells(cluster6, idents = c('4'))

Idents(cluster6) <- cluster6$predicted.id
pred4 <- WhichCells(cluster6, idents = '4')

int4 <- intersect(sub4, pred4) #4289

#cluster 5
Idents(cluster6) <- cluster6$seurat_clusters
sub5 <- WhichCells(cluster6, idents = c('5'))

Idents(cluster6) <- cluster6$predicted.id
pred5 <- WhichCells(cluster6, idents = '5')

int5 <- intersect(sub5, pred5) #3651

#cluster 6 
Idents(cluster6) <- cluster6$seurat_clusters
sub6 <- WhichCells(cluster6, idents = c('6'))

Idents(cluster6) <- cluster6$predicted.id
pred6 <- WhichCells(cluster6, idents = '6')

int6 <- intersect(sub6, pred6) #2464

#cluster 7 
Idents(cluster6) <- cluster6$seurat_clusters
sub7 <- WhichCells(cluster6, idents = c('4', '8'))

Idents(cluster6) <- cluster6$predicted.id
pred7 <- WhichCells(cluster6, idents = '7')

int7 <- intersect(sub7, pred7) #2280

#cluster 8 
Idents(cluster6) <- cluster6$seurat_clusters
sub8 <- WhichCells(cluster6, idents = c('7'))

Idents(cluster6) <- cluster6$predicted.id
pred8 <- WhichCells(cluster6, idents = '8')

int8 <- intersect(sub8, pred8) #1839

#cluster 9 
Idents(cluster6) <- cluster6$seurat_clusters
sub9 <- WhichCells(cluster6, idents = c('9'))

Idents(cluster6) <- cluster6$predicted.id
pred9 <- WhichCells(cluster6, idents = '9')

int9 <- intersect(sub9, pred9) #919

#cluster 10
Idents(cluster6) <- cluster6$seurat_clusters
sub10 <- WhichCells(cluster6, idents = c('10'))

Idents(cluster6) <- cluster6$predicted.id
pred10 <- WhichCells(cluster6, idents = '10')

int10 <- intersect(sub10, pred10) #847

#cluster 11
Idents(cluster6) <- cluster6$seurat_clusters
sub11 <- WhichCells(cluster6, idents = c('11'))

Idents(cluster6) <- cluster6$predicted.id
pred11 <- WhichCells(cluster6, idents = '11')

int11 <- intersect(sub11, pred11) #673

#cluster 12
Idents(cluster6) <- cluster6$seurat_clusters
sub12 <- WhichCells(cluster6, idents = c('12'))

Idents(cluster6) <- cluster6$predicted.id
pred12 <- WhichCells(cluster6, idents = '12')

int12 <- intersect(sub12, pred12) #413

#cluster 13
Idents(cluster6) <- cluster6$seurat_clusters
sub13 <- WhichCells(cluster6, idents = c('13'))

Idents(cluster6) <- cluster6$predicted.id
pred13 <- WhichCells(cluster6, idents = '13')

int13 <- intersect(sub13, pred13) #340

#cluster 14
Idents(cluster6) <- cluster6$seurat_clusters
sub14 <- WhichCells(cluster6, idents = c('14'))

Idents(cluster6) <- cluster6$predicted.id
pred14 <- WhichCells(cluster6, idents = '14')

int14 <- intersect(sub14, pred14) #228

#cluster 15
Idents(cluster6) <- cluster6$seurat_clusters
sub15 <- WhichCells(cluster6, idents = c('15'))

Idents(cluster6) <- cluster6$predicted.id
pred15 <- WhichCells(cluster6, idents = '15')

int15 <- intersect(sub15, pred15) #64

#total 
#pre:55641
#53418
total6 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total6,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query6_low_res_barcodes_10_28_24.csv', row.names=F)

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

cluster7 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub_res20_4_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #14306

#cluster 1
Idents(cluster7) <- cluster7$seurat_clusters
sub1 <- WhichCells(cluster7, idents = c('1'))

Idents(cluster7) <- cluster7$predicted.id
pred1 <- WhichCells(cluster7, idents = '1')

int1 <- intersect(sub1, pred1) #10565

#cluster 2
Idents(cluster7) <- cluster7$seurat_clusters
sub2 <- WhichCells(cluster7, idents = c('2'))

Idents(cluster7) <- cluster7$predicted.id
pred2 <- WhichCells(cluster7, idents = '2')

int2 <- intersect(sub2, pred2) #5364

#cluster 3
Idents(cluster7) <- cluster7$seurat_clusters
sub3 <- WhichCells(cluster7, idents = c('3'))

Idents(cluster7) <- cluster7$predicted.id
pred3 <- WhichCells(cluster7, idents = '3')

int3 <- intersect(sub3, pred3) #5607

#cluster 4
Idents(cluster7) <- cluster7$seurat_clusters
sub4 <- WhichCells(cluster7, idents = c('4'))

Idents(cluster7) <- cluster7$predicted.id
pred4 <- WhichCells(cluster7, idents = '4')

int4 <- intersect(sub4, pred4) #4202

#cluster 5
Idents(cluster7) <- cluster7$seurat_clusters
sub5 <- WhichCells(cluster7, idents = c('5'))

Idents(cluster7) <- cluster7$predicted.id
pred5 <- WhichCells(cluster7, idents = '5')

int5 <- intersect(sub5, pred5) #3419

#cluster 6 
Idents(cluster7) <- cluster7$seurat_clusters
sub6 <- WhichCells(cluster7, idents = c('6'))

Idents(cluster7) <- cluster7$predicted.id
pred6 <- WhichCells(cluster7, idents = '6')

int6 <- intersect(sub6, pred6) #2465

#cluster 7 
Idents(cluster7) <- cluster7$seurat_clusters
sub7 <- WhichCells(cluster7, idents = c('7'))

Idents(cluster7) <- cluster7$predicted.id
pred7 <- WhichCells(cluster7, idents = '7')

int7 <- intersect(sub7, pred7) #2346

#cluster 8 
Idents(cluster7) <- cluster7$seurat_clusters
sub8 <- WhichCells(cluster7, idents = c('8'))

Idents(cluster7) <- cluster7$predicted.id
pred8 <- WhichCells(cluster7, idents = '8')

int8 <- intersect(sub8, pred8) #1883

#cluster 9 
Idents(cluster7) <- cluster7$seurat_clusters
sub9 <- WhichCells(cluster7, idents = c('9'))

Idents(cluster7) <- cluster7$predicted.id
pred9 <- WhichCells(cluster7, idents = '9')

int9 <- intersect(sub9, pred9) #924

#cluster 10 
Idents(cluster7) <- cluster7$seurat_clusters
sub10 <- WhichCells(cluster7, idents = c('10'))

Idents(cluster7) <- cluster7$predicted.id
pred10 <- WhichCells(cluster7, idents = '10')

int10 <- intersect(sub10, pred10) #850

#cluster 11
Idents(cluster7) <- cluster7$seurat_clusters
sub11 <- WhichCells(cluster7, idents = c('11'))

Idents(cluster7) <- cluster7$predicted.id
pred11 <- WhichCells(cluster7, idents = '11')

int11 <- intersect(sub11, pred11) #804

#cluster 12
Idents(cluster7) <- cluster7$seurat_clusters
sub12 <- WhichCells(cluster7, idents = c('12'))

Idents(cluster7) <- cluster7$predicted.id
pred12 <- WhichCells(cluster7, idents = '12')

int12 <- intersect(sub12, pred12) #428

#cluster 13
Idents(cluster7) <- cluster7$seurat_clusters
sub13 <- WhichCells(cluster7, idents = c('13'))

Idents(cluster7) <- cluster7$predicted.id
pred13 <- WhichCells(cluster7, idents = '13')

int13 <- intersect(sub13, pred13) #362

#cluster 14
Idents(cluster7) <- cluster7$seurat_clusters
sub14 <- WhichCells(cluster7, idents = c('14'))

Idents(cluster7) <- cluster7$predicted.id
pred14 <- WhichCells(cluster7, idents = '14')

int14 <- intersect(sub14, pred14) #224

#cluster 15
Idents(cluster7) <- cluster7$seurat_clusters
sub15 <- WhichCells(cluster7, idents = c('4'))

Idents(cluster7) <- cluster7$predicted.id
pred15 <- WhichCells(cluster7, idents = '15')

int15 <- intersect(sub15, pred15) #65

#total 
#pre:55641
#53814
total7 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total7,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query7_low_res_barcodes_10_28_24.csv', row.names=F)

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

cluster8 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub2_res20_4_low_res_cluster_pred.rds")

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

int0 <- intersect(sub0, pred0) #11681

#cluster 1
Idents(cluster8) <- cluster8$seurat_clusters
sub1 <- WhichCells(cluster8, idents = c('1'))

Idents(cluster8) <- cluster8$predicted.id
pred1 <- WhichCells(cluster8, idents = '1')

int1 <- intersect(sub1, pred1) #10313

#cluster 2
Idents(cluster8) <- cluster8$seurat_clusters
sub2 <- WhichCells(cluster8, idents = c('0', '5'))

Idents(cluster8) <- cluster8$predicted.id
pred2 <- WhichCells(cluster8, idents = '2')

int2 <- intersect(sub2, pred2) #7269

#cluster 3 
Idents(cluster8) <- cluster8$seurat_clusters
sub3 <- WhichCells(cluster8, idents = c('2'))

Idents(cluster8) <- cluster8$predicted.id
pred3 <- WhichCells(cluster8, idents = '3')

int3 <- intersect(sub3, pred3) #5720

#cluster 4
Idents(cluster8) <- cluster8$seurat_clusters
sub4 <- WhichCells(cluster8, idents = c('3'))

Idents(cluster8) <- cluster8$predicted.id
pred4 <- WhichCells(cluster8, idents = '4')

int4 <- intersect(sub4, pred4) #4149

#cluster 5
Idents(cluster8) <- cluster8$seurat_clusters
sub5 <- WhichCells(cluster8, idents = c('4'))

Idents(cluster8) <- cluster8$predicted.id
pred5 <- WhichCells(cluster8, idents = '5')

int5 <- intersect(sub5, pred5) #3488

#cluster 6
Idents(cluster8) <- cluster8$seurat_clusters
sub6 <- WhichCells(cluster8, idents = c('6'))

Idents(cluster8) <- cluster8$predicted.id
pred6 <- WhichCells(cluster8, idents = '6')

int6 <- intersect(sub6, pred6) #2452

#cluster 7 
Idents(cluster8) <- cluster8$seurat_clusters
sub7 <- WhichCells(cluster8, idents = c('3', '8'))

Idents(cluster8) <- cluster8$predicted.id
pred7 <- WhichCells(cluster8, idents = '7')

int7 <- intersect(sub7, pred7) #2412

#cluster 8 
Idents(cluster8) <- cluster8$seurat_clusters
sub8 <- WhichCells(cluster8, idents = c('7'))

Idents(cluster8) <- cluster8$predicted.id
pred8 <- WhichCells(cluster8, idents = '8')

int8 <- intersect(sub8, pred8) #1916

#cluster 9 
Idents(cluster8) <- cluster8$seurat_clusters
sub9 <- WhichCells(cluster8, idents = c('9'))

Idents(cluster8) <- cluster8$predicted.id
pred9 <- WhichCells(cluster8, idents = '9')

int9 <- intersect(sub9, pred9) #906

#cluster 10 
Idents(cluster8) <- cluster8$seurat_clusters
sub10 <- WhichCells(cluster8, idents = c('10'))

Idents(cluster8) <- cluster8$predicted.id
pred10 <- WhichCells(cluster8, idents = '10')

int10 <- intersect(sub10, pred10) #875

#cluster 11
Idents(cluster8) <- cluster8$seurat_clusters
sub11 <- WhichCells(cluster8, idents = c('11'))

Idents(cluster8) <- cluster8$predicted.id
pred11 <- WhichCells(cluster8, idents = '11')

int11 <- intersect(sub11, pred11) #768

#cluster 12
Idents(cluster8) <- cluster8$seurat_clusters
sub12 <- WhichCells(cluster8, idents = c('13'))

Idents(cluster8) <- cluster8$predicted.id
pred12 <- WhichCells(cluster8, idents = '12')

int12 <- intersect(sub12, pred12) #415

#cluster 13
Idents(cluster8) <- cluster8$seurat_clusters
sub13 <- WhichCells(cluster8, idents = c('14'))

Idents(cluster8) <- cluster8$predicted.id
pred13 <- WhichCells(cluster8, idents = '13')

int13 <- intersect(sub13, pred13) #315

#cluster 14
Idents(cluster8) <- cluster8$seurat_clusters
sub14 <- WhichCells(cluster8, idents = c('15'))

Idents(cluster8) <- cluster8$predicted.id
pred14 <- WhichCells(cluster8, idents = '14')

int14 <- intersect(sub14, pred14) #221

#cluster 15
Idents(cluster8) <- cluster8$seurat_clusters
sub15 <- WhichCells(cluster8, idents = c('16'))

Idents(cluster8) <- cluster8$predicted.id
pred15 <- WhichCells(cluster8, idents = '15')

int15 <- intersect(sub15, pred15) #60

#total 
#pre:55641
#52960
total8 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total8,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query8_low_res_barcodes_10_28_24.csv', row.names=F)

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

cluster9 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub_res20_5_low_res_cluster_pred.rds")

cluster9$predicted.id <- factor(cluster9$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster9, reduction = 'umap', label = TRUE) + DimPlot(cluster9, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

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

int0 <- intersect(sub0, pred0) #13842

#cluster 1
Idents(cluster9) <- cluster9$seurat_clusters
sub1 <- WhichCells(cluster9, idents = c('1'))

Idents(cluster9) <- cluster9$predicted.id
pred1 <- WhichCells(cluster9, idents = '1')

int1 <- intersect(sub1, pred1) #10524

#cluster 2
Idents(cluster9) <- cluster9$seurat_clusters
sub2 <- WhichCells(cluster9, idents = c('2'))

Idents(cluster9) <- cluster9$predicted.id
pred2 <- WhichCells(cluster9, idents = '2')

int2 <- intersect(sub2, pred2) #5528

#cluster 3
Idents(cluster9) <- cluster9$seurat_clusters
sub3 <- WhichCells(cluster9, idents = c('3'))

Idents(cluster9) <- cluster9$predicted.id
pred3 <- WhichCells(cluster9, idents = '3')

int3 <- intersect(sub3, pred3) #5658

#cluster 4 
Idents(cluster9) <- cluster9$seurat_clusters
sub4 <- WhichCells(cluster9, idents = c('4'))

Idents(cluster9) <- cluster9$predicted.id
pred4 <- WhichCells(cluster9, idents = '4')

int4 <- intersect(sub4, pred4) #4291

#cluster 5
Idents(cluster9) <- cluster9$seurat_clusters
sub5 <- WhichCells(cluster9, idents = c('5'))

Idents(cluster9) <- cluster9$predicted.id
pred5 <- WhichCells(cluster9, idents = '5')

int5 <- intersect(sub5, pred5) #3417

#cluster 6
Idents(cluster9) <- cluster9$seurat_clusters
sub6 <- WhichCells(cluster9, idents = c('6'))

Idents(cluster9) <- cluster9$predicted.id
pred6 <- WhichCells(cluster9, idents = '6')

int6 <- intersect(sub6, pred6) #2456

#cluster 7 
Idents(cluster9) <- cluster9$seurat_clusters
sub7 <- WhichCells(cluster9, idents = c('4', '8'))

Idents(cluster9) <- cluster9$predicted.id
pred7 <- WhichCells(cluster9, idents = '7')

int7 <- intersect(sub7, pred7) #2255

#cluster 8
Idents(cluster9) <- cluster9$seurat_clusters
sub8 <- WhichCells(cluster9, idents = c('7'))

Idents(cluster9) <- cluster9$predicted.id
pred8 <- WhichCells(cluster9, idents = '8')

int8 <- intersect(sub8, pred8) #1933

#cluster 9 
Idents(cluster9) <- cluster9$seurat_clusters
sub9 <- WhichCells(cluster9, idents = c('9'))

Idents(cluster9) <- cluster9$predicted.id
pred9 <- WhichCells(cluster9, idents = '9')

int9 <- intersect(sub9, pred9) #943

#cluster 10 
Idents(cluster9) <- cluster9$seurat_clusters
sub10 <- WhichCells(cluster9, idents = c('10'))

Idents(cluster9) <- cluster9$predicted.id
pred10 <- WhichCells(cluster9, idents = '10')

int10 <- intersect(sub10, pred10) #883

#cluster 11
Idents(cluster9) <- cluster9$seurat_clusters
sub11 <- WhichCells(cluster9, idents = c('11'))

Idents(cluster9) <- cluster9$predicted.id
pred11 <- WhichCells(cluster9, idents = '11')

int11 <- intersect(sub11, pred11) #757

#cluster 12
Idents(cluster9) <- cluster9$seurat_clusters
sub12 <- WhichCells(cluster9, idents = c('12'))

Idents(cluster9) <- cluster9$predicted.id
pred12 <- WhichCells(cluster9, idents = '12')

int12 <- intersect(sub12, pred12) #415

#cluster 13
Idents(cluster9) <- cluster9$seurat_clusters
sub13 <- WhichCells(cluster9, idents = c('13'))

Idents(cluster9) <- cluster9$predicted.id
pred13 <- WhichCells(cluster9, idents = '13')

int13 <- intersect(sub13, pred13) #351

#cluster 14
Idents(cluster9) <- cluster9$seurat_clusters
sub14 <- WhichCells(cluster9, idents = c('9'))

Idents(cluster9) <- cluster9$predicted.id
pred14 <- WhichCells(cluster9, idents = '14')

int14 <- intersect(sub14, pred14) #231

#cluster 15
Idents(cluster9) <- cluster9$seurat_clusters
sub15 <- WhichCells(cluster9, idents = c('4'))

Idents(cluster9) <- cluster9$predicted.id
pred15 <- WhichCells(cluster9, idents = '15')

int15 <- intersect(sub15, pred15) #64

#total 
#pre:55641
#53548
total9 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total9,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query9_low_res_barcodes_10_28_24.csv', row.names=F)

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

cluster10 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/colon_sub2_res20_5_low_res_cluster_pred.rds")

cluster10$predicted.id <- factor(cluster10$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

#visualize the umap 
DimPlot(cluster10, reduction = 'umap', label = TRUE) + DimPlot(cluster10, reduction = 'umap', group.by = 'predicted.id', label = TRUE)

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

int0 <- intersect(sub0, pred0) #13355

#cluster 1
Idents(cluster10) <- cluster10$seurat_clusters
sub1 <- WhichCells(cluster10, idents = c('1'))

Idents(cluster10) <- cluster10$predicted.id
pred1 <- WhichCells(cluster10, idents = '1')

int1 <- intersect(sub1, pred1) #10526

#cluster 2
Idents(cluster10) <- cluster10$seurat_clusters
sub2 <- WhichCells(cluster10, idents = c('3'))

Idents(cluster10) <- cluster10$predicted.id
pred2 <- WhichCells(cluster10, idents = '2')

int2 <- intersect(sub2, pred2) #5024

#cluster 3
Idents(cluster10) <- cluster10$seurat_clusters
sub3 <- WhichCells(cluster10, idents = c('2'))

Idents(cluster10) <- cluster10$predicted.id
pred3 <- WhichCells(cluster10, idents = '3')

int3 <- intersect(sub3, pred3) #5692

#cluster 4
Idents(cluster10) <- cluster10$seurat_clusters
sub4 <- WhichCells(cluster10, idents = c('4'))

Idents(cluster10) <- cluster10$predicted.id
pred4 <- WhichCells(cluster10, idents = '4')

int4 <- intersect(sub4, pred4) #4331

#cluster 5
Idents(cluster10) <- cluster10$seurat_clusters
sub5 <- WhichCells(cluster10, idents = c('5'))

Idents(cluster10) <- cluster10$predicted.id
pred5 <- WhichCells(cluster10, idents = '5')

int5 <- intersect(sub5, pred5) #3676

#cluster 6
Idents(cluster10) <- cluster10$seurat_clusters
sub6 <- WhichCells(cluster10, idents = c('6'))

Idents(cluster10) <- cluster10$predicted.id
pred6 <- WhichCells(cluster10, idents = '6')

int6 <- intersect(sub6, pred6) #2481

#cluster 7
Idents(cluster10) <- cluster10$seurat_clusters
sub7 <- WhichCells(cluster10, idents = c('4', '8'))

Idents(cluster10) <- cluster10$predicted.id
pred7 <- WhichCells(cluster10, idents = '7')

int7 <- intersect(sub7, pred7) #2404

#cluster 8 
Idents(cluster10) <- cluster10$seurat_clusters
sub8 <- WhichCells(cluster10, idents = c('7'))

Idents(cluster10) <- cluster10$predicted.id
pred8 <- WhichCells(cluster10, idents = '8')

int8 <- intersect(sub8, pred8) #1962

#cluster 9 
Idents(cluster10) <- cluster10$seurat_clusters
sub9 <- WhichCells(cluster10, idents = c('9'))

Idents(cluster10) <- cluster10$predicted.id
pred9 <- WhichCells(cluster10, idents = '9')

int9 <- intersect(sub9, pred9) #923

#cluster 10 
Idents(cluster10) <- cluster10$seurat_clusters
sub10 <- WhichCells(cluster10, idents = c('10'))

Idents(cluster10) <- cluster10$predicted.id
pred10 <- WhichCells(cluster10, idents = '10')

int10 <- intersect(sub10, pred10) #836

#cluster 11
Idents(cluster10) <- cluster10$seurat_clusters
sub11 <- WhichCells(cluster10, idents = c('11'))

Idents(cluster10) <- cluster10$predicted.id
pred11 <- WhichCells(cluster10, idents = '11')

int11 <- intersect(sub11, pred11) #787

#cluster 12
Idents(cluster10) <- cluster10$seurat_clusters
sub12 <- WhichCells(cluster10, idents = c('12'))

Idents(cluster10) <- cluster10$predicted.id
pred12 <- WhichCells(cluster10, idents = '12')

int12 <- intersect(sub12, pred12) #415

#cluster 13
Idents(cluster10) <- cluster10$seurat_clusters
sub13 <- WhichCells(cluster10, idents = c('13'))

Idents(cluster10) <- cluster10$predicted.id
pred13 <- WhichCells(cluster10, idents = '13')

int13 <- intersect(sub13, pred13) #327

#cluster 14
Idents(cluster10) <- cluster10$seurat_clusters
sub14 <- WhichCells(cluster10, idents = c('14'))

Idents(cluster10) <- cluster10$predicted.id
pred14 <- WhichCells(cluster10, idents = '14')

int14 <- intersect(sub14, pred14) #236

#cluster 15
Idents(cluster10) <- cluster10$seurat_clusters
sub15 <- WhichCells(cluster10, idents = c('8'))

Idents(cluster10) <- cluster10$predicted.id
pred15 <- WhichCells(cluster10, idents = '15')

int15 <- intersect(sub15, pred15) #63

#total 
#pre:55641
#53038
total10 <- paste(c(int0, int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12, int13, int14, int15))
write.csv(total10,file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query10_low_res_barcodes_10_28_24.csv', row.names=F)

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

#---------Identify stably assigned cells---------------#

#load in barcodes 
total1 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query1_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)
total2 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query2_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)
total3 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query3_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)
total4 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query4_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)
total5 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query5_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)
total6 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query6_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)
total7 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query7_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)
total8 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query8_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)
total9 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query9_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)
total10 <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/query10_low_res_barcodes_10_28_24.csv', what = "", sep = ",", skip = 1)

#merge the barcodes into a list 
#533183 barcodes 
total_merge <- c(total1, total2, total3, total4, total5, total6, total7, total8, total9, total10)

#counts the number of times each barcode is present
value_counts <- table(total_merge)
table(value_counts)

#define a threshold - more times a barcode is present, more stable it is
threshold <- 4

#filter values that meet your threshold 
filtered_values <- names(value_counts[value_counts >= threshold])  

#1: 111018 - keeps ~99.76% of data
#2: 110091 - keeps ~98.9% of data 
#3: 108318 - keeps ~97.34% of data
#4: 105066 - keeps ~94.41% of data ----select threshold of 4 (biggest jump after threshold of 4)
#5: 98690 - keeps ~88.68% of data

#save the filtered barcodes 
write.csv(filtered_values, file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/filtered_low_res_barcodes_10_29_24.csv', row.names = F)

#subset the colon seurat object on the server 

#load SO
colon <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_ileum_low_res_ref_10_21_24.rds")

#load in stably assigned cells 
filtered_values <- scan('/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/filtered_low_res_barcodes_10_29_24.csv', what = "", sep = ",", skip = 1)

#subset colon based on stable barcodes 
colon_sub <- subset(colon, cells = filtered_values)

#save 
saveRDS(colon_sub, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_low_res_stable_cell_10_29_24.rds")

#----------highlight the unstably assigned cells---------#

#load SO
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_rpca_10_17_24.rds")

#colon barcodes 
barcodes <- colnames(colon)

unstable <- setdiff(barcodes, filtered_values)

DimPlot(colon, reduction = "umap", cells.highlight = unstable, sizes.highlight = 0.1, raster = F) + scale_color_manual(labels = c("Stably Assigned Cells", "Unstably Assigned Cells"), values = c("grey", "red"))


#-----------Low resolution ARI results------------#
#10/30/2024

#merged ARI scores
df <- read.csv('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/df_ri_helm_batch1_13_colon_low_res_dune_merge_10_30_24.csv', row.names = NULL, sep = ",")

df <- df[,-1]
names(df) <- gsub("^paste0\\.", "", names(df))

columns <- colnames(df)

#plot the ARI - without paste0
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/figures/RI_helm_batch1_13_colon_low_res_all_param_merge_10_30_24.pdf", width = 15, height = 15)
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

#top 10 median ARI scores 
median_ari <- as.data.frame(ari_values$row_median)
rownames(median_ari) <- rownames(ari_values)
top_n(median_ari, 10)


#-----------Assess top 3 clustering results-------------#

#------cluster 1---------#

#parameters: res: 0.25, PC: 30, HVG: 2000

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_30_2000_rpca_11_2_24.rds")

colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
colon$batch <- factor(colon$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(colon, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(colon, features = 'nCount_RNA', pt.size = 0)
VlnPlot(colon, features = 'percent.mt', pt.size = 0)

#visualize PCA results
DimPlot(colon, reduction = 'pca', group.by = "seurat_clusters", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "batch", raster = F)

#visualize UMAP
DimPlot(colon, reduction = 'umap', group.by = "seurat_clusters", label = TRUE, raster = F)
DimPlot(colon, reduction = 'umap', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'umap', group.by = "batch", raster = F)

#proportion of cells per cluster 
colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
Idents(colon) <- colon$seurat_clusters
pt1 <- table(Idents(colon), colon$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18'))

library(randomcoloR)
no_of_colors <- 19

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
Idents(colon) <- colon$batch
pt2 <- table(Idents(colon), colon$seurat_clusters)
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

#------cluster 2--------#

#parameters: res: 0.50, PC: 30, HVG: 2000

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_50_30_2000_rpca_11_2_24.rds")

colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
colon$batch <- factor(colon$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(colon, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(colon, features = 'nCount_RNA', pt.size = 0)
VlnPlot(colon, features = 'percent.mt', pt.size = 0)

#visualize PCA results
DimPlot(colon, reduction = 'pca', group.by = "seurat_clusters", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "batch", raster = F)

#visualize UMAP
DimPlot(colon, reduction = 'umap', group.by = "seurat_clusters", label = TRUE, raster = F)
DimPlot(colon, reduction = 'umap', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'umap', group.by = "batch", raster = F)

#proportion of cells per cluster 
colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
Idents(colon) <- colon$seurat_clusters
pt1 <- table(Idents(colon), colon$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23'))

library(randomcoloR)
no_of_colors <- 24

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
Idents(colon) <- colon$batch
pt2 <- table(Idents(colon), colon$seurat_clusters)
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


#------cluster 3--------#

#parameters: res: 0.25, PC: 15, HVG: 5000

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_15_5000_rpca_11_2_24.rds")

colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
colon$batch <- factor(colon$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(colon, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(colon, features = 'nCount_RNA', pt.size = 0)
VlnPlot(colon, features = 'percent.mt', pt.size = 0)

#visualize PCA results
DimPlot(colon, reduction = 'pca', group.by = "seurat_clusters", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "batch", raster = F)

#visualize UMAP
DimPlot(colon, reduction = 'umap', group.by = "seurat_clusters", label = TRUE, raster = F)
DimPlot(colon, reduction = 'umap', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'umap', group.by = "batch", raster = F)

#proportion of cells per cluster 
colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
Idents(colon) <- colon$seurat_clusters
pt1 <- table(Idents(colon), colon$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

library(randomcoloR)
no_of_colors <- 16

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
Idents(colon) <- colon$batch
pt2 <- table(Idents(colon), colon$seurat_clusters)
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

#--------choose cluster 1--------#

#subset doublet cluster - cluster 18 

#do this on server 
colon <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_30_2000_rpca_11_2_24.rds")

#check Idents
Idents(colon) <- colon$seurat_clusters

#subset cluster 18 
colon_sub <- subset(colon, ident = "18", invert = T)

#save the clusters
saveRDS(colon_sub, file = "/storage/home/swashburn30/Helmsley/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_30_2000_sub_11_4_24.rds")

#--------Assess colon cluster 1 subset-----------#

#11/5/2024

#parameters: res: 0.25, PC: 30, HVG: 2000 - subset

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_30_2000_rpca_sub_cluster_11_4_24.rds")

colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
colon$batch <- factor(colon$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

#visualize distribution of data (clusters)
VlnPlot(colon, features = 'nFeature_RNA', pt.size = 0)
VlnPlot(colon, features = 'nCount_RNA', pt.size = 0)
VlnPlot(colon, features = 'percent.mt', pt.size = 0)

#visualize PCA results
DimPlot(colon, reduction = 'pca', group.by = "seurat_clusters", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'pca', group.by = "batch", raster = F)

#visualize UMAP
DimPlot(colon, reduction = 'umap', group.by = "seurat_clusters", label = TRUE, raster = F)
DimPlot(colon, reduction = 'umap', group.by = "orig.ident", raster = F)
DimPlot(colon, reduction = 'umap', group.by = "batch", raster = F)

#proportion of cells per cluster 
colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
Idents(colon) <- colon$seurat_clusters
pt1 <- table(Idents(colon), colon$orig.ident)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)

pt1$Var1 <- factor(pt1$Var1, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18'))

library(randomcoloR)
no_of_colors <- 19

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
Idents(colon) <- colon$batch
pt2 <- table(Idents(colon), colon$seurat_clusters)
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


#--------Preliminary annotation of colon cluster 1 subset-----------#


DefaultAssay(colon) <- "RNA"

Idents(colon) <- colon$seurat_clusters

#-----T cells-----#

#naive/central memory CD4 T cell
DotPlot(colon, features = c("CD3D", 'CD4', 'CCR7', 'SELL'))

#tissue resident memory CD4 T cell
DotPlot(colon, features = c("CD3D", 'CD4', 'ITGAE', 'ITGA1', 'SPRY1'))

#tissue resident memory T helper 17 
DotPlot(colon, features = c("CD3D", 'CD4', 'ITGAE', 'ITGA1', 'SPRY1', 'IL17A', 'CCR6', 'CCL20'))

#naive t follicular helper 
DotPlot(colon, features = c("CD3D", 'CD4', 'SELL', 'CXCR5', 'BCL6', 'CD40LG'))

#t follicular helper cell 
DotPlot(colon, features = c("CD3D", 'CD4', 'CXCR5', 'BCL6', 'PDCD1', 'CD40LG'))

#regulatory t cell 
DotPlot(colon, features = c("CD3D", 'CD4', 'FOXP3', 'CTLA4', 'IL2RA'))

#regulatory t cell IL10+
DotPlot(colon, features = c("CD3D", 'CD4', 'FOXP3', 'CTLA4', 'IL2RA', 'IL10'))

#naive/central memory CD8 t cell 
DotPlot(colon, features = c("CD3D", 'CD8A', 'CD8B', 'CCR7', 'SELL'))

#tissue resident memory CD8 t cell 
DotPlot(colon, features = c("CD3D", 'CD8A', 'CD8B', 'ITGAE', 'ITGA1', 'SPRY1'))

#tissue resident memory/effector memory Cd8 t cell 
DotPlot(colon, features = c("CD3D", 'CD8A', 'CD8B', 'GZMK', 'CRTAM', 'EOMES'))

#naive gamma delta t cell 
DotPlot(colon, features = c("CD3D", 'SELL', 'KLRC2', 'TRDC', 'TRGC1', 'KIR2DL4', 'CCR7'))

#gamma delta t cell 
DotPlot(colon, features = c("CD3D", 'KLRC2', 'TRDC', 'TRGC1', 'KIR2DL4'))

#MAIT (mucosal associated invariant t cell)
DotPlot(colon, features = c("CD3D", 'SLC4A10', 'TRAV1-2'))

#CD16+ NK cell 
DotPlot(colon, features = c('KLRF1', 'NKG7', 'FCGR3A', 'GNLY', 'GZMB'))

#CD56 bright natural killer cell 
DotPlot(colon, features = c('KLRF1', 'NKG7', 'XCL1', 'IL2RB', 'NCR1', 'FCER1G', 'NCAM1'))

#ILC3 innate lymphoid cell type 3 (cd3d negative)
DotPlot(colon, features = c('IL7R', 'RORC', 'KIT', 'LST1', 'PCDH9', 'IL1R1', 'IL23R'))

#cycling t or NK cell
DotPlot(colon, features = c('CD3D', 'MKI67', 'TOP2A'))

#--------myeloid------------#

#type 1 conventional DC 
DotPlot(colon, features = c('HLA-DRA', 'HLA-DPA1', 'CLEC9A', 'XCR1', 'BATF3', 'CADM1', 'RAB7B'))

#type 2 conventional DC
DotPlot(colon, features = c('HLA-DRA', 'HLA-DPA1', 'CLEC10A', 'FCER1A', 'CD1C'))

#migratory DC
DotPlot(colon, features = c('HLA-DRA', 'HLA-DPA1', 'CCR7', 'LAMP3'))

#plasmacytoid DC
DotPlot(colon, features = c('HLA-DRA', 'HLA-DPA1', 'IRF7', 'CLEC4C', 'JCHAIN', 'LILRA4', 'GZMB'))

#langerhans DC
DotPlot(colon, features = c('HLA-DRA', 'HLA-DPA1', 'ITGAX', 'IL22RA2', 'CD207', 'RUNX3'))

#monocyte
DotPlot(colon, features = c('FCN1', 'S100A8', 'S100A9', 'IL1B', 'EREG', 'NAMPT', 'PLAUR', 'VCAN', 'FPR1', 'CD3D00E'))

#LYVE1 macrophage
DotPlot(colon, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'LYVE1', 'RNASE1', 'FOLR2'))

#MMP9 macrophage
DotPlot(colon, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'MMP9', 'PLA2G2D', 'ADAMDEC1'))

#TREM2 macrophage 
DotPlot(colon, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'TREM2', 'ACP5', 'CTSD', 'CSTB'))

#CD5L macrophage
DotPlot(colon, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'CD5L', 'VCAM1', 'CXCL12', 'PDK4', 'RBP7'))

#macrophage
DotPlot(colon, features = c('CD163', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'CD209'))

#mast cell 
DotPlot(colon, features = c('CD69', 'KIT', 'TPSB2', 'TPSAB1'))

#eosinophil/basophil
DotPlot(colon, features = c('GATA2', 'CNRIP1', 'PRG2', 'GIHCG', 'CLC'))

#erythrocytes
DotPlot(colon, features = c('GATA1', 'HBZ', 'HBE1', 'HBG1'))

#monocyte/neutrophil progenitor
DotPlot(colon, features = c('FCN1', 'S100A8', 'S100A9', 'MPO', 'RETN', 'RNASE2', 'PCLAF'))

#megakaryocyte/platelet
DotPlot(colon, features = c('GATA1', 'TAL1', 'MMRN1', 'CMTM5', 'MPIG6B', 'ITGA2B', 'PF4'))

#----------b cells-------------#

#pre b cell 
DotPlot(colon, features = c('CD19', 'HLA-DRA', 'CD79B', 'SPIB', 'TCL1A', 'CD3D7'))

#pro b cell 
DotPlot(colon, features = c('CD19', 'HLA-DRA', 'CD79B', 'IGLL1', 'RAG1', 'DNTT', 'VPREB3'))

#naive b cell 
DotPlot(colon, features = c('CD19', 'HLA-DRA', 'CD79A', 'MS4A1', 'SELL', 'TCL1A', 'IGHD'))

#memory b cell 
DotPlot(colon, features = c('CD19', 'HLA-DRA', 'CD79A', 'MS4A1', 'CD27', 'TNFSF13B'))

#germinal center b cell 1
DotPlot(colon, features = c('CD19', 'HLA-DRA', 'CD79A', 'MS4A1', 'MKI67', 'AICDA', 'BCL6', 'SUGCT'))

#plasmablast
DotPlot(colon, features = c('MZB1', 'JCHAIN', 'XBP1'))

#IgM plasma cell 
DotPlot(colon, features = c('MZB1', 'JCHAIN', 'IGHM'))

#IgA plasma cell 
DotPlot(colon, features = c('MZB1', 'JCHAIN', 'IGHA1', 'IGHA2'))

#IgG plasma cell 
DotPlot(colon, features = c('MZB1', 'JCHAIN', 'IGHG3', 'IGHG1', 'IGHG2', 'IGHG4'))

#---------endothelial cells-----------#

#capillary endothelial cell 
DotPlot(colon, features = c('PECAM1', 'CD3D6', 'RGCC', 'COL4A1', 'COL4A2', 'IL32', 'MCAM', 'MYO1B'))

#arterial endothelial cell 
DotPlot(colon, features = c('PECAM1', 'CD3D6', 'GJA4', 'HEY1', 'CXCL12', 'SEMA3G', 'IGFBP3', 'FBLN2', 'FBLN5', 'ELN', 'BTNL9', 'ALPL'))

#venous endothelial cell
DotPlot(seur_obj, features = c('PECAM1', 'CD3D6', 'ACKR1', 'CCL14', 'SELE', 'TNFRSF6B'))

#lymphatic endothelial cell 
DotPlot(seur_obj, features = c('PECAM1', 'CD3D6', 'CCL21', 'TFF3', 'PROX2', 'NTS'))

#cycling endothelial cell 
DotPlot(seur_obj, features = c('PECAM1', 'CD3D6', 'MKI67', 'TOP2A'))

#-------------mesenchymal cells-------------#

#vascular smooth muscle cell 
DotPlot(colon, features = c('VIM', 'DCN', 'TAGLN', 'ACTA2', 'TPM2', 'MYH11', 'RERGL', 'MUSTN1', 'LBH', 'NET1', 'MAP3K20'))

#pericyte
DotPlot(colon, features = c('VIM', 'DCN', 'TAGLN', 'ACTA2', 'TPM2', 'MYH11', 'COX4I2', 'HIGD1B', 'RGS5', 'NDUFA4L2'))

#immune recruiting pericyte
DotPlot(seur_obj, features = c('VIM', 'DCN', 'TAGLN', 'ACTA2', 'TPM2', 'MYH11', 'GPC3', 'COL14A1', 'ECRG4', 'ID4', 'FHL2', 'CXCL12'))

#myofibroblast
DotPlot(colon, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'ACTG2', 'HHIP', 'SOSTDC1', 'NPNT'))

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
DotPlot(seur_obj, features = c('VIM', 'DCN', 'VCAN', 'PDGFRA', 'KCNN3', 'THBS4', 'FNDC1', 'PPFIBP1'))

#mesothelium 
DotPlot(seur_obj, features = c('VIM', 'DCN', 'UPK3B', 'MSLN', 'SLPI', 'PLAT', 'KRT19'))

#-----------epithelium-------------#

#---------large intestine----------#

#colonocyte
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'CA1', 'GPT'))

#BEST4 colonocyte
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'BEST4', 'CA7', 'OTOP2'))

#enteroendocrine cell 
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'CHGA', 'PCSK1N', 'SCT', 'SCGN', 'NEUROD1'))

#goblet cell 
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'MUC2', 'TFF3', 'FCGBP', 'ZG16'))

#mature colonocyte
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'SLC26A3', 'AQP8', 'CEACAM7'))

#deep crypt secretory cell 
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'MUC17', 'TFF1', 'CD55', 'TM4SF1', 'DUOX2', 'DUOXA2'))

#tuft cell 
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'SH2D6', 'LRMP', 'MATK', 'FYB1', 'HPGDS', 'POU2F3', 'TRPM5'))

#enterocyte
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'FABP2', 'APOA1', 'ALDOB'))

#transit amplifying cell 
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'MKI67', 'TOP2A', 'PCLAF', 'PCNA'))

#stem cell 
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'LGR5', 'RGMB', 'ASCL2', 'OLFM4'))

#enteroendocrine progenitor
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'SOX4', 'SCGN', 'NEUROD1'))

#proliferating deep crypt secretory cell
DotPlot(colon, features = c('CDH1', 'KRT19', 'EPCAM', 'MUC17', 'TFF1', 'MKI67', 'PCLAF', 'PCNA', 'TOP2A'))

#comparisons 

#cluster 1 vs. 8 
c1_8 <- FindMarkers(colon, ident.1 = "1", ident.2 = "8")

#cluster 0 vs. 2
c0_2 <- FindMarkers(colon, ident.1 = "0", ident.2 = "2")

#cluster 5 vs. cluster 12
c5_12 <- FindMarkers(colon, ident.1 = "5", ident.2 = "12")

#-------preliminary colon annotation-----------#

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_30_2000_rpca_sub_cluster_11_4_24.rds")

colon$orig.ident <- factor(colon$orig.ident, levels = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"))
colon$batch <- factor(colon$batch, levels = c("batch1", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "batch13"))

colon@meta.data$cell_typev1 <- colon@meta.data$seurat_clusters
colon$cell_typev1 <- plyr::mapvalues(
  x = colon$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18'),
  to = c('Colonocyte', 'Naive B Cell', 'Mature Colonocyte', 'IgA Plasma Cell', 'Stem/TA Cell', 'CD4 T Cell', 'Macrophage', 'NK/CD8 T Cell', 'Memory B Cell', 'Goblet Cell', 'Mesenchymal Cell', 'Monocyte', 'Regulatory T Cell', 'IgG Plasma Cell', 'Cycling B Cell', 'Tuft Cell', 'Mast Cell', 'Endothelial Cell', 'BEST4 Colonocyte')
)

colon@meta.data$cell_typev2 <- colon@meta.data$seurat_clusters
colon$cell_typev2 <- plyr::mapvalues(
  x = colon$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18'),
  to = c('Colonocyte', 'Naive B Cell', 'Colonocyte', 'IgA Plasma Cell', 'Stem/TA Cell', 'CD4 T Cell', 'Macrophage', 'NK/CD8 T Cell', 'Memory B Cell', 'Goblet Cell', 'Mesenchymal Cell', 'Monocyte', 'Regulatory T Cell', 'IgG Plasma Cell', 'Cycling B Cell', 'Tuft Cell', 'Mast Cell', 'Endothelial Cell', 'Colonocyte')
)

colon$cell_typev1 <- factor(colon$cell_typev1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "IgA Plasma Cell", "IgG Plasma Cell", "CD4 T Cell", "Regulatory T Cell", "NK/CD8 T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem/TA Cell", "Colonocyte", "Mature Colonocyte", "BEST4 Colonocyte", "Goblet Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

colon$cell_typev2 <- factor(colon$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "IgA Plasma Cell", "IgG Plasma Cell", "CD4 T Cell", "Regulatory T Cell", "NK/CD8 T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem/TA Cell", "Colonocyte", "Goblet Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

#t cells 
blue <- c("#34bfe5", "#4a93ce", "#b0daf2") 

#green = b cells
green <- c("#05fc95", "#9de587", "#1da54f", "#61e291", "#2bd147")

#orange = myeloid 
orange <- c("#fcab85", "#d3b050", "#f75f13", "#cc5414")

#pink = epithelial
pink <- c("#f41f8d", "#ffb7de", "#ed6add", "#fd8cff", "#cc0076", "#ffbff9", "#f271d6", "#f279b5")

#purple = mesenchymal and endothelial 
purple <- c("#cd7af9", "#996cc9")

#UMAPs
DimPlot(colon, reduction = "umap", group.by = "cell_typev1", cols = c(green[1], green[3], green[2], green[4], green[5], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], pink[1], pink[2], pink[3], pink[4], pink[5], pink[8], purple[1], purple[2]), raster = F)

DimPlot(colon, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[3], green[2], green[4], green[5], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], pink[1], pink[2], pink[5], pink[8], purple[1], purple[2]), raster = F)


#stacked bar plots 

Idents(colon) <- colon$cell_typev1
pt2 <- table(Idents(colon), colon$orig.ident)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

pt2$Var1 <- factor(pt2$Var1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "IgA Plasma Cell", "IgG Plasma Cell", "CD4 T Cell", "Regulatory T Cell", "NK/CD8 T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem/TA Cell", "Colonocyte", "Mature Colonocyte", "BEST4 Colonocyte", "Goblet Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[3], green[2], green[4], green[5], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], pink[1], pink[2], pink[3], pink[4], pink[5], pink[8], purple[1], purple[2])) +
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

Idents(colon) <- colon$cell_typev2
pt2 <- table(Idents(colon), colon$orig.ident)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

pt2$Var1 <- factor(pt2$Var1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "IgA Plasma Cell", "IgG Plasma Cell", "CD4 T Cell", "Regulatory T Cell", "NK/CD8 T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem/TA Cell", "Colonocyte", "Goblet Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

plot(ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[3], green[2], green[4], green[5], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], pink[1], pink[2], pink[5], pink[8], purple[1], purple[2])) +
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


#broad cell annotation 
colon@meta.data$broad_cell <- colon@meta.data$cell_typev2
colon$broad_cell <- plyr::mapvalues(
  x = colon$cell_typev2,
  from = c('Colonocyte', 'Naive B Cell', 'IgA Plasma Cell', 'Stem/TA Cell', 'CD4 T Cell', 'Macrophage', 'NK/CD8 T Cell', 'Memory B Cell', 'Goblet Cell', 'Mesenchymal Cell', 'Monocyte', 'Regulatory T Cell', 'IgG Plasma Cell', 'Cycling B Cell', 'Tuft Cell', 'Mast Cell', 'Endothelial Cell'),
  to = c('Epithelial Cell', 'Immune Cell', 'Immune Cell', 'Epithelial Cell', 'Immune Cell', 'Immune Cell', 'Immune Cell', 'Immune Cell', 'Epithelial Cell', 'Mesenchymal Cell', 'Immune Cell', 'Immune Cell', 'Immune Cell', 'Immune Cell', 'Epithelial Cell', 'Immune Cell', 'Endothelial Cell')
)

#donor information
colon@meta.data$donor <- colon@meta.data$orig.ident
colon$donor <- plyr::mapvalues(
  x = colon$orig.ident,
  from = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"),
  to = c('donor1', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor15_s1', 'donor15_s2', 'donor16_s1', 'donor16_s2', 'donor17', 'donor18', 'donor19', 'donor20', 'donor21', 'donor22', 'donor23', 'donor24', 'donor25', 'donor26', 'donor27', 'donor28_s1', 'donor28_s2', 'donor29', 'donor30', 'donor31', 'donor32', 'donor33', 'donor34')
)


#cell proportion analysis 
library(speckle)
library(limma)
library(ggplot2)

props <- getTransformedProps(colon$cell_typev2, colon$orig.ident, transform="logit")

colon_metadata <- data.frame(colon[["orig.ident"]], colon[["batch"]], colon[["macro_IF"]], colon[["micro_IF"]], colon[["ancestry"]], colon[["sex"]], colon[["donor"]], stringsAsFactors = TRUE)

colon_metadata <- unique(colon_metadata)

colon_metadata <- colon_metadata %>% arrange(orig.ident)

#create design matrix to take into account group and individual information
macro_IF <- colon_metadata$macro_IF
donors <- colon_metadata$donor
micro_IF <- colon_metadata$micro_IF
sex <- colon_metadata$sex
ancestry <- colon_metadata$ancestry
batch <- colon_metadata$batch


#proportion difference across macro IF status
mod <- model.matrix(~macro_IF)

dupcor <- duplicateCorrelation(props$TransformedProps, design = mod, block = donors)

#consensus correlation 0.106864

#fit linear model accounting donor as a random effect 
fit1 <- lmFit(props$TransformedProps, design = mod, block = donors, correlation = dupcor$consensus.correlation)
fit1 <- eBayes(fit1)

summary(decideTests(fit1))
topTable(fit1)

#monocytes significantly associated with MacroIF status; adj p value = 0.0415

#proportion difference across micro IF status
mod2 <- model.matrix(~micro_IF)

dupcor2 <- duplicateCorrelation(props$TransformedProps, design = mod2, block = donors)

#consensus correlation -0.07474989

#fit linear model accounting donor as a random effect 
fit2 <- lmFit(props$TransformedProps, design = mod2, block = donors, correlation = dupcor2$consensus.correlation)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))
topTable(fit2)

#monocytes are significantly associated with micro IF status adj p value = 5.283797e-05  

#proportion difference across batch 
mod3 <- model.matrix(~batch)

dupcor3 <- duplicateCorrelation(props$TransformedProps, design = mod3, block = donors)

#consensus correlation -0.07474989

#fit linear model accounting donor as a random effect 
fit3 <- lmFit(props$TransformedProps, design = mod3, block = donors, correlation = dupcor3$consensus.correlation)
fit3 <- eBayes(fit3)

summary(decideTests(fit3))
topTable(fit3)

#adj p value is not significant for all batches

#proportion difference across ancestry 
mod4 <- model.matrix(~ancestry)

dupcor4 <- duplicateCorrelation(props$TransformedProps, design = mod4, block = donors)

#consensus correlation -0.07474989

#fit linear model accounting donor as a random effect 
fit4 <- lmFit(props$TransformedProps, design = mod4, block = donors, correlation = dupcor4$consensus.correlation)
fit4 <- eBayes(fit4)

summary(decideTests(fit4))
topTable(fit4)

#ancestry is not significant

#proportion difference across sex 
mod5 <- model.matrix(~sex)

dupcor5 <- duplicateCorrelation(props$TransformedProps, design = mod5, block = donors)

#consensus correlation -0.04201999

#fit linear model accounting donor as a random effect 
fit5 <- lmFit(props$TransformedProps, design = mod5, block = donors, correlation = dupcor5$consensus.correlation)
fit5 <- eBayes(fit5)

summary(decideTests(fit5))
topTable(fit5)

#sex is not significant 

#look at broad cell type annotations 

props <- getTransformedProps(colon$broad_cell, colon$orig.ident, transform="logit")

#proportion difference across macro IF status
mod <- model.matrix(~macro_IF)

dupcor <- duplicateCorrelation(props$TransformedProps, design = mod, block = donors)

#consensus correlation -0.2790446

#fit linear model accounting donor as a random effect 
fit1 <- lmFit(props$TransformedProps, design = mod, block = donors, correlation = dupcor$consensus.correlation)
fit1 <- eBayes(fit1)

summary(decideTests(fit1))
topTable(fit1)

#broad cell proportions are not associated with IF status 

#proportion difference across micro IF status
mod2 <- model.matrix(~micro_IF)

dupcor2 <- duplicateCorrelation(props$TransformedProps, design = mod2, block = donors)

#consensus correlation -0.07474989

#fit linear model accounting donor as a random effect 
fit2 <- lmFit(props$TransformedProps, design = mod2, block = donors, correlation = dupcor2$consensus.correlation)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))
topTable(fit2)

#broad cell proportions are not associated with IF status 

#proportion difference across batch 
mod3 <- model.matrix(~batch)

dupcor3 <- duplicateCorrelation(props$TransformedProps, design = mod3, block = donors)

#consensus correlation -0.07474989

#fit linear model accounting donor as a random effect 
fit3 <- lmFit(props$TransformedProps, design = mod3, block = donors, correlation = dupcor3$consensus.correlation)
fit3 <- eBayes(fit3)

summary(decideTests(fit3))
topTable(fit3)

#endothelial cell is marginally associated wtih batch (adjusted p value = 0.05159023)

#proportion difference across ancestry 
mod4 <- model.matrix(~ancestry)

dupcor4 <- duplicateCorrelation(props$TransformedProps, design = mod4, block = donors)

#consensus correlation -0.07474989

#fit linear model accounting donor as a random effect 
fit4 <- lmFit(props$TransformedProps, design = mod4, block = donors, correlation = dupcor4$consensus.correlation)
fit4 <- eBayes(fit4)

summary(decideTests(fit4))
topTable(fit4)

#ancestry is not significant

#proportion difference across sex 
mod5 <- model.matrix(~sex)

dupcor5 <- duplicateCorrelation(props$TransformedProps, design = mod5, block = donors)

#consensus correlation -0.04201999

#fit linear model accounting donor as a random effect 
fit5 <- lmFit(props$TransformedProps, design = mod5, block = donors, correlation = dupcor5$consensus.correlation)
fit5 <- eBayes(fit5)

summary(decideTests(fit5))
topTable(fit5)

#sex is not significant

#cell type v1 annotation 

props <- getTransformedProps(colon$cell_typev1, colon$orig.ident, transform="logit")

colon_metadata <- data.frame(colon[["orig.ident"]], colon[["batch"]], colon[["macro_IF"]], colon[["micro_IF"]], colon[["ancestry"]], colon[["sex"]], colon[["donor"]], stringsAsFactors = TRUE)

colon_metadata <- unique(colon_metadata)

colon_metadata <- colon_metadata %>% arrange(orig.ident)

#create design matrix to take into account group and individual information
macro_IF <- colon_metadata$macro_IF
donors <- colon_metadata$donor
micro_IF <- colon_metadata$micro_IF
sex <- colon_metadata$sex
ancestry <- colon_metadata$ancestry
batch <- colon_metadata$batch


#proportion difference across macro IF status
mod <- model.matrix(~macro_IF)

dupcor <- duplicateCorrelation(props$TransformedProps, design = mod, block = donors)

#consensus correlation 0.106864

#fit linear model accounting donor as a random effect 
fit1 <- lmFit(props$TransformedProps, design = mod, block = donors, correlation = dupcor$consensus.correlation)
fit1 <- eBayes(fit1)

summary(decideTests(fit1))
topTable(fit1)

#monocytes significantly associated with MacroIF status; adj p value = 0.04914667

#proportion difference across micro IF status
mod2 <- model.matrix(~micro_IF)

dupcor2 <- duplicateCorrelation(props$TransformedProps, design = mod2, block = donors)

#consensus correlation -0.07474989

#fit linear model accounting donor as a random effect 
fit2 <- lmFit(props$TransformedProps, design = mod2, block = donors, correlation = dupcor2$consensus.correlation)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))
topTable(fit2)

#monocytes are significantly associated with micro IF status adj p value = 6.675208e-05  

#proportion difference across batch 
mod3 <- model.matrix(~batch)

dupcor3 <- duplicateCorrelation(props$TransformedProps, design = mod3, block = donors)

#consensus correlation -0.07474989

#fit linear model accounting donor as a random effect 
fit3 <- lmFit(props$TransformedProps, design = mod3, block = donors, correlation = dupcor3$consensus.correlation)
fit3 <- eBayes(fit3)

summary(decideTests(fit3))
topTable(fit3)

#adj p value is not significant for all batches

#proportion difference across ancestry 
mod4 <- model.matrix(~ancestry)

dupcor4 <- duplicateCorrelation(props$TransformedProps, design = mod4, block = donors)

#consensus correlation -0.07474989

#fit linear model accounting donor as a random effect 
fit4 <- lmFit(props$TransformedProps, design = mod4, block = donors, correlation = dupcor4$consensus.correlation)
fit4 <- eBayes(fit4)

summary(decideTests(fit4))
topTable(fit4)

#ancestry is not significant

#proportion difference across sex 
mod5 <- model.matrix(~sex)

dupcor5 <- duplicateCorrelation(props$TransformedProps, design = mod5, block = donors)

#consensus correlation -0.04201999

#fit linear model accounting donor as a random effect 
fit5 <- lmFit(props$TransformedProps, design = mod5, block = donors, correlation = dupcor5$consensus.correlation)
fit5 <- eBayes(fit5)

summary(decideTests(fit5))
topTable(fit5)

#sex is not significant 


#----------Hierarchical Clustering-----------#


#based on cell proportions 

colon$cell_typev1 <- factor(colon$cell_typev1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "IgA Plasma Cell", "IgG Plasma Cell", "CD4 T Cell", "Regulatory T Cell", "NK/CD8 T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem/TA Cell", "Colonocyte", "Mature Colonocyte", "BEST4 Colonocyte", "Goblet Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

props <- getTransformedProps(colon$cell_typev1, colon$donor, transform="logit")

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
                 cutree_cols = 4,
                 color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
hmap


save_pheatmap_pdf(hmap, "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/colon/heatmap_batch1_13_colon_hierarchical_cluster_cellprop_v1_11_25_24.pdf", width = 8, height = 5)

#----------scITD group comparison------------#

#11/27/2024
#subset donor 4, 5, 8, 9 - was not kept in scITD analysis 
Idents(colon) <- colon$donor
colon_sub <- subset(colon, idents = c("donor4", "donor5", "donor8", "donor9"), invert = T)

#group 1 = positive score
#group 2 = negative score

colon_sub@meta.data$group <- colon_sub@meta.data$orig.ident
colon_sub$group <- plyr::mapvalues(
  x = colon_sub$orig.ident,
  from = c("helm_sam2", "helm_sam8", "helm_sam33", "helm_sam36", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"),
  to = c('group1', 'group1', "group1", 'group1', 'group1', 'group2', 'group1', 'group1', 'group2', 'group2', 'group1', 'group2', 'group1', 'group1', 'group2', 'group1', 'group1', 'group1', 'group2', 'group1', 'group1', 'group2', 'group1', 'group1', 'group2', 'group2', 'group2', 'group2', 'group2', 'group1', 'group2', 'group1')
)

#save the colon seurat object with group information 
saveRDS(colon_sub, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_groupInfo_11_27_24.rds")

#prep for dreamlet - have to convert SO to single cell experiment 

#code below 
DefaultAssay(colon) <- "RNA"
colon.sce <- as.SingleCellExperiment(colon)
saveRDS(colon.sce, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_group_sce_11_27_24.rds")

#create sce object of all of colon donors 
colon <- readRDS(file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_cell_annotation_11_26_24.rds")
DefaultAssay(colon) <- "RNA"
colon.sce <- as.SingleCellExperiment(colon)
saveRDS(colon.sce, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_sce_12_4_24.rds")

#---------group 1 colon donors for cell chat----------#

#load in colon group seurat object
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_groupInfo_11_27_24.rds")

#set idents to group 
Idents(colon) <- colon$group

#subset group 1
colon_g1 <- subset(colon, idents = "group1")

#save group 1
saveRDS(colon_g1, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_group1_12_5_24.rds")

#subset group 2
colon_g2 <- subset(colon, idents = "group2")

#save group 2
saveRDS(colon_g2, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_group2_12_5_24.rds")


#-----------myeloid cell activation ------------#
metadata <- pb@colData
metadata <- as.data.frame(metadata)

#load in myeloid genes 
myl <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/colon_pseuodbulk_myeloid_cell.csv")

rownames(myl) <- myl$X

#relaxed threshold 
mac_gene <- c("SLC11A1", "TLR2", "LRRK2", "C5AR1", "ITGAM", "CCL3", "FCGR3A", "CD93", "PTPRC", "TYROBP", "JAK2")

myl$gene <- rownames(myl)
myl <- myl[,-1]

myl_filt <- myl[myl$gene %in% mac_gene,]

myl_filt <- myl_filt[,-33] #11 genes 

myl_filt <- t(myl_filt)


#perform PCA 
myl_pca <- prcomp(myl_filt, scale = T)



all(row.names(myl_pca) == row.names(metadata)) #check that samples match ##TRUE

#variance explained by each PC
pc_eigenvalues <- myl_pca$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues #PC1: 61.7% of variance; PC2: 11.1% of variance 

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- myl_pca$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

pc_scores

#add sample name column 
metadata$sample <- rownames(metadata)
pc_score_metadata <- full_join(pc_scores, metadata, by = ("sample"))

pc_score_metadata$group <- factor(pc_score_metadata$group, levels = c("group1", "group2"))

plot <- pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(group))) +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Group")) +
  theme(text = element_text(size = 18)) + 
  xlab("PC1 (61.7%)") +
  ylab("PC2 (11.1%)") +
  scale_color_manual(labels = c("Group 1", "Group 2"), values = c('#F8766D', '#619CFF'))

png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/figures/PCA_colon_group_myl_mac_activation_DEG_relax_12_10_24.png", width = 1800, height = 1800, res = 300)
plot(plot)
dev.off()


#create PCA plot 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(group))) +
  geom_point(size = 3) 
#guides(color = guide_legend(title = "Group")) +
#theme(text = element_text(size = 10)) + 
#scale_color_manual(labels = c("Group 1", "Group 2"), values = c('#F8766D', '#619CFF'))

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(batch))) +
  geom_point() 

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(sex))) +
  geom_point() 

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(ancestry))) +
  geom_point() 


plot <- pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(ancestry))) +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Donor Ancestry")) +
  theme(text = element_text(size = 18)) + 
  xlab("PC1 (35.6%)") +
  ylab("PC2 (20.2%)") +
  scale_color_manual(labels = c("AA", "Asian", "EA"), values = c('#52B2BF', '#9DC183', "#D21F3C"))

png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_9/figures/batch1_9_update_3_20_24/GCA/pca_colon_g1_g2_myl_mac_inflam_ancestry_8_29_24.png", width = 2000, height = 1800, res = 300)
plot(plot)
dev.off()


pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(Macro_IF))) +
  geom_point() 

pc_score_metadata$macro_IF <- factor(pc_score_metadata$macro_IF, levels = c("NI", "IF"))
plot <- pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(macro_IF))) +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Macro IF Status")) +
  theme(text = element_text(size = 18)) + 
  xlab("PC1 (61.7%)") +
  ylab("PC2 (11.1%)") +
  scale_color_manual(labels = c("NI", "IF"), values = c('#0E4D92', '#A94064'))

png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/figures/PCA_colon_group_myl_mac_activation_macroIFstatus_DEG_relax_12_10_24.png", width = 2000, height = 1800, res = 300)
plot(plot)
dev.off()

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(micro_IF))) +
  geom_point() 

pc_score_metadata$micro_IF <- factor(pc_score_metadata$micro_IF, levels = c("NI", "IF"))
plot <- pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(micro_IF))) +
  geom_point(size = 3) +
  guides(color = guide_legend(title = "Micro IF Status")) +
  theme(text = element_text(size = 18)) + 
  xlab("PC1 (61.7%)") +
  ylab("PC2 (11.1%)") +
  scale_color_manual(labels = c("NI", "IF"), values = c('#0E4D92', '#A94064'))

png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/figures/PCA_colon_group_myl_mac_activation_DEG_relax_microIFstatus_12_10_24.png", width = 2000, height = 1800, res = 300)
plot(plot)
dev.off()



t.test(pc_score_metadata$PC1~pc_score_metadata$group)#p value = 8.801e-07
t.test(pc_score_metadata$PC2~pc_score_metadata$group) #p value = 0.7731
t.test(pc_score_metadata$PC1~pc_score_metadata$macro_IF) #p value = p-value = 0.003684
t.test(pc_score_metadata$PC2~pc_score_metadata$macro_IF) #p-value = 0.05735
t.test(pc_score_metadata$PC1~pc_score_metadata$micro_IF) #p-value = 0.0001235
t.test(pc_score_metadata$PC2~pc_score_metadata$micro_IF) #p value = 0.5319
batch_aov <- aov(pc_score_metadata$PC1~pc_score_metadata$batch)
summary(batch_aov) #p value =0.612

batch_aov <- aov(pc_score_metadata$PC2~pc_score_metadata$batch)
summary(batch_aov) #p value =0.28


an_aov <- aov(pc_score_metadata$PC2~pc_score_metadata$ancestry)
summary(an_aov) #p value = 0.0189

an_aov <- aov(pc_score_metadata$PC1~pc_score_metadata$ancestry)
summary(an_aov) #p value = 0.501

#PC1 is associated with group, macro, and micro IF status 

#-----------Macrophage activation module score and IF module score-----------_#
#perform on server 

library(UCell)

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_groupInfo_11_27_24.rds")

DefaultAssay(colon) <- "RNA"

colon$group <- factor(colon$group, levels = c("group1", "group2"))

#intersection between all mac genes in genes in colon data 

mac_genes <- read.csv(file = "/storage/home/swashburn30/Helmsley/Batch1_13/DEG_analysis/macrophage_act_genes_all.csv")

mac_genes <- mac_genes[,-1]

colon_gene <- rownames(colon)

mac_colon <- intersect(mac_genes, colon_gene)

colon <- AddModuleScore_UCell(colon, features = list(mac_colon), name = "mac_sig")

#save the colon SO with module score 
saveRDS(colon, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_group_mac_act_score_12_7_24.rds")


#----------IF module score---------#

if_gene <- read.csv(file = "/storage/home/swashburn30/Helmsley/Batch1_13/DEG_analysis/regulation_IF_gobp_genes.csv")

if_colon <- intersect(if_gene$GENE_SYMBOLS, colon_gene)

colon <- AddModuleScore_UCell(colon, features = list(if_colon), name = "IF_sig")

#both mac activation and IF regulation 
saveRDS(colon, file = "/storage/home/swashburn30/Helmsley/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_group_UCell_score_12_7_24.rds")


#---------macrophage activation and response to inflammation status relationship----------# 

colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_group_UCell_score_12_7_24.rds")

colon$donor <- factor(colon$donor, levels = c("donor1", "donor3", "donor6", "donor7", "donor10", "donor12", "donor13", "donor15_s2", "donor16_s2", "donor17", "donor19", "donor20", "donor21", "donor23", "donor24", "donor26", "donor27", "donor32", "donor34", "donor11", "donor14", "donor15_s1", "donor16_s1", "donor18", "donor22", "donor25", "donor28_s1", "donor28_s2", "donor29", "donor30", "donor31", "donor33"))

#calculate the average score for each module per donor
colon_metadata <- colon@meta.data

# Calculate average module scores for each donor
average_scores <- colon_metadata %>%
  group_by(donor) %>%
  summarize(across(starts_with("signature"), ~ mean(.x, na.rm = TRUE)))

average_scores$group <- NA

average_scores$group <- c("group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group1", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2", "group2")

# Change the point size, and shape

#IF vs. Macrophage
plot <- ggplot(average_scores, aes(x=signature_1mac_sig, y=signature_1IF_sig, color=group)) +
  geom_point(size=3) +
  geom_smooth(method = "lm", se = TRUE) +
  guides(color = guide_legend(title = "Group")) +
  theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  xlab("Mean Macrophage Activation Score") +
  ylab("Mean Regulation of Inflammatory Response Score") +
  scale_color_manual(labels = c("Group 1", "Group 2"), values = c('#D65F5F', '#7986CB'))

plot
#get the slope for each group 
#test for interaction between group and macrophage activation score 
model <- lm(signature_1IF_sig ~ signature_1mac_sig * group, data = average_scores)
summary(model)

#slope of group 1 is 0.471105 (signature_1mac_sig) and slope of group 2 is 0.800566 (signature_1mac_sig + signature_1mac_sig:groupgroup2)


#save the plot as a pdf 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure5_12_16_24/ScatterPlot_helm_batch1_13_colon_myeloid_cell_mac_activation_vs_inflam_response_score_12_16_24.pdf", width = 7, height = 7)
plot(plot)
dev.off()


