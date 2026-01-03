#Helmsley GCA Batch 1-13 Supplementary Figures

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
library(speckle)
library(limma)
library(factoextra)
library(tibble)


#----------Supplementary Figure 1----------#


#12/17/2024

#####Supp. Figure 1b#######

#load in ileum SO 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_udpate_11_26_24.rds")

#colors to use for batch 
batch_colors <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000", "#BBCDF1", "#1CF815", "#9B97E1", "#51CFC6", "#2319A5", "#B34C90", "#CC6677")

#order by donor (batch)
ileum$donor <- factor(ileum$donor, levels = c('donor1', 'donor2', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor17', 'donor18', 'donor20', 'donor21', 'donor23', 'donor24', 'donor26', 'donor27', 'donor29', 'donor30', 'donor31', 'donor33', 'donor34'))

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/VlnPlot_helm_batch1_13_ileum_nfeature_rna_12_17_24.pdf", width = 10, height = 7)
plot(VlnPlot(ileum, features = "nFeature_RNA", pt.size = 0, group.by = "donor", cols = c(batch_colors[1], batch_colors[1], batch_colors[1], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[3], batch_colors[3], batch_colors[4],batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[5], batch_colors[6], batch_colors[6], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[8], batch_colors[9], batch_colors[10], batch_colors[10], batch_colors[11], batch_colors[11], batch_colors[12])) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank()) + labs(title = "Gene Count"))
dev.off()

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/VlnPlot_helm_batch1_13_ileum_percentmt_12_17_24.pdf", width = 10, height = 7)
plot(VlnPlot(ileum, features = "percent.mt", pt.size = 0, group.by = "donor", cols = c(batch_colors[1], batch_colors[1], batch_colors[1], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[3], batch_colors[3], batch_colors[4],batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[5], batch_colors[6], batch_colors[6], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[8], batch_colors[9], batch_colors[10], batch_colors[10], batch_colors[11], batch_colors[11], batch_colors[12])) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank()) + labs(title = "Mitochondria Percent"))
dev.off()

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/VlnPlot_helm_batch1_13_ileum_ncount_rna_12_17_24.pdf", width = 10, height = 7)
plot(VlnPlot(ileum, features = "nCount_RNA", pt.size = 0, group.by = "donor", cols = c(batch_colors[1], batch_colors[1], batch_colors[1], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[3], batch_colors[3], batch_colors[4],batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[5], batch_colors[6], batch_colors[6], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[8], batch_colors[9], batch_colors[10], batch_colors[10], batch_colors[11], batch_colors[11], batch_colors[12])) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank()) + labs(title = "Transcript Count"))
dev.off()


#load in colon SO
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_cell_annotation_11_26_24.rds")

#order the donor 
colon$donor <- factor(colon$donor, levels = c('donor1', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor15_s1', 'donor15_s2', 'donor16_s1', 'donor16_s2', 'donor17', 'donor18', 'donor19', 'donor20', 'donor21', 'donor22', 'donor23', 'donor24', 'donor25', 'donor26', 'donor27', 'donor28_s1', 'donor28_s2', 'donor29', 'donor30', 'donor31', 'donor32', 'donor33', 'donor34'))

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/VlnPlot_helm_batch1_13_colon_nfeature_rna_12_17_24.pdf", width = 10.5, height = 7)
plot(VlnPlot(colon, features = 'nFeature_RNA', group.by = "donor", cols = c(batch_colors[1], batch_colors[1], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[3], batch_colors[3], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[6], batch_colors[6], batch_colors[6], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[8], batch_colors[8], batch_colors[9], batch_colors[9], batch_colors[9], batch_colors[10], batch_colors[10], batch_colors[11], batch_colors[11], batch_colors[11], batch_colors[12]), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank()) + labs(title = "Gene Count"))
dev.off()


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/VlnPlot_helm_batch1_13_colon_percentmt_12_17_24.pdf", width = 10.5, height = 7)
plot(VlnPlot(colon, features = 'percent.mt', group.by = "donor", cols = c(batch_colors[1], batch_colors[1], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[3], batch_colors[3], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[6], batch_colors[6], batch_colors[6], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[8], batch_colors[8], batch_colors[9], batch_colors[9], batch_colors[9], batch_colors[10], batch_colors[10], batch_colors[11], batch_colors[11], batch_colors[11], batch_colors[12]), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank()) + labs(title = "Mitochondria Percent"))
dev.off()


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/VlnPlot_helm_batch1_13_colon_ncount_rna_12_17_24.pdf", width = 10.5, height = 7)
plot(VlnPlot(colon, features = 'nCount_RNA', group.by = "donor", cols = c(batch_colors[1], batch_colors[1], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[3], batch_colors[3], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[6], batch_colors[6], batch_colors[6], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[8], batch_colors[8], batch_colors[9], batch_colors[9], batch_colors[9], batch_colors[10], batch_colors[10], batch_colors[11], batch_colors[11], batch_colors[11], batch_colors[12]), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank()) + labs(title = "Transcript Count"))
dev.off()


#load in rectum SO 
rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_cell_annotation_12_8_24.rds")

#order the rectum data 
rectum$donor <- factor(rectum$donor, levels = c("donor1", "donor3", "donor4", "donor5", "donor6", "donor7", "donor8", "donor9", 'donor10', "donor11", "donor12", "donor13", "donor14", "donor15", "donor16", "donor17", "donor18", "donor19", "donor20", "donor21", "donor22", "donor24", "donor25", "donor26", "donor27", "donor28", "donor29", "donor30", "donor31", "donor32", "donor33", "donor34"))


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/VlnPlot_helm_batch1_13_rectum_nfeature_rna_12_17_24.pdf", width = 10, height = 7)
plot(VlnPlot(rectum, features = 'nFeature_RNA', group.by = "donor", cols = c(batch_colors[1], batch_colors[1], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[3], batch_colors[3], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[6], batch_colors[6], batch_colors[6], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[8], batch_colors[8], batch_colors[9], batch_colors[9], batch_colors[10], batch_colors[10], batch_colors[11], batch_colors[11], batch_colors[11], batch_colors[12]), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank()) + labs(title = "Gene Count"))
dev.off()

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/VlnPlot_helm_batch1_13_rectum_percentmt_12_17_24.pdf", width = 10, height = 7)
plot(VlnPlot(rectum, features = 'percent.mt', group.by = "donor", cols = c(batch_colors[1], batch_colors[1], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[3], batch_colors[3], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[6], batch_colors[6], batch_colors[6], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[8], batch_colors[8], batch_colors[9], batch_colors[9], batch_colors[10], batch_colors[10], batch_colors[11], batch_colors[11], batch_colors[11], batch_colors[12]), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank()) + labs(title = "Mitochondria Percent"))
dev.off()

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/VlnPlot_helm_batch1_13_rectum_ncount_rna_12_17_24.pdf", width = 10, height = 7)
plot(VlnPlot(rectum, features = 'nCount_RNA', group.by = "donor", cols = c(batch_colors[1], batch_colors[1], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[2], batch_colors[3], batch_colors[3], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[4], batch_colors[5], batch_colors[5], batch_colors[5], batch_colors[6], batch_colors[6], batch_colors[6], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[7], batch_colors[8], batch_colors[8], batch_colors[9], batch_colors[9], batch_colors[10], batch_colors[10], batch_colors[11], batch_colors[11], batch_colors[11], batch_colors[12]), pt.size = 0) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank()) + labs(title = "Transcript Count"))
dev.off()

#----------Supplementary Figure 2----------#

#####Supp. Figure 2a#######

ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_rPCA_10_17_24.rds")

filt_barcodes <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/filtered_low_res_barcodes_10_24_24.csv', what = "", sep = ",", skip = 1)

barcodes <- colnames(ileum)

unstable <- setdiff(barcodes, filt_barcodes)

unstable_cells <- subset(ileum, cells = unstable)

#label unstable cells 
unstable_cells$stable_assignment <- "unstable"

unstable_cells@meta.data$donor <- unstable_cells@meta.data$orig.ident
unstable_cells$donor <- plyr::mapvalues(
  x = unstable_cells$orig.ident,
  from = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'),
  to = c('donor1', 'donor2', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor17', 'donor18', 'donor20', 'donor21', 'donor23', 'donor24', 'donor26', 'donor27', 'donor29', 'donor30', 'donor31', 'donor33', 'donor34')
)


#read in ileum SO - final SO 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_udpate_11_26_24.rds")

#add stable lable 
ileum$stable_assignment <- "stable"

#merge unstable and ileum 
ileum_unstable_merge <- merge(ileum, y = unstable_cells)

#set ident to donor
Idents(ileum_unstable_merge) <- ileum_unstable_merge$donor

#n feature RNA

VlnPlot(ileum_unstable_merge, features = c("nFeature_RNA"))

#color by stability assessment 
stability_color <- c("grey", "red")

#order by donor (batch)
ileum_unstable_merge$donor <- factor(ileum_unstable_merge$donor, levels = c('donor1', 'donor2', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor17', 'donor18', 'donor20', 'donor21', 'donor23', 'donor24', 'donor26', 'donor27', 'donor29', 'donor30', 'donor31', 'donor33', 'donor34'))

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/VlnPlot_helm_batch1_13_ileum_nfeature_stable_unstable_9_23_25.pdf", width = 5, height = 12)
plot(VlnPlot(ileum_unstable_merge, features = "nFeature_RNA", pt.size = 0, group.by = "donor", split.by = "stable_assignment", cols = stability_color) + geom_boxplot(width=0.2, position = position_dodge(width = 0.9),outlier.shape = NA,) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y = element_blank()) + labs(title = "Gene count") + coord_flip())
dev.off()


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/VlnPlot_helm_batch1_13_ileum_percentmt_stable_unstable_9_23_25.pdf", width = 5, height = 12)
plot(VlnPlot(ileum_unstable_merge, features = "percent.mt", pt.size = 0, group.by = "donor", split.by = "stable_assignment", cols = stability_color) + geom_boxplot(width=0.2, position = position_dodge(width = 0.9),outlier.shape = NA,) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y = element_blank()) + labs(title = "Mitochondria percent") + coord_flip())
dev.off()

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/VlnPlot_helm_batch1_13_ileum_ncount_stable_unstable_9_23_25.pdf", width = 5, height = 12)
plot(VlnPlot(ileum_unstable_merge, features = "nCount_RNA", pt.size = 0, group.by = "donor", split.by = "stable_assignment", cols = stability_color) + geom_boxplot(width=0.2, position = position_dodge(width = 0.9),outlier.shape = NA,) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(title = "UMI count") + coord_flip())
dev.off()

#colon 
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_rpca_10_17_24.rds")

filt_barcodes <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/filtered_low_res_barcodes_10_29_24.csv', what = "", sep = ",", skip = 1)

barcodes <- colnames(colon)

unstable <- setdiff(barcodes, filt_barcodes)

unstable_cells <- subset(colon, cells = unstable)

#stable assignment 
unstable_cells$stable_assignment <- "unstable"

unstable_cells@meta.data$donor <- unstable_cells@meta.data$orig.ident
unstable_cells$donor <- plyr::mapvalues(
  x = unstable_cells$orig.ident,
  from = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"),
  to = c('donor1', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor15_s1', 'donor15_s2', 'donor16_s1', 'donor16_s2', 'donor17', 'donor18', 'donor19', 'donor20', 'donor21', 'donor22', 'donor23', 'donor24', 'donor25', 'donor26', 'donor27', 'donor28_s1', 'donor28_s2', 'donor29', 'donor30', 'donor31', 'donor32', 'donor33', 'donor34')
)


####Supp. Figure 2b####
#load in filtered colon SO
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_cell_annotation_11_26_24.rds")

#add stable lable 
colon$stable_assignment <- "stable"

#merge unstable and colon 
colon_unstable_merge <- merge(colon, y = unstable_cells)

#order the donor 
colon_unstable_merge$donor <- factor(colon_unstable_merge$donor, levels = c('donor1', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor15_s1', 'donor15_s2', 'donor16_s1', 'donor16_s2', 'donor17', 'donor18', 'donor19', 'donor20', 'donor21', 'donor22', 'donor23', 'donor24', 'donor25', 'donor26', 'donor27', 'donor28_s1', 'donor28_s2', 'donor29', 'donor30', 'donor31', 'donor32', 'donor33', 'donor34'))


#set ident to donor
Idents(colon_unstable_merge) <- colon_unstable_merge$donor


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/VlnPlot_helm_batch1_13_colon_nfeature_stable_unstable_9_23_25.pdf", width = 5, height = 12.5)
plot(VlnPlot(colon_unstable_merge, features = "nFeature_RNA", pt.size = 0, group.by = "donor", split.by = "stable_assignment", cols = stability_color) + geom_boxplot(width=0.2, position = position_dodge(width = 0.9),outlier.shape = NA,) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(title = "Gene count") + coord_flip())
dev.off()


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/VlnPlot_helm_batch1_13_colon_percentmt_stable_unstable_9_23_25.pdf", width = 5, height = 12.5)
plot(VlnPlot(colon_unstable_merge, features = 'percent.mt', pt.size = 0, group.by = "donor", split.by = "stable_assignment", cols = stability_color) + theme(legend.position = "none") + geom_boxplot(width=0.2, position = position_dodge(width = 0.9),outlier.shape = NA,) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(title = "Mitochondria Percent") + coord_flip())
dev.off()


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/VlnPlot_helm_batch1_13_colon_ncount_rna_stable_unstable_9_23_25.pdf", width = 5, height = 12.5)
plot(VlnPlot(colon_unstable_merge, features = 'nCount_RNA', group.by = "donor", split.by = "stable_assignment", cols = stability_color, pt.size = 0) + geom_boxplot(width=0.2, position = position_dodge(width = 0.9),outlier.shape = NA,) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(title = "UMI Count") + coord_flip())
dev.off()


####Supp. Figure 2C####

#rectum 

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_rpca_10_17_24.rds")

#load in stably assigned cells 
filtered_values <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/filtered_low_res_barcodes_11_4_24.csv', what = "", sep = ",", skip = 1)

barcodes <- colnames(rectum)

unstable <- setdiff(barcodes, filtered_values)

unstable_cells <- subset(rectum, cells = unstable)

unstable_cells@meta.data$donor <- unstable_cells@meta.data$orig.ident
unstable_cells$donor <- plyr::mapvalues(
  x = unstable_cells$orig.ident,
  from = c("helm_sam3", "helm_sam9", "helm_sam25", "helm_sam30", "helm_sam34", "helm_sam37", "helm_sam46", "helm_sam49", "helm_sam51", "helm_sam55", "helm_sam58", "helm_sam61", "helm_sam64", "helm_sam67", "helm_sam73", "helm_sam79", "helm_sam85", "helm_sam87", "helm_sam90", "helm_sam93", "helm_sam97", "helm_sam103", "helm_sam111", "helm_sam114", "helm_sam121", "helm_sam130", "helm_sam136", "helm_sam149", "helm_sam155", "helm_sam158", "helm_sam164", "helm_sam173"),
  to = c("donor1", "donor3", "donor4", "donor5", "donor6", "donor7", "donor8", "donor9", 'donor10', "donor11", "donor12", "donor13", "donor14", "donor15", "donor16", "donor17", "donor18", "donor19", "donor20", "donor21", "donor22", "donor24", "donor25", "donor26", "donor27", "donor28", "donor29", "donor30", "donor31", "donor32", "donor33", "donor34")
)

unstable_cells$stable_assignment <- "unstable"

#load in filtered rectum SO 
rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_cell_annotation_12_8_24.rds")

rectum$stable_assignment <- "stable"

#merge 
rectum_unstable_merge <- merge(rectum, y = unstable_cells)



#order the rectum data 
rectum_unstable_merge$donor <- factor(rectum_unstable_merge$donor, levels = c("donor1", "donor3", "donor4", "donor5", "donor6", "donor7", "donor8", "donor9", 'donor10', "donor11", "donor12", "donor13", "donor14", "donor15", "donor16", "donor17", "donor18", "donor19", "donor20", "donor21", "donor22", "donor24", "donor25", "donor26", "donor27", "donor28", "donor29", "donor30", "donor31", "donor32", "donor33", "donor34"))

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/VlnPlot_helm_batch1_13_rectum_nfeature_rna_stable_ustable_9_23_25.pdf", width = 5, height = 12)
plot(VlnPlot(rectum_unstable_merge, features = 'nFeature_RNA', group.by = "donor", split.by = "stable_assignment", cols = stability_color, pt.size = 0) + geom_boxplot(width=0.2, position = position_dodge(width = 0.9),outlier.shape = NA,) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(title = "Gene Count") + coord_flip())
dev.off()

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/VlnPlot_helm_batch1_13_rectum_percentmt_stable_unstable_9_23_25.pdf", width = 5, height = 12)
plot(VlnPlot(rectum_unstable_merge, features = 'percent.mt', group.by = "donor", split.by = "stable_assignment", cols = stability_color, pt.size = 0) + theme(legend.position = "none") + geom_boxplot(width=0.2, position = position_dodge(width = 0.9),outlier.shape = NA,) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(title = "Mitochondria Percent") + coord_flip())
dev.off()

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/VlnPlot_helm_batch1_13_rectum_ncount_rna_stable_unstable_9_23_25.pdf", width = 5, height = 12)
plot(VlnPlot(rectum_unstable_merge, features = 'nCount_RNA', group.by = "donor", split.by = "stable_assignment", cols = stability_color, pt.size = 0) + theme(legend.position = "none") + geom_boxplot(width=0.2, position = position_dodge(width = 0.9),outlier.shape = NA,) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(title = "UMI Count") + coord_flip())
dev.off()


######Supp. Figure 2e####
#-------stacked bar plot of stably assigned vs. unstably assigned cells----------#

#load in rPCA seurat object 

ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_rPCA_10_17_24.rds")

#ileum barcodes 
barcodes <- colnames(ileum)

#load in stably assigned cells 
filtered_values <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/filtered_low_res_barcodes_10_24_24.csv', what = "", sep = ",", skip = 1)


unstable <- setdiff(barcodes, filtered_values)

DimPlot(ileum, reduction = "umap", cells.highlight = unstable, sizes.highlight = 0.1) + scale_color_manual(labels = c("Stably Assigned Cells", "Unstably Assigned Cells"), values = c("grey", "red"))

#add meta data labeling stable and unstable cells 
ileum$stable_assignment <- ifelse(
  colnames(ileum) %in% unstable,
  "unstable",   # label for cells in the vector
  "stable"     # label for cells not in the vector
)

#add broad cell type labeling 
ileum@meta.data$cell_type_broad <- ileum@meta.data$seurat_clusters
ileum$cell_type_broad <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"),
  to = c("Naive B Cell", 'Stem/Paneth Cell', 'Enterocyte', 'CD4 T Cell', 'Plasma Cell', 'Macrophage', 'NK/CD8 T Cell', 'Goblet Cell', 'Monocyte', 'Memory B Cell', 'Mesenchymal Cell', "Endothelial Cell", "Mast Cell", "Tuft Cell", "Doublet Cell", "Plasmablast")
)

#create new label 
ileum$celltype_stable <- paste(ileum$cell_type_broad, ileum$stable_assignment, sep = "_")



pt <- table(ileum$cell_type_broad, ileum$stable_assignment)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "NK/CD8 T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem/Paneth Cell", "Enterocyte", "Goblet Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell", "Doublet Cell"))

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure8_7_30_25/StackedBarPlot_batch1_13_ileum_stable_unstable_8_6_25.pdf", width = 10, height = 8)
plot(ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
       scale_fill_manual(values = c("grey", "red")) +
       theme_bw(base_size=16) +
       geom_col(position = "fill", width = 0.5) +
       xlab("Broad Cell Type") +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
       theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
       ylab("Proportion") +
       
       #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
       theme(legend.title = element_blank()) +
       theme(legend.text=element_text(size=16)) +
       theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
       theme(axis.text.y = element_text(size=16)))
dev.off()

#-----------------Supplementary Figure 3----------------#
#####Supp. Figure 3a#######
#ileum high resolution stably assigned cell
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

#load in ileum high resolution stably assigned cell clustering results 
ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_13_ileum_75_15_2000_rpca_11_11_24.rds')

ileum@meta.data$cell_typev2 <- ileum@meta.data$seurat_clusters
ileum$cell_typev2 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36'),
  to = c("Naive B Cell", "Memory B Cell", "Enterocyte", "Stem/Paneth Cell", "CD4 T Cell", "Enterocyte", "Stem/Paneth Cell", 'CD8 T Cell', "Plasma Cell", "Enterocyte", "Monocyte", "Plasma Cell", "Enterocyte", "Goblet Cell", "Stem Cell", "Regulatory T Cell", "Ribosomal Cluster", "Goblet Cell", "Inflammatory Macrophage", "Macrophage", "Ambiguous T Cell", "Inflammatory Monocyte", "TA Cell", "NK/T Cell", "Mesenchymal Cell", "Germinal Center B Cell", "Endothelial Cell", "Mesenchymal Cell", "Mast Cell", "Tuft Cell", "Doublet Cluster", "Plasma Cell", "Naive B Cell", "ILC3", "Enteroendocrine", "Plasmablast", "Cycling NK/T Cell")
)

ileum$cell_typev2 <- factor(ileum$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Germinal Center B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/T Cell", "Cycling NK/T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Inflammatory Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Stem Cell", "TA Cell", "Enterocyte", "Goblet Cell", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))

ileum@meta.data$donor <- ileum@meta.data$orig.ident
ileum$donor <- plyr::mapvalues(
  x = ileum$orig.ident,
  from = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'),
  to = c('donor1', 'donor2', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor17', 'donor18', 'donor20', 'donor21', 'donor23', 'donor24', 'donor26', 'donor27', 'donor29', 'donor30', 'donor31', 'donor33', 'donor34')
)


DimPlot(ileum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[4], green[3], green[6], blue[1], blue[5], blue[6], blue[2], blue[7], blue[8], blue[3], orange[1], orange[5], orange[2], orange[3], orange[4], pink[1], pink[11], pink[6], pink[2], pink[3], pink[5], pink[10], purple[1], purple[2], red[1], red[2])) + theme(plot.title = element_blank())

#save as high resolution pdf 
pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure2_12_18_24/UMAP_batch1_13_ileum_high_res_cell_type_12_18_24.pdf", width = 14, height = 8)
plot(DimPlot(ileum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[4], green[3], green[6], blue[1], blue[5], blue[6], blue[2], blue[7], blue[8], blue[3], orange[1], orange[5], orange[2], orange[3], orange[4], pink[1], pink[11], pink[6], pink[2], pink[3], pink[5], pink[10], purple[1], purple[2], red[1], red[2])) + theme(plot.title = element_blank()) + theme(text = element_text(size = 16)) + theme(axis.text.x=element_text(size=16)) + theme(axis.text.y=element_text(size=16)))
dev.off()

ileum$donor <- factor(ileum$donor, levels = c('donor1', 'donor2', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor17', 'donor18', 'donor20', 'donor21', 'donor23', 'donor24', 'donor26', 'donor27', 'donor29', 'donor30', 'donor31', 'donor33', 'donor34'))

Idents(ileum) <- ileum$cell_typev2
pt <- table(Idents(ileum), ileum$donor)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Germinal Center B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/T Cell", "Cycling NK/T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Inflammatory Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Stem Cell", "TA Cell", "Enterocyte", "Goblet Cell", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure2_12_18_24/StackedBarPlot_batch1_13_ileum_high_res_cell_type_12_18_24.pdf", width = 14, height = 8)
plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[4], green[3], green[6], blue[1], blue[5], blue[6], blue[2], blue[7], blue[8], blue[3], orange[1], orange[5], orange[2], orange[3], orange[4], pink[1], pink[11], pink[6], pink[2], pink[3], pink[5], pink[10], purple[1], purple[2], red[1], red[2])
       )) +
  theme_bw(base_size=16) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  theme(axis.text.y = element_text(size=16))
dev.off()

#####Supp. Figure 3b#######

#ileum reference clustering results 

ileum <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/reference/seurat_object/helm_batch1_13_ileum_75_30_2000_ref_rpca_11_11_24.rds')

ileum@meta.data$donor <- ileum@meta.data$orig.ident
ileum$donor <- plyr::mapvalues(
  x = ileum$orig.ident,
  from = c("helm_sam1", "helm_sam4", "helm_sam7", "helm_sam24", "helm_sam29", "helm_sam32", "helm_sam35", "helm_sam44", "helm_sam47", "helm_sam50", "helm_sam53", "helm_sam56", "helm_sam59", "helm_sam62", "helm_sam77", "helm_sam83", "helm_sam88", "helm_sam91", "helm_sam98", "helm_sam101", "helm_sam112", 'helm_sam119', 'helm_sam134', 'helm_sam147', 'helm_sam153', 'helm_sam162', 'helm_sam171'),
  to = c('donor1', 'donor2', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor17', 'donor18', 'donor20', 'donor21', 'donor23', 'donor24', 'donor26', 'donor27', 'donor29', 'donor30', 'donor31', 'donor33', 'donor34')
)

ileum@meta.data$cell_typev2 <- ileum@meta.data$seurat_clusters
ileum$cell_typev2 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35'),
  to = c("Naive B Cell", "Enterocyte", 'Memory B Cell', 'Enterocyte', 'Stem/Paneth Cell', "Stem/Paneth Cell", 'CD4 T Cell', 'Macrophage', "Plasma Cell", 'Monocyte', 'Stem Cell', 'Plasma Cell', 'Goblet Cell', 'Enterocyte', 'Goblet Cell', 'Ribosomal Cluster', 'Regulatory T Cell', 'NK/CD8 T Cell', 'CD8 T Cell', 'Ambiguous T Cell', 'Inflammatory Monocyte', 'TA Cell', 'Cycling B Cell', 'CD4 T Cell', 'Mesenchymal Cell', 'CD8 T Cell', 'Mesenchymal Cell', 'Endothelial Cell', 'Tuft Cell', 'Mast Cell', 'Enterocyte', 'Naive B Cell', 'Doublet Cluster', 'Enteroendocrine', 'ILC3', 'Plasmablast')
)

ileum$cell_typev2 <- factor(ileum$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/CD8 T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Stem Cell", "TA Cell", "Enterocyte", "Goblet Cell", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))

DimPlot(ileum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[4], green[3], green[6], blue[1], blue[5], blue[6], blue[2], blue[8], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[11], pink[6], pink[2], pink[3], pink[5], pink[10], purple[1], purple[2], red[1], red[2]))

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure2_12_18_24/UMAP_batch1_13_ileum_ref_cell_type_12_18_24.pdf", width = 12, height = 8)
plot(DimPlot(ileum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[4], green[3], green[6], blue[1], blue[5], blue[6], blue[2], blue[8], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[11], pink[6], pink[2], pink[3], pink[5], pink[10], purple[1], purple[2], red[1], red[2])) + theme(plot.title = element_blank()) + theme(text = element_text(size = 16)) + theme(axis.text.x=element_text(size=16)) + theme(axis.text.y=element_text(size=16)))
dev.off()


ileum$donor <- factor(ileum$donor, levels = c('donor1', 'donor2', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor17', 'donor18', 'donor20', 'donor21', 'donor23', 'donor24', 'donor26', 'donor27', 'donor29', 'donor30', 'donor31', 'donor33', 'donor34'))

Idents(ileum) <- ileum$cell_typev2
pt <- table(Idents(ileum), ileum$donor)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "Regulatory T Cell", "CD8 T Cell", "NK/CD8 T Cell", "ILC3", "Ambiguous T Cell", "Macrophage", "Monocyte", "Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Stem Cell", "TA Cell", "Enterocyte", "Goblet Cell", "Tuft Cell", "Enteroendocrine", "Mesenchymal Cell", "Endothelial Cell", "Ribosomal Cluster", "Doublet Cluster"))

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure2_12_18_24/StackedBarPlot_batch1_13_ileum_ref_cell_type_12_18_24.pdf", width = 14, height = 8)
plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[4], green[3], green[6], blue[1], blue[5], blue[6], blue[2], blue[8], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[11], pink[6], pink[2], pink[3], pink[5], pink[10], purple[1], purple[2], red[1], red[2])
       )) +
  theme_bw(base_size=16) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  ylab("Proportion") +
  
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  theme(axis.text.y = element_text(size=16))
dev.off()

#####Supp. Figure 3c#######

DefaultAssay(ileum) <- "RNA"

#Dotplot of marker genes for cell annotation - ileum 
marker_genes <- c("MS4A1", "IGHD", 'CD27', 'MEF2B', "MZB1", "IGHG2", "XBP1", "HSP90B1", "CD4", "TRBC1", "CD8A", "NKG7", "IL7R", "CD3D", 'APOE', 'C1QB', 'S100A8', 'NAMPT', 'TIMP1', "IDO1", 'KIT', 'TPSB2', 'OLFM4', 'DEFA5', 'FABP1', 'APOA4', 'MUC2', 'TFF3', 'MKI67', 'TOP2A', 'LRMP', 'HPGDS', 'VIM', 'COL1A1', 'PECAM1', 'SPARCL1')

#flip axis 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure1_12_17_24/DotPlot_helm_batch1_13_ileum_marker_genes_12_18_24.pdf", width = 10, height = 14)
plot(DotPlot(ileum, features = marker_genes, group.by = "cell_typev2") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(legend.text=element_text(size=16)) + coord_flip())
dev.off()

#Dotplot of marker genes for cell annotation - colon
DefaultAssay(colon) <- "RNA"

marker_genes <- c("MS4A1", "IGHD", 'CD27', 'MEF2B', "CD79B", "POU2F2", "JCHAIN", "IGHA1", "MZB1", "IGHG2", "CD4", "TRBC1", 'CTLA4', 'IL2RA', "CD8A", "NKG7", 'APOE', 'C1QB', 'S100A8', 'NAMPT', 'KIT', 'TPSB2', 'OLFM4', 'TUBA1B', 'CA2', 'PHGR1', 'SLC26A3', 'CEACAM7', "BEST4", "OTOP2", 'MUC2', 'TFF3', 'LRMP', 'HPGDS', 'VIM', 'COL1A1', 'PECAM1', 'SPARCL1')

#flip axis 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure2_12_18_24/DotPlot_helm_batch1_13_colon_marker_genes_12_17_24.pdf", width = 10, height = 14)
plot(DotPlot(colon, features = marker_genes, group.by = "cell_typev1") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(legend.text=element_text(size=16)) + coord_flip())
dev.off()

#Dotplot of marker genes for cell annotation - rectum 
DefaultAssay(rectum) <- "RNA"

marker_genes <- c("MS4A1", "IGHD", 'CD27', 'MEF2B', "CD79B", "POU2F2", "MZB1", "IGHG2", "HLA-DPB1", "VPREB3", "CD4", "TRBC1", "CD8A", "NKG7", "GNLY", "GZMA", 'IL7R', 'RORC', "CD3D", "IL32", 'APOE', 'C1QB', 'S100A8', 'NAMPT', 'KIT', 'TPSB2', 'OLFM4', "LEFTY1", "TOP2A", 'TUBA1B', 'CA2', 'SELENBP1', 'SLC26A3', 'CEACAM7', 'MUC2', 'TFF3', 'LRMP', 'HPGDS', "CHGA", "NEUROD1", 'VIM', 'COL1A1', 'PECAM1', 'SPARCL1')

#flip axis 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure2_12_18_24/DotPlot_helm_batch1_13_rectum_marker_genes_flip_12_18_24.pdf", width = 10, height = 14)
plot(DotPlot(rectum, features = marker_genes, group.by = "cell_typev2") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(color="black", size=16), axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(legend.text=element_text(size=16)) + coord_flip())
dev.off()

#----------Supplementary Figure 4----------#

#####Supp. Figure 4a#######

#use ileum, colon, and rectum SO from supp. figure 1

#ileum
props <- getTransformedProps(ileum$cell_typev2, ileum$donor, transform="logit")

props_ileum <- props$Proportions
props_ileum <- as.data.frame(props_ileum)

ileum_wide <- props_ileum %>%
  select(clusters, sample, Freq) %>%
  pivot_wider(names_from = sample, values_from = Freq, values_fill = 0)

ileum_wide <- as.data.frame(ileum_wide)

rownames(ileum_wide) <- ileum_wide$clusters

ileum_wide <- ileum_wide[,-1]

# Define epi cell types
epi_cells <- c("Stem/Paneth Cell", "Enterocyte", "Goblet Cell", "TA Cell", "Tuft Cell")

# Sum rows for epi cells
epi_cells_sum <- colSums(ileum_wide[rownames(ileum_wide) %in% epi_cells, ])

# Add the new row
ileum_wide["Epi Cell", ] <- epi_cells_sum



#colon
props <- getTransformedProps(colon$cell_typev1, colon$donor, transform="logit")

props_colon <- props$Proportions
props_colon <- as.data.frame(props_colon)

colon_wide <- props_colon %>%
  select(clusters, sample, Freq) %>%
  pivot_wider(names_from = sample, values_from = Freq, values_fill = 0)

colon_wide <- as.data.frame(colon_wide)

rownames(colon_wide) <- colon_wide$clusters

colon_wide <- colon_wide[,-1]

# Define immune cell types
epi_cells <- c("Stem/TA Cell", "Colonocyte", "Mature Colonocyte", "BEST4 Colonocyte", "Goblet Cell", "Tuft Cell")

# Sum rows for epithelial cells
epi_cells_sum <- colSums(colon_wide[rownames(colon_wide) %in% epi_cells, ])

# Add the new row
colon_wide["Epi Cell", ] <- epi_cells_sum

#rectum
props <- getTransformedProps(rectum$cell_typev2, rectum$donor, transform="logit")

props_rectum <- props$Proportions
props_rectum <- as.data.frame(props_rectum)

rectum_wide <- props_rectum %>%
  select(clusters, sample, Freq) %>%
  pivot_wider(names_from = sample, values_from = Freq, values_fill = 0)

rectum_wide <- as.data.frame(rectum_wide)

rownames(rectum_wide) <- rectum_wide$clusters

rectum_wide <- rectum_wide[,-1]

# Define immune cell types
epi_cells <- c("Stem Cell", "TA Cell", "Colonocyte", "Mature Colonocyte", "Goblet Cell", "Tuft Cell", "Enteroendocrine Cell")

# Sum rows for immune cells
epi_cells_sum <- colSums(rectum_wide[rownames(rectum_wide) %in% epi_cells, ])

# Add the new row
rectum_wide["Epi Cell", ] <- epi_cells_sum


#----merge ileum and colon data-------#

# Subset the epithelial cell row for ileum and colon data
epithelial_ileum <- ileum_wide["Epi Cell", , drop = FALSE]
epithelial_colon <- colon_wide["Epi Cell", , drop = FALSE]

# Identify matching donor names between ileum and colon
matched_donors <- intersect(colnames(epithelial_ileum), colnames(epithelial_colon))

# Subset to matched donors
epithelial_ileum_matched <- as.numeric(epithelial_ileum[, matched_donors])
epithelial_colon_matched <- as.numeric(epithelial_colon[, matched_donors])

# Combine the data into a dataframe
ileum_colon <- data.frame(
  Donor = matched_donors,
  Ileum = epithelial_ileum_matched,
  Colon = epithelial_colon_matched
)

num_colors <- 26
palette_colors <- colorRampPalette(brewer.pal(8, "Set2"))(num_colors)

colors <- c("#FFC0CB", "#E6E6FA", "#FFDAB9", "#98FF98", "#B0E0E6", "#FFFFE0", "#F08080", "#C8A2C8", "#9FE2BF", "#89CFF0", "#87CEEB", "#FFF5E1", "#77DD77", "#F8C8D0", "#D8B7DD", "#FFCC99", "#A8E6CF", "#E0B0FF", "#66B2B2", "#CCCCFF", "#D3D3D3", "#F9E79F", "#D6A1C2", "#B2AC88", "#7FFFD4", "#FFDAB9")
# Scatter plot of ileum vs colon epithelial cell proportions

plot <- ggplot(ileum_colon, aes(x = Ileum, y = Colon)) +
  geom_point(size = 4, color = colors) + # Points
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  #geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +  # Linear regression line
  labs(
    x = "Ileum Proportion of Epithelial Cells",
    y = "Colon Proportion of Epithelial Cells"
  ) +
  theme(plot.background = element_rect(fill = "white")) +
  theme_minimal(base_size = 16)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure3_12_18_24/ScatterPlot_batch1_13_ileum_vs_colon_epi_prop_12_18_24.pdf", width = 5, height = 5)
plot(plot)
dev.off()

#"donor1"  "donor3"  "donor4"  "donor5"  "donor6"  "donor7"  "donor8" 
#"donor9"  "donor10" "donor11" "donor12" "donor13" "donor14" "donor17"
#"donor18" "donor20" "donor21" "donor23" "donor24" "donor26" "donor27"
#"donor29" "donor30" "donor31" "donor33" "donor34"

#-------merge colon and rectum data--------#

#add columns to rectum data to match colon data (15s1/s2, 16_s1/s2, and 28_s1/s2)
rectum_wide$donor15_s1 <- rectum_wide$donor15
rectum_wide$donor15_s2 <- rectum_wide$donor15
rectum_wide$donor16_s1 <- rectum_wide$donor16
rectum_wide$donor16_s2 <- rectum_wide$donor16
rectum_wide$donor28_s1 <- rectum_wide$donor28
rectum_wide$donor28_s2 <- rectum_wide$donor28


# Subset the epithelial cell row for ileum and colon data
epithelial_colon <- colon_wide["Epi Cell", , drop = FALSE]
epithelial_rectum <- rectum_wide["Epi Cell", , drop = FALSE]

# Identify matching donor names between ileum and rectum (same as colon)
matched_donors <- intersect(colnames(epithelial_colon), colnames(epithelial_rectum))

# Subset to matched donors
epithelial_colon_matched <- as.numeric(epithelial_colon[, matched_donors])
epithelial_rectum_matched <- as.numeric(epithelial_rectum[, matched_donors])

# Combine the data into a dataframe
colon_rectum <- data.frame(
  Donor = matched_donors,
  Colon = epithelial_colon_matched,
  Rectum = epithelial_rectum_matched
)

colors_c_r <- c("#FFC0CB", "#E6E6FA", "#FFDAB9", "#98FF98", "#B0E0E6", "#FFFFE0", "#F08080", "#C8A2C8", "#9FE2BF", "#89CFF0", "#87CEEB", "#FFF5E1", "#77DD77", "#FFA07A", "#FFA07A", "#D9A6D1", "#D9A6D1", "#F8C8D0", "#D8B7DD", "#A3D8A7", "#FFCC99", "#A8E6CF", "#B8A4D3", "#66B2B2", "#A2C2E8", "#CCCCFF", "#D3D3D3", "#AFEEEE", "#AFEEEE", "#F9E79F", "#D6A1C2", "#B2AC88", "#F6A8B8", "#7FFFD4", "#FFDAB9")


# Scatter plot of colon and rectal epithelial cell proportions

plot <- ggplot(colon_rectum, aes(x = Colon, y = Rectum)) +
  geom_point(size = 4, color = colors_c_r) + # Points
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  #geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +  # Linear regression line
  labs(
    x = "Colon Proportion of Epithelial Cells",
    y = "Rectum Proportion of Epithelial Cells"
  ) +
  theme(plot.background = element_rect(fill = "white")) +
  theme_minimal(base_size = 16)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure3_12_18_24/ScatterPlot_batch1_13_colon_vs_rectum_epi_prop_12_18_24.pdf", width = 5, height = 5)
plot(plot)
dev.off()

#---------merge colon and rectum data--------#

# Subset the epithelial cell row for ileum and rectum data
epithelial_ileum <- ileum_wide["Epi Cell", , drop = FALSE]
epithelial_rectum <- rectum_wide["Epi Cell", , drop = FALSE]

# Identify matching donor names between ileum and rectum (same as colon)
matched_donors <- intersect(colnames(epithelial_ileum), colnames(epithelial_rectum))

# Subset to matched donors
epithelial_ileum_matched <- as.numeric(epithelial_ileum[, matched_donors])
epithelial_rectum_matched <- as.numeric(epithelial_rectum[, matched_donors])

# Combine the data into a dataframe
ileum_rectum <- data.frame(
  Donor = matched_donors,
  Ileum = epithelial_ileum_matched,
  Rectum = epithelial_rectum_matched
)

num_colors <- 26
palette_colors <- colorRampPalette(brewer.pal(8, "Set2"))(num_colors)

colors_i_r <- c("#FFC0CB", "#E6E6FA", "#FFDAB9", "#98FF98", "#B0E0E6", "#FFFFE0", "#F08080", "#C8A2C8", "#9FE2BF", "#89CFF0", "#87CEEB", "#FFF5E1", "#77DD77", "#F8C8D0", "#D8B7DD", "#FFCC99", "#A8E6CF", "#E0B0FF", "#CCCCFF", "#D3D3D3", "#F9E79F", "#D6A1C2", "#B2AC88", "#7FFFD4", "#FFDAB9")
# Scatter plot of ileum vs colon epithelial cell proportions

plot <- ggplot(ileum_rectum, aes(x = Ileum, y = Rectum)) +
  geom_point(size = 4, color = colors_i_r) + # Points
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  #geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +  # Linear regression line
  labs(
    x = "Ileum Proportion of Epithelial Cells",
    y = "Rectum Proportion of Epithelial Cells"
  ) +
  theme(plot.background = element_rect(fill = "white")) +
  theme_minimal(base_size = 16)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure3_12_18_24/ScatterPlot_batch1_13_ileum_vs_rectum_epi_prop_12_18_24.pdf", width = 5, height = 5)
plot(plot)
dev.off()

#####Supp. Figure 4b#######

#ileum

ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_udpate_11_26_24.rds")

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

#####Supp. Figure 4c#######

#colon

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

#####Supp. Figure 4d#######

#rectum 

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



#test if group membership is statistically significant

#create contingency table 

group_counts <- matrix(
  c(3, 7, 17, 0, 0, 0, 0, 0, 0, 0,   # Ileum
    0, 0, 0, 5, 22, 1, 8, 0, 0, 0,   # Colon
    0, 0, 0, 0, 0, 0, 0, 13, 9, 10), # Rectum
  nrow = 3, byrow = TRUE
)

# Add row and column names
rownames(group_counts) <- c("Ileum", "Colon", "Rectum")
colnames(group_counts) <- paste0("Group ", 1:10)
group_counts

# Perform the chi-squared test
chisq_test <- chisq.test(group_counts)
chisq_test

#X-squared = 190, df = 18, p-value < 2.2e-16

# Perform Fisher's exact test
fisher_test <- fisher.test(group_counts)
fisher_test

#p-value < 2.2e-16
#alternative hypothesis: two.sided

#alt. groupings based on interpretation 

group_counts_alt <- matrix(
  c(3, 7, 17, 0, 0,    # Ileum
    1, 5, 22, 8, 0,   # Colon
    0, 9, 13, 0, 10), # Rectum
  nrow = 3, byrow = TRUE
)

# Add row and column names
rownames(group_counts_alt) <- c("Ileum", "Colon", "Rectum")
colnames(group_counts_alt) <- paste0("Group ", 1:5)
group_counts

# Perform the chi-squared test
chisq_test <- chisq.test(group_counts_alt)
chisq_test

#X-squared = 40.943, df = 8, p-value = 2.137e-06

# Perform Fisher's exact test
fisher_test <- fisher.test(group_counts_alt)
fisher_test

#p-value = 1.46e-06
#alternative hypothesis: two.sided


####Supp. Figure 4e####
monocyte_props <- colon@meta.data %>%
  dplyr::group_by(donor) %>%
  dplyr::mutate(total_cells = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(donor, cell_typev1, macro_IF) %>%
  dplyr::summarize(
    n_cells = n(),
    total_cells = unique(total_cells),
    proportion = n_cells / total_cells,
    .groups = "drop"
  ) %>%
  dplyr::filter(cell_typev1 == "Monocyte")


plot_s4e <- ggplot(monocyte_props, aes(x = macro_IF, y = proportion, fill = macro_IF)) +
  #geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot() +
  labs(
    fill = "Macro IF status",
    x = "Macro IF status",
    y = "Proportion of Monocytes"
  ) +
  scale_fill_manual(values = c("IF" = "#CC79A7", "NI" = "#009E73")) +
  theme_minimal(base_size = 16) 
#theme(plot.background = element_rect(fill = "white")) +
#theme_minimal(base_size = 16)


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure3_12_18_24/BoxPlot_colon_monocyte_prop_macro_IF_6_25_24.pdf", width = 5, height = 5)
plot(plot_s4e)
dev.off()

####Supp. Figure 4f####
#proportion of monocytes - micro IF status 
total_cells_df <- as_tibble(colon@meta.data) %>%
  group_by(donor) %>%
  summarize(total_cells = n(), .groups = "drop")

# Step 2: Get monocyte counts per donor and micro_IF
monocyte_props_micro <- colon@meta.data %>%
  as.data.frame() %>%
  dplyr::filter(cell_typev1 == "Monocyte") %>%
  dplyr::group_by(donor, micro_IF) %>%
  dplyr::summarize(n_cells = n(), .groups = "drop") %>%
  dplyr::left_join(total_cells_df, by = "donor") %>%
  dplyr::mutate(proportion = n_cells / total_cells)

plot_s4f <- ggplot(monocyte_props_micro, aes(x = micro_IF, y = proportion, fill = micro_IF)) +
  #geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot() +
  labs(
    fill = "Micro IF status",
    x = "Micro IF status",
    y = "Proportion of Monocytes"
  ) +
  scale_fill_manual(values = c("IF" = "#CC79A7", "NI" = "#009E73")) +
  theme_minimal(base_size = 16) 
#theme(plot.background = element_rect(fill = "white")) +
#theme_minimal(base_size = 16)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure3_12_18_24/BoxPlot_colon_monocyte_prop_micro_IF_7_30_24.pdf", width = 5, height = 5)
plot(plot_s4f)
dev.off()



#----------Supplementary Figure 5----------#

#####Supp. Figure 5a#######

#enterocyte volcano plot 

ent <- topTable(res.dl[["Enterocyte"]], coef = "groupgroup2", number = Inf)

ent$group <- 'NS'
ent$ID <- rownames(ent)
ent$group[ent$logFC < (0) & ent$adj.P.Val<0.05] <-'DOWN'
ent$group[ent$logFC > (0) & ent$adj.P.Val<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
ent$label <- ifelse(ent$group == "NS", NA, ent$group) #might have to do this 

ent$label[ent$group !='NS'] <- ent$ID[ent$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=ent,mapping=aes(x=logFC,y=-log10(adj.P.Val),color=group,label=label))+
  geom_point()+
  #xlim(-0.228,0.228) +
  #labs(title='volcano plot of significance against DEG recurrence and non-recurrence') +
  xlab("log2FC") +
  labs() +
  geom_label_repel(na.rm = TRUE, show.legend = F, box.padding = 0.5, size = 5) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16), legend.position = "right") +
  #geom_vline(xintercept=c(-1, 1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  labs(color = ent$group) + 
  scale_color_manual(values=c("DOWN"= "#F8766D", "NS"="#999999", "UP"="#619CFF"), name = "Group", labels = c("DOWN" = "group 1", "UP" = "group 2"))

volcano_gaussian

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure4_12_18_24/VolcanoPlot_helm_batch1_13_ileum_group_enterocyte_1_2_25.pdf", width = 7, height = 7)
plot(volcano_gaussian)
dev.off()


#####Supp. Figure 5b#######

ent_g2 <- rownames(ent_DEG_up)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = ent_g2, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_ent_g2 <- data.frame(ego)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure4_12_18_24/BarPlot_helm_batch1_13_ileum_group2_enterocyte_pathway_1_2_25.pdf", width = 10, height = 10)
barplot(ego, showCategory = 12, font.size = 16)
dev.off()


#####Supp. Figure 5c#######

cellchat_g1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/data/helm_batch1_13_ileum_group1_cellchat_obj_11_23_24.rds")

cellchat_g2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/data/helm_batch1_13_ileum_group2_cellchat_obj_11_23_24.rds")

object.list <- list(G1 = cellchat_g1, G2 = cellchat_g2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure4_12_18_24/BarPlot_helm_batch1_13_ileum_cellchat_pathways_1_2_25.pdf", width = 9, height = 10)
plot(rankNet(cellchat, mode = "comparison", measure = "weight", color.use = c('#F8766D', '#619CFF'), sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, font.size = 12))
dev.off()


#####Supp. Figure 5d#######

#metallopeptidase activity genes

ileum_donor <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_scITD_group_donor_avg_met_pep_counts.csv")

rownames(ileum_donor) <- ileum_donor$X

colnames(ileum_donor) <- c("genes", "donor1", "donor2", "donor3", "donor4", "donor5", "donor6", "donor7", "donor8", "donor9", "donor10", "donor11", "donor12", "donor14", "donor17", "donor18", "donor20", "donor21", "donor23", "donor24", "donor26", "donor27", "donor30", "donor34")

ileum_donor <- subset(ileum_donor, select = -genes)

# Example donor metadata
donor_metadata <- data.frame(
  donor = colnames(ileum_donor),
  group = c("group1", "group1", "group1", "group1", "group1", "group2", "group2", "group1", "group1", "group2", "group1", "group2", "group1", "group2", "group1", "group1", "group2", "group1", "group2", "group2", "group1", "group2", "group2") # Example groups
)

# Annotate donors
group_colors <- c("group1" = "#F8766D", "group2" = "#619CFF")
donor_annot <- HeatmapAnnotation(Group = donor_metadata$group,  col = list(Group = group_colors))


mat_scaled_sub <- t(scale(t(ileum_donor)))


plot_mp <- Heatmap(mat_scaled_sub, top_annotation = donor_annot, show_heatmap_legend = TRUE, cluster_columns = TRUE, name = "Z-score", row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 16), heatmap_legend_param = list(labels_gp = gpar(fontsize = 16)))

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure4_12_18_24/Heatmapt_batch1_13_ileum_group_met_pep_bias_gene.pdf", width = 10, height = 8)
plot(plot_mp)
dev.off()

#####Supp. Figure 5e#######

ileum_neut_chem <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_scITD_group_donor_avg_counts_neut_chem_bias.csv")

rownames(ileum_neut_chem) <- ileum_neut_chem$X

colnames(ileum_neut_chem) <- c("genes", "donor1", "donor2", "donor3", "donor4", "donor5", "donor6", "donor7", "donor8", "donor9", "donor10", "donor11", "donor12", "donor14", "donor17", "donor18", "donor20", "donor21", "donor23", "donor24", "donor26", "donor27", "donor30", "donor34")

ileum_neut_chem <- subset(ileum_neut_chem, select = -genes)

# Example donor metadata
donor_metadata <- data.frame(
  donor = colnames(ileum_donor),
  group = c("group1", "group1", "group1", "group1", "group1", "group2", "group2", "group1", "group1", "group2", "group1", "group2", "group1", "group2", "group1", "group1", "group2", "group1", "group2", "group2", "group1", "group2", "group2") # Example groups
)

mat_scaled_sub_neut <- t(scale(t(ileum_neut_chem)))


Heatmap(mat_scaled_sub_neut, top_annotation = donor_annot, show_heatmap_legend = TRUE, cluster_columns = TRUE, name = "Z-score")

plot_neut <- Heatmap(mat_scaled_sub_neut, top_annotation = donor_annot, show_heatmap_legend = TRUE, cluster_columns = TRUE, name = "Z-score", row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 16), heatmap_legend_param = list(labels_gp = gpar(fontsize = 16)))

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure4_12_18_24/Heatmapt_batch1_13_ileum_group_neut_chem_bias_gene_update.pdf", width = 10, height = 8)
plot(plot_neut)
dev.off()


#----------Supplementary Figure 6----------#

#####Supp. Figure 6a#######
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_scITD_group_UCell_score_12_7_24.rds")

colon$donor <- factor(colon$donor, levels = c("donor1", "donor3", "donor6", "donor7", "donor10", "donor12", "donor13", "donor15_s2", "donor16_s2", "donor17", "donor19", "donor20", "donor21", "donor23", "donor24", "donor26", "donor27", "donor32", "donor34", "donor11", "donor14", "donor15_s1", "donor16_s1", "donor18", "donor22", "donor25", "donor28_s1", "donor28_s2", "donor29", "donor30", "donor31", "donor33"))

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure5_12_18_24/VlnPlot_helm_batch1_13_colon_mac_activation_donor_12_18_24.pdf", width = 10, height = 8)
plot(VlnPlot(colon, features = "signature_1mac_sig", group.by = "donor", pt.size = 0, cols = c('#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#D65F5F', '#7986CB', '#7986CB', '#7986CB', '#7986CB', '#7986CB', '#7986CB', '#7986CB', '#7986CB', '#7986CB', '#7986CB', '#7986CB', '#7986CB', '#7986CB')) + theme(legend.position = "none") + theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)))
dev.off()

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure5_12_18_24/VlnPlot_helm_batch1_13_colon_mac_activation_celltype_12_18_24.pdf", width = 11, height = 8)
plot(VlnPlot(colon, features = "signature_1mac_sig", group.by = "cell_typev1", pt.size = 0,  split.by = "group", cols = c('#D65F5F', '#7986CB'), sort = "decreasing") + theme(legend.position = "none") + theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)))
dev.off()


#####Supp. Figure 6b#######

metadata <- pb@colData
metadata <- as.data.frame(metadata)

colon_counts <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/helm_batch1_13_colon_group_donor_avg_counts_12_7_24.csv")

rownames(colon_counts) <- colon_counts$X

# Remove "RNA" from column names
colnames(colon_counts) <- gsub("RNA.", "", colnames(colon_counts))

colon_counts$gene <- colon_counts$X

#1/14/2025 - update with relaxed threshold DEGs
mac_gene <- c("SLC11A1", "TLR2", "LRRK2", "C5AR1", "ITGAM", "CCL3", "FCGR3A", "CD93", "PTPRC", "TYROBP", "JAK2")


colon_counts <- colon_counts[,-1]

colon_filt <- colon_counts[colon_counts$gene %in% mac_gene,]

colon_filt <- colon_filt[,-33] #8 genes; update- 11 genes

colon_filt <- t(colon_filt)


#perform PCA 
colon_pca <- prcomp(colon_filt, scale = T)



all(row.names(colon_pca) == row.names(metadata)) #check that samples match ##TRUE

#variance explained by each PC
pc_eigenvalues <- colon_pca$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues #PC1: 76.2% of variance; PC2: 7.21 % of variance 

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- colon_pca$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

pc_scores

#add sample name column 
metadata$sample <- rownames(metadata)
pc_score_metadata <- full_join(pc_scores, metadata, by = ("sample"))

pc_score_metadata$group <- factor(pc_score_metadata$group, levels = c("group1", "group2"))


#relaxed threshold
plot <- pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(group))) +
  geom_point(size = 3) +
  theme_bw() +
  guides(color = guide_legend(title = "Group")) +
  theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  xlab("PC1 (76.2%)") +
  ylab("PC2 (7.21%)") +
  scale_color_manual(labels = c("Group 1", "Group 2"), values = c('#D65F5F', '#7986CB'))

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure5_12_18_24/PCA_batch1_13_colon_group_mac_act_DEG_all_cell_3_6_25.pdf", width = 6, height = 5)
plot(plot)
dev.off()


#relaxed threshold

pc_score_metadata$macro_IF <- factor(pc_score_metadata$macro_IF, levels = c("NI", "IF"))
plot <- pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(macro_IF))) +
  geom_point(size = 3) +
  theme_bw() +
  guides(color = guide_legend(title = "Macro IF Status")) +
  theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  xlab("PC1 (76.2%)") +
  ylab("PC2 (7.21%)") +
  scale_color_manual(labels = c("NI", "IF"), values = c('#0E4D92', '#A94064'))
plot

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure5_12_18_24/PCA_batch1_13_colon_group_mac_act_DEG_all_cell_macroIF_3_6_25.pdf", width = 6.5, height = 5)
plot(plot)
dev.off()


#####Supp. Figure 6c#######

#colon cellchat figure 

cellchat_g1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/colon/data/helm_batch1_13_colon_group1_cellchat_obj_12_5_24.rds")

cellchat_g2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/colon/data/helm_batch1_13_colon_group2_cellchat_obj_12_5_24.rds")

object.list <- list(G1 = cellchat_g1, G2 = cellchat_g2)

#load in cellchat object 
cellchat <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/colon/data/helm_batch1_13_colon_group_comb_cellchat_obj_12_5_24.rds")

#pathway figure
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure5_12_18_24/BarPlot_helm_batch1_13_colon_group_cellchat_pahtways_12_18_24.pdf", width = 10, height = 12)
plot(rankNet(cellchat, mode = "comparison", measure = "weight", color.use = c('#D65F5F', '#7986CB'), sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, font.size = 12))
dev.off()


#----------Supplementary Figure 7----------#

#####Supp. Figure 7a#######

#myeloid cell pathway enrichment in group 1

myeloid_DEG_up <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_myeloid_cell_batch_cor_group1_relax_12_10_24.csv")

myeloid_g1 <- myeloid_DEG_up$X

#boferroni - more stringent threshold 
ego <- enrichGO(gene = myeloid_g1, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_myl_g1 <- data.frame(ego)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure6_12_19_24/BarPlot_helm_batch1_13_rectum_group1_myeloid_cell_pathway_12_19_24.pdf", width = 10, height = 8)
barplot(ego, showCategory = 8, font.size = 16)
dev.off()

#####Supp. Figure 7b#######

#cell chat barplot

cellchat_g1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/data/helm_batch1_13_rectum_group1_cellchat_obj_12_11_24.rds")

cellchat_g2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/data/helm_batch1_13_rectum_group2_cellchat_obj_12_11_24.rds")

object.list <- list(G1 = cellchat_g1, G2 = cellchat_g2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure6_12_19_24/BarPlot_helm_batch1_13_rectum_group_cellchat_pathway_12_19_24.pdf", width = 10, height = 12)
plot(rankNet(cellchat, mode = "comparison", measure = "weight", color.use = c("#F06292", "#6A9FB5"), sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, font.size = 12))
dev.off()

#####Supp. Figure 7c#######

#TNF pathway enriched in group 2 chord plot 

#change colors to those used in UMAP 

colors_rectum <- c("Naive B Cell" = green[1], "Memory B Cell" = green[2], "Cycling B Cell" = green[4],"Plasma Cell" = green[3], "Ambig. B Cell" = green[7], "CD4 T Cell" = blue[1], "CD8 T Cell" = blue[6], "NK/CD8 T Cell" = blue[2], "ILC3" = blue[8], "Ambig. T Cell" = blue[3], "Macrophage" = orange[1], "Monocyte" = orange[2], "Mast Cell" = orange[4],"Stem Cell" =  pink[1], "Colonocyte" = pink[2], "Mature Colonocyte" = pink[10], "Goblet Cell" = pink[3], "TA Cell" = pink[6], "Tuft Cell" = pink[5], "Enteroendocrine Cell" = pink[12],"Mesenchymal Cell" = purple[1], "Endothelial Cell" = purple[2])

chordplot <- netVisual_chord_cell(cellchat_g2, lab.cex = 1, signaling = "TNF", color.use = colors_rectum, title.name = paste0("TNF signaling network"))

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure6_12_19_24/ChordPlot_helm_batch1_13_rectum_group2_cellchat_TNF_2_3_25.pdf", width = 8, height = 10)
netVisual_chord_cell(cellchat_g2, lab.cex = 1, signaling = "TNF", color.use = colors_rectum, title.name = paste0("TNF signaling network"))
dev.off()


#####Supp. Figure 7d#######

#IFN gamma module score 

#load in rectum with IFN and TNF response 
rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_group_scITD_mod_score_12_11_24.rds")

#re-order donors by group 
rectum$donor <- factor(rectum$donor, levels = c("donor4", "donor10", "donor14", "donor16", "donor18", "donor21", "donor22", "donor25", "donor28", "donor30", "donor31", "donor33", "donor34", "donor1", "donor5", "donor6", "donor7", "donor11", "donor12", "donor13", "donor15", "donor17", "donor19", "donor20", "donor26", "donor27", "donor29", "donor32"))

#group by donor
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure6_12_19_24/VlnPlot_helm_batch1_13_rectum_group_ifn_gamma_mod_score_donor_12_19_24.pdf", width = 10, height = 8)
plot(VlnPlot(rectum, features = "signature_1ifn_sig", group.by = "donor", pt.size = 0, cols = c("#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5")) + theme(legend.position = "none") + theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) + ggtitle("Response to interferon gamma module"))
dev.off()

#group by cell type 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure6_12_19_24/VlnPlot_helm_batch1_13_rectum_group_ifn_gamma_mod_score_celltype_12_19_24.pdf", width = 11, height = 8)
plot(VlnPlot(rectum, features = "signature_1ifn_sig", group.by = "cell_typev2", pt.size = 0, split.by = "group", cols = c('#F06292', '#6A9FB5'), sort = "decreasing") + theme(legend.position = "none") + theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) + ggtitle("Response to interferon gamma module"))
dev.off()


#####Supp. Figure 7e#######

#TNF module score 

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure6_12_19_24/VlnPlot_helm_batch1_13_rectum_group_tnf_mod_score_donor_12_19_24.pdf", width = 10, height = 8)
plot(VlnPlot(rectum, features = "signature_1tnf_sig", group.by = "donor", pt.size = 0, cols = c("#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#F06292", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5", "#6A9FB5")) + theme(legend.position = "none") + theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) + ggtitle("Response to TNF module"))
dev.off()

#group by cell type 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure6_12_19_24/VlnPlot_helm_batch1_13_rectum_group_tnf_mod_score_celltype_12_19_24.pdf", width = 11, height = 8)
plot(VlnPlot(rectum, features = "signature_1tnf_sig", group.by = "cell_typev2", pt.size = 0, split.by = "group", cols = c('#F06292', '#6A9FB5'), sort = "decreasing") + theme(legend.position = "none") + theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) + ggtitle("Response to TNFs module"))
dev.off()

#------------Supplementary Figure 8-------------#

#combined GWAS analysis
#read in marginal results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_gca_conti_spe_2025.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_marg = P) %>% 
  dplyr::arrange(P_marg)

n_sig <- ceiling(sqrt(nrow(df_cond))) #3
n_clusters <- nrow(df_marg)

#combine marginal and conditional results 
df_comb <- left_join(df_cond, df_marg, by = "VARIABLE")
df_comb <- df_comb %>% 
  mutate(VarCode = rep(c("a", "b"),times=nrow(df_comb)/2)) %>%
  pivot_wider(names_from = VarCode, values_from=c(VARIABLE, P_cond, P_marg)) %>%
  mutate(PS_a=log10(P_cond_a)/log10(P_marg_a), PS_b=log10(P_cond_b)/log10(P_marg_b))
# assert order
stopifnot(df_comb$P_marg_a <= df_comb$P_marg_b)

#reverse the order of forward selection if the criterion below is reached 
select_list <- df_marg$VARIABLE[1:n_sig]
reverse_list <- df_comb$VARIABLE_a[df_comb$PS_a < 0.2 & df_comb$PS_b >= 0.2]
select_list <- select_list[select_list %in% reverse_list == F][-1]

# The independent list starts with the most marginally significant cell type
indep_list <- c(toString(df_marg[1,1])) 
for (cur_var in select_list) {
  for (indep_var in indep_list) {
    # drop if the condition indicates so with any cell types in indep_list
    df_cur = df_comb[df_comb$VARIABLE_a==indep_var & df_comb$VARIABLE_b==cur_var,]
    if ((df_cur$PS_a >= 0.8 & df_cur$PS_b >= 0.8) | 
        (df_cur$PS_a >= 0.5 & df_cur$PS_b >= 0.5 & df_cur$P_cond_b < 0.001)) {
      if (tail(indep_list,1) == indep_var) indep_list <- c(indep_list, cur_var)
    } else break
  }
}
cat(indep_list)

#regulatory t cell and monocyte
bonferroni_threshold <- 0.05 / nrow(df_marg)

df_marg$cleaned_column <- gsub("^RNA\\.", "", df_marg$VARIABLE) # remove "RNA."
df_marg$cleaned_column <- gsub("\\.", " ", df_marg$cleaned_column)     # replace "." with space

df_marg$cleaned_column <- factor(df_marg$cleaned_column, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell","Plasma Cell", "IgA Plasma Cell", "IgG Plasma Cell", "Plasmablast", "Ambig  B Cell", "CD4 T Cell", 'CD8 T Cell', "NK CD8 T Cell", "ILC3", "Regulatory T Cell", "Ambig  T Cell", "Macrophage", "Monocyte", "Pro Inflammatory Monocyte", "Mast Cell", "Stem Cell", "Stem Paneth Cell", "Stem TA Cell", "TA Cell", "Enterocyte", "Colonocyte", "Mature Colonocyte", "BEST4 Colonocyte", "Goblet Cell", "Tuft Cell", "Enteroendocrine Cell", "Mesenchymal Cell", "Endothelial Cell"))
cols <- c(green[1], green[2], green[4], green[3], green[3], green[5], green[6], green[7], blue[1], blue[6], blue[2], blue[8], blue[5], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[1], pink[1], pink[6], pink[2], pink[2], pink[10], pink[4], pink[3], pink[5], pink[12], purple[1], purple[2])

sfig_8 <- df_marg %>% 
  # create the plot
  ggplot(aes(x = cleaned_column, y = -log10(P_marg), colour = factor(cleaned_column))) +
  geom_point(aes(size = -log10(P_marg))) +
  geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed", color = "red", size = 0.8) +
  theme_classic() +
  scale_size_area(max_size = 9) +
  #guides(color = guide_legend(title = "Cell type")) +
  theme(legend.position = "none", axis.text.x = element_text(color="black", size=16, angle = 60, hjust = 1),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  scale_color_manual(values = cols)+
  #scale_x_discrete(labels = c("RNA.Naive.B.Cell" = "Naive B Cell", "RNA.Memory.B.Cell" = "Memory B Cell", "Cycling_B_Cell" = "Cycling B Cell", "Plasma_Cell" = "Plasma Cell", "Ambig__B_Cell" = "Ambig. B Cell", "CD4_T_Cell" = "CD4 T Cell", "CD8_T_Cell" = 'CD8 T Cell', "NK_CD8_T_Cell" = "NK/CD8 T Cell", "ILC3" = "ILC3", "Ambig__T_Cell" = "Ambig. T Cell", "Ambig__T_Cell", "Macrophage" = "Macrophage", "Monocyte" = "Monocyte", "Mast_Cell" = "Mast Cell", "Stem_Cell" = "Stem Cell", "TA_Cell" = "TA Cell", "Colonocyte" = "Colonocyte", "Mature_Colonocyte" = "Mature Colonocyte", "Goblet_Cell" = "Goblet Cell", "Tuft_Cell" = "Tuft Cell", "Enteroendocrine_Cell" = "Enteroendocrine Cell", "Mesenchymal_Cell" = "Mesenchymal Cell", "Endothelial_Cell" = "Endothelial Cell")) +
  xlab("Cell types") +
  ylab("-log10(p value)")


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure9_7_30_25/ScatterPlot_all_topcell_bonferroni.pdf", width = 9, height = 7)
sfig_8
dev.off()


#----------Supp. Figure 9-------------#

#Venn diagram of group assignment 

library(ggVennDiagram)

#######Fig 9a########

#pro-inflammatory group 
G_proinflam <- list(Ileum = c("donor6", "donor7", "donor10", "donor12", "donor17", "donor21", "donor24", "donor26", "donor30", "donor34"), Colon = c("donor11", "donor14", "donor15", "donor16", "donor18", "donor22", "donor25", "donor28", "donor28", "donor29", "donor30", "donor31", "donor33"), Rectum = c("donor4", "donor10", "donor14", "donor16", "donor18", "donor21", "donor22", "donor25", "donor28", "donor30", "donor31", "donor33", "donor34"))

venn_proinflam <- ggVennDiagram(G_proinflam) + theme(legend.position = "none") + scale_fill_gradient(low = "white", high = "#C41E3A")

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure7_1_15_25/VennDiagram_pro_inflam_groups_1_15_25.pdf", width = 4, height = 4)
plot(venn_proinflam)
dev.off()

#######Fig 9b########
#other group

G_other <- list(Ileum = c("donor1", "donor2", "donor3", "donor4", "donor5", "donor8", "donor9", "donor11", "donor14", "donor18", "donor20", "donor23", "donor27"), Colon = c("donor1", "donor3", "donor6", "donor7", "donor10", "donor12", "donor13", "donor15_s2", "donor16_s2", "donor17", "donor19", "donor20", "donor21", "donor23", "donor24", "donor26", "donor27", "donor32", "donor34"), Rectum = c("donor1", "donor5", "donor6", "donor7", "donor11", "donor12", "donor13", "donor15", "donor17", "donor19", "donor20", "donor26", "donor27", "donor29", "donor32"))

venn_other <- ggVennDiagram(G_other) + theme(legend.position = "none") + scale_fill_gradient(low = "white", high = "#6082B6")

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/supplementary_figures/supp_figure7_1_15_25/VennDiagram_other_groups_1_15_25.pdf", width = 4, height = 4)
plot(venn_other)
dev.off()

