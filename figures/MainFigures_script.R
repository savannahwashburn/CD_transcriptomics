#Helmsley GCA Batch 1-13 Manuscript Figures

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

#---------------------Figure 1----------------------#
#12/13/2024

#####Figure 1a#######

#reference ileum clustering 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_rPCA_10_17_24.rds")

#save as png
png(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure1_12_13_24/UMAP_ileum_reference_12_13_24.png", width = 2500, height = 2000, res = 300)
plot(DimPlot(ileum, reduction = "umap", label = T, label.size = 5) + theme(text = element_text(size = 16)) + theme(axis.text.x=element_text(size=16)) + theme(axis.text.y=element_text(size=16))) 
dev.off()


#####Figure 1b#######

#use query 5 
cluster5 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/ileum_sub_res20_3_low_res_cluster_pred.rds")

cluster5$predicted.id <- factor(cluster5$predicted.id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

cluster5$query_label <- paste(cluster5$seurat_clusters, "query", sep = "_")
cluster5$ref_label <- paste(cluster5$predicted.id, "ref", sep = "_")


library(ArchR)
cm1 <- confusionMatrix(paste0(cluster5$ref_label), paste0(cluster5$query_label))
cm1 <- cm1 / Matrix::rowSums(cm1)

#change order of rows 
rownames(cm1)

new_order_row <- c("9_ref", "15_ref", "8_ref", "4_ref", "11_ref", "10_ref", "5_ref", "13_ref", "6_ref", "3_ref", "14_ref", "7_ref", "12_ref", "0_ref", "1_ref", "2_ref")

cm1 <- cm1[new_order_row, ]

#change order of columns 
colnames(cm1)

new_order_col <- c("10_query", "9_query", "4_query", "12_query", "11_query", "5_query", "13_query", "6_query", "3_query", "15_query", "7_query", "14_query", "0_query", "2_query", "8_query", "1_query")

cm1 <- cm1[ ,new_order_col]


cm1_plot <- pheatmap::pheatmap(
  mat = as.matrix(cm1), 
  color = paletteContinuous("whiteBlue"),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  display_numbers = FALSE,
  num_color = "white",
  fontsize = 16
)

#save as pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(cm1_plot, "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure1_12_13_24/ConfusionMatrixHeatmap_query5_cluster_predid_noclust_2_3_25.pdf")

#####Figure 1c#######

#highlight stably assigned and unstably assigned cells - use the ileum seurat object loaded in figure 1a
#load in filtered barcodes 

filt_barcodes <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/barcodes/filtered_low_res_barcodes_10_24_24.csv', what = "", sep = ",", skip = 1)

#ileum barcodes 
barcodes <- colnames(ileum)

unstable <- setdiff(barcodes, filt_barcodes)

png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure1_12_13_24/UMAP_ileum_highlight_unstable_12_13_24.png", width = 3000, height = 2000, res = 300)
plot(DimPlot(ileum, reduction = "umap", cells.highlight = unstable, sizes.highlight = 0.1) + scale_color_manual(labels = c("Stably assigned\ncells", "Unstably assigned\ncells"), values = c("grey", "red")) + theme(text = element_text(size = 18)) + theme(axis.text.x=element_text(size=18)) + theme(axis.text.y=element_text(size=18)))
dev.off()

#####Figure 1d#######

#jaccard index heatmap low vs. high resolution 

ileum_low <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/low_res/seurat_object/helm_batch1_13_ileum_50_15_2000_rpca_10_29_24.rds')

ileum_high <- readRDS(file = '/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/ileum/high_res/seurat_object/helm_batch1_13_ileum_75_15_2000_rpca_11_11_24.rds')


Idents(ileum_low) <- ileum_low$cell_typev2
Idents(ileum_high) <- ileum_high$cell_typev2

clusters_low <- unique(Idents(ileum_low))
clusters_high <- unique(Idents(ileum_high))


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
    low_cells <- WhichCells(ileum_low, idents = as.character(low_cluster))
    high_cells <- WhichCells(ileum_high, idents = as.character(high_cluster))
    
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

save_pheatmap_pdf <- function(x, filename, width=15, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(jaccard_plot, "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure1_12_13_24/HeatMap_jaccard_index_low_vs_high_res_12_13_24.pdf")

#####Figure 1e#######

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
                              color = colorRampPalette(c("purple", "white", "red"))(50), # Color gradient from white to blue
                              cluster_rows = TRUE,  # Disable clustering if you don't want it
                              cluster_cols = TRUE,
                              display_numbers = TRUE,
                              number_color = "black", 
                              fontsize = 16,
                              fontsize_number = 12)  # Optionally display the actual Jaccard index values

save_pheatmap_pdf <- function(x, filename, width=15.5, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(spearman_cor_plot, "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure1_12_13_24/HeatMap_spearman_cor_log2FC_low_vs_high_res_12_13_24.pdf")

#---------------------Figure 2----------------------#

#load ileum seurat object 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_11_19_24.rds")

#use cell type v2 for analysis 
ileum@meta.data$cell_typev2 <- ileum@meta.data$seurat_clusters
ileum$cell_typev2 <- plyr::mapvalues(
  x = ileum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21'),
  to = c("Naive B Cell", 'Stem/Paneth Cell', 'Enterocyte', 'CD4 T Cell', 'NK/CD8 T Cell', 'Enterocyte', 'Macrophage', "Plasma Cell", "Plasma Cell", 'Monocyte', 'Enterocyte', 'Goblet Cell', 'Goblet Cell', 'Ambig. T Cell', 'Mesenchymal Cell', "Pro-Inflammatory Monocyte", "TA Cell", "Memory B Cell", "Endothelial Cell", 'Mast Cell', "Tuft Cell", "Plasmablast")
)

ileum$cell_typev2 <- factor(ileum$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "NK/CD8 T Cell", "Ambig. T Cell", "Macrophage", "Monocyte", "Pro-Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Enterocyte", "Goblet Cell", "TA Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))


#Create UMAP 

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


DimPlot(ileum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[3], green[6], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[6], pink[5], purple[1], purple[2])) + theme(plot.title = element_blank())

#save high resolution png 
png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure2_12_13_24/UMAP_batch1_13_ileum_celltype_v2_12_13_24.png", width = 2700, height = 1800, res = 300)
plot(DimPlot(ileum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[3], green[6], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[6], pink[5], purple[1], purple[2])) + theme(plot.title = element_blank()) + theme(text = element_text(size = 16)) + theme(axis.text.x=element_text(size=16)) + theme(axis.text.y=element_text(size=16)))
dev.off()


#cell type v2 - group by hierarchical clusterin results
ileum$donor <- factor(ileum$donor, levels = c('donor29', 'donor13', 'donor31', 'donor33', 'donor11', 'donor17', 'donor23', 'donor18', 'donor4', 'donor9', 'donor1', 'donor8', 'donor20', 'donor6', 'donor7', 'donor27', 'donor5', 'donor12', 'donor26', 'donor30', 'donor10', 'donor21', 'donor24', 'donor14', 'donor34', 'donor2', 'donor3'))
Idents(ileum) <- ileum$cell_typev2
pt <- table(Idents(ileum), ileum$donor)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Plasma Cell", "Plasmablast", "CD4 T Cell", "NK/CD8 T Cell", "Ambig. T Cell", "Macrophage", "Monocyte", "Pro-Inflammatory Monocyte", "Mast Cell", "Stem/Paneth Cell", "Enterocyte", "Goblet Cell", "TA Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure2_12_13_24/StackedBarPlot_batch1_13_ileum_celltype_v2_hierarch_order_12_13_24.png", width = 3200, height = 2000, res = 300)
plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[6], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[6], pink[5], purple[1], purple[2])
       )) +
  theme_bw(base_size=10) +
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

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure2_12_13_24/StackedBarPlot_batch1_13_ileum_celltype_v2_hierarch_order_2_3_25.pdf", width = 12, height = 8)
plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[3], green[6], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[6], pink[5], purple[1], purple[2])
       )) +
  theme_bw(base_size=10) +
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

#---------load in colon seurat object---------#

#2/3/2025 - update with new colors 

#load in colon clustering annotation
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/seurat_object/helm_batch1_13_colon_25_30_2000_rpca_sub_cluster_11_4_24.rds")

colon$cell_typev1 <- plyr::mapvalues(
  x = colon$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18'),
  to = c('Colonocyte', 'Naive B Cell', 'Mature Colonocyte', 'IgA Plasma Cell', 'Stem/TA Cell', 'CD4 T Cell', 'Macrophage', 'NK/CD8 T Cell', 'Memory B Cell', 'Goblet Cell', 'Mesenchymal Cell', 'Monocyte', 'Regulatory T Cell', 'IgG Plasma Cell', 'Cycling B Cell', 'Tuft Cell', 'Mast Cell', 'Endothelial Cell', 'BEST4 Colonocyte')
)

colon$cell_typev1 <- factor(colon$cell_typev1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "IgA Plasma Cell", "IgG Plasma Cell", "CD4 T Cell", "Regulatory T Cell", "NK/CD8 T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem/TA Cell", "Colonocyte", "Mature Colonocyte", "BEST4 Colonocyte", "Goblet Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

#donor information
colon@meta.data$donor <- colon@meta.data$orig.ident
colon$donor <- plyr::mapvalues(
  x = colon$orig.ident,
  from = c("helm_sam2", "helm_sam8", "helm_sam26", "helm_sam31", "helm_sam33", "helm_sam36", "helm_sam45", "helm_sam48", "helm_sam52", "helm_sam54", "helm_sam57", "helm_sam60", "helm_sam63", "helm_sam65", "helm_sam66", "helm_sam71", "helm_sam72", "helm_sam78", "helm_sam84", "helm_sam86", "helm_sam89", "helm_sam92", "helm_sam96", "helm_sam99", "helm_sam102", "helm_sam110", "helm_sam113", "helm_sam120", "helm_sam128", "helm_sam129", "helm_sam135", "helm_sam148", "helm_sam154", "helm_sam157", "helm_sam163", "helm_sam172"),
  to = c('donor1', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8', 'donor9', 'donor10', 'donor11', 'donor12', 'donor13', 'donor14', 'donor15_s1', 'donor15_s2', 'donor16_s1', 'donor16_s2', 'donor17', 'donor18', 'donor19', 'donor20', 'donor21', 'donor22', 'donor23', 'donor24', 'donor25', 'donor26', 'donor27', 'donor28_s1', 'donor28_s2', 'donor29', 'donor30', 'donor31', 'donor32', 'donor33', 'donor34')
)


DimPlot(colon, reduction = "umap", group.by = "cell_typev1", cols = c(green[1], green[2], green[4], green[3], green[5], blue[1], blue[5], blue[2], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink [4], pink[3], pink[5], purple[1], purple[2]), raster = FALSE) + theme(plot.title = element_blank())

#save high resolution png 
png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure2_12_13_24/UMAP_batch1_13_colon_celltype_v1_2_3_25.png", width = 2500, height = 1800, res = 300)
plot(DimPlot(colon, reduction = "umap", group.by = "cell_typev1", cols = c(green[1], green[2], green[4], green[3], green[5], blue[1], blue[5], blue[2], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink [4], pink[3], pink[5], purple[1], purple[2]), raster = FALSE) + theme(plot.title = element_blank()) + theme(text = element_text(size = 16)) + theme(axis.text.x=element_text(size=16)) + theme(axis.text.y=element_text(size=16)))
dev.off()


#cell type v2 - group by hierarchical clustering results
colon$donor <- factor(colon$donor, levels = c('donor1', 'donor3', 'donor8', "donor16_s2", 'donor7', 'donor11', 'donor25', 'donor13', 'donor32', 'donor19', 'donor5', 'donor26', 'donor33', 'donor12', 'donor15_s2', 'donor24', 'donor21', 'donor22', 'donor27', 'donor6', 'donor17', 'donor28_s2', 'donor9', 'donor28_s1', 'donor10', 'donor23', 'donor34', 'donor4', 'donor15_s1', 'donor30', 'donor16_s1', 'donor20', 'donor14','donor29', 'donor18', 'donor31'))

Idents(colon) <- colon$cell_typev1
pt <- table(Idents(colon), colon$donor)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "IgA Plasma Cell", "IgG Plasma Cell", "CD4 T Cell", "Regulatory T Cell", "NK/CD8 T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem/TA Cell", "Colonocyte", "Mature Colonocyte", "BEST4 Colonocyte", "Goblet Cell", "Tuft Cell", "Mesenchymal Cell", "Endothelial Cell"))

png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure2_12_13_24/StackedBarPlot_batch1_13_colon_celltype_v1_hierarch_cluster_group_2_3_25.png", width = 3400, height = 2200, res = 300)
plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[4], green[3], green[5], blue[1], blue[5], blue[2], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink [4], pink[3], pink[5], purple[1], purple[2])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  ylab("Proportion") +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  theme(axis.text.y = element_text(size=16))
dev.off()

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure2_12_13_24/StackedBarPlot_batch1_13_colon_celltype_v1_hierarch_cluster_group_2_3_25.pdf", width = 12, height = 8)
plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[4], green[3], green[5], blue[1], blue[5], blue[2], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink [4], pink[3], pink[5], purple[1], purple[2])
       )) +
  theme_bw(base_size=10) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  ylab("Proportion") +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  theme(axis.text.y = element_text(size=16))
dev.off()

#-------load in rectum seurat object-------#
rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_cell_annotation_12_8_24.rds")

rectum@meta.data$cell_typev2 <- rectum@meta.data$seurat_clusters
rectum$cell_typev2 <- plyr::mapvalues(
  x = rectum$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27'),
  to = c('Colonocyte', 'Colonocyte', 'Naive B Cell', 'Stem Cell', 'TA Cell', 'Colonocyte', 'Mature Colonocyte', 'Memory B Cell', 'CD4 T Cell', 'Goblet Cell', 'Plasma Cell', 'Macrophage', 'Mature Colonocyte', 'Cycling B Cell', 'CD8 T Cell', 'Goblet Cell', 'Plasma Cell', 'Mature Colonocyte', 'Ambig. T Cell', 'NK/CD8 T Cell', 'Enteroendocrine Cell', 'Tuft Cell', 'Mesenchymal Cell', 'Monocyte', 'Mast Cell', 'Endothelial Cell', 'Ambig. B Cell', 'ILC3')
)

rectum$cell_typev2 <- factor(rectum$cell_typev2, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell", "Ambig. B Cell", "CD4 T Cell", "CD8 T Cell", "NK/CD8 T Cell", "ILC3", "Ambig. T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem Cell", "Colonocyte", "Mature Colonocyte", "Goblet Cell", "TA Cell", "Tuft Cell", "Enteroendocrine Cell", "Mesenchymal Cell", "Endothelial Cell"))

DimPlot(rectum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[4], green[3], green[7], blue[1], blue[6], blue[2], blue[8], blue[3], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink[3], pink[6], pink[5], pink[12], purple[1], purple[2]), raster = FALSE) + theme(plot.title = element_blank())

#save high resolution png 
png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure2_12_13_24/UMAP_batch1_13_rectum_celltype_v2_2_3_25.png", width = 3000, height = 1800, res = 300)
plot(DimPlot(rectum, reduction = "umap", group.by = "cell_typev2", cols = c(green[1], green[2], green[4], green[3], green[7], blue[1], blue[6], blue[2], blue[8], blue[3], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink[3], pink[6], pink[5], pink[12], purple[1], purple[2]), raster = FALSE) + theme(plot.title = element_blank()) + theme(text = element_text(size = 16)) + theme(axis.text.x=element_text(size=16)) + theme(axis.text.y=element_text(size=16)))
dev.off()


#cell type v2 - group by hierarchical clustering results
rectum$donor <- factor(rectum$donor, levels = c('donor33', 'donor17', 'donor26', "donor29", 'donor11', 'donor12', 'donor15', 'donor31', 'donor4', 'donor30', 'donor13', 'donor5', 'donor20', 'donor24', 'donor28', 'donor16', 'donor22', 'donor10', 'donor34', 'donor25', 'donor18', 'donor19', 'donor8', 'donor1', 'donor3', 'donor9', 'donor14', 'donor6', 'donor32', 'donor21', 'donor7', 'donor27'))

Idents(rectum) <- rectum$cell_typev2
pt <- table(Idents(rectum), rectum$donor)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("Naive B Cell", "Memory B Cell", "Cycling B Cell", "Plasma Cell", "Ambig. B Cell", "CD4 T Cell", "CD8 T Cell", "NK/CD8 T Cell", "ILC3", "Ambig. T Cell", "Macrophage", "Monocyte", "Mast Cell", "Stem Cell", "Colonocyte", "Mature Colonocyte", "Goblet Cell", "TA Cell", "Tuft Cell", "Enteroendocrine Cell", "Mesenchymal Cell", "Endothelial Cell"))

#save as pdf to edit the legend (one column) 
pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure2_12_13_24/StackedBarPlot_batch1_13_rectum_celltype_v2_hierarch_cluster_group_2_3_25.pdf", width = 14, height = 8)
plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[4], green[3], green[7], blue[1], blue[6], blue[2], blue[8], blue[3], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink[3], pink[6], pink[5], pink[12], purple[1], purple[2])
       )) +
  theme_bw(base_size=16) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  ylab("Proportion") +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  theme(axis.text.y = element_text(size=16))
dev.off()

png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure2_12_13_24/StackedBarPlot_batch1_13_rectum_celltype_v2_hierarch_cluster_group_2_3_25.png", width = 4000, height = 2000, res = 300)
plot(ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
       scale_fill_manual(values = c(green[1], green[2], green[4], green[3], green[7], blue[1], blue[6], blue[2], blue[8], blue[3], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink[3], pink[6], pink[5], pink[12], purple[1], purple[2])
       )) +
  theme_bw(base_size=16) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  ylab("Proportion") +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  #scale_fill_discrete(breaks=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  theme(axis.text.y = element_text(size=16))
dev.off()

#------------------Figure 3------------------#

#####Figure 3a#######

#ileum donor matrix heatmap 

ileum_cor <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/ileum/data/helm_batch1_13_ileum_batch_cor_container_11_20_24.rds")

ileum_metadata <- ileum@meta.data

ileum_metadata <- ileum_metadata[, -c(2:4, 7, 8, 13:15, 17)]

ileum_metadata <- unique(ileum_metadata)

# Define custom colors
batch_colors <- c("batch1" = "#648FFF", "batch3" = "#785EF0", "batch4" = "#DC267F", "batch5" = "#FE6100", "batch6" = "#FFB000", "batch7" = "#BBCDF1", "batch8" = "#1CF815", "batch9" = "#9B97E1", "batch10" = "#51CFC6", "batch11" = "#2319A5", "batch12" = "#B34C90", "batch13" = "#CC6677")
macroIF_colors <- c("IF" = "#E66100", "NI" = "#5D3A9B")
microIF_colors <- c("IF" = "#F88379", "NI" = "#AA336A", "N/A" = "#a3bac3")
sex_colors <- c("F" = "#dab1da", "M" = "yellow")
ancestry_colors <- c("EA" = "#c9def4", "AA" = "#f5ccd4", "Asian" = "#b8a4c9")



# Create annotation object
ha <- HeatmapAnnotation(
  batch = ileum_metadata$batch,
  macro_IF = ileum_metadata$macro_IF,
  micro_IF = ileum_metadata$micro_IF,
  sex = ileum_metadata$sex,
  ancestry = ileum_metadata$ancestry,
  col = list(batch = batch_colors, macro_IF = macroIF_colors, micro_IF = microIF_colors, sex = sex_colors, ancestry = ancestry_colors)
)

# get donor scores-metadata associations
ileum_cor <- get_meta_associations(ileum_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='rsq')

# plot donor scores
ileum_cor <- plot_donor_matrix(ileum_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                               show_donor_ids = TRUE,
                               add_meta_associations='rsq') 

# show the donor scores heatmap
ileum_rsq <- ileum_cor$plots$donor_matrix
Heatmap(ileum_rsq, top_annotation = ha)


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure3_12_14_24/DonorMatrixHeatmap_batch1_13_ileum_rsq_12_14_24.pdf", width = 10, height = 7)
draw(ileum_rsq)
dev.off()


ileum_cor <- get_meta_associations(ileum_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='pval')


# plot donor scores
ileum_cor <- plot_donor_matrix(ileum_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                               show_donor_ids = TRUE,
                               add_meta_associations='pval')

# show the donor scores heatmap
ileum_pval <- ileum_cor$plots$donor_matrix

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure3_12_14_24/DonorMatrixHeatmap_batch1_13_ileum_pval_12_14_24.pdf", width = 10, height = 7)
draw(ileum_pval)
dev.off()

#####Figure 3b#######


batch_colors <- c("batch1" = "#648FFF", "batch3" = "#785EF0", "batch4" = "#DC267F", "batch5" = "#FE6100", "batch6" = "#FFB000", "batch7" = "#BBCDF1", "batch8" = "#1CF815", "batch9" = "#9B97E1", "batch10" = "#51CFC6", "batch11" = "#2319A5", "batch12" = "#B34C90", "batch13" = "#CC6677")
macroIF_colors <- c("IF" = "#E66100", "NI" = "#5D3A9B")
microIF_colors <- c("IF" = "#F88379", "NI" = "#AA336A", "N/A" = "#a3bac3")
sex_colors <- c("F" = "#dab1da", "M" = "yellow")
ancestry_colors <- c("EA" = "#c9def4", "AA" = "#f5ccd4", "Asian" = "#b8a4c9")


colon_cor <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/colon/data/helm_batch1_13_colon_batch_cor_container_11_27_24.rds")

# get donor scores-metadata associations
colon_cor <- get_meta_associations(colon_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='rsq')

# plot donor scores
colon_cor <- plot_donor_matrix(colon_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                               show_donor_ids = TRUE,
                               add_meta_associations='rsq')

# show the donor scores heatmap
colon_rsq <- colon_cor$plots$donor_matrix

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure3_12_14_24/DonorMatrixHeatmap_batch1_13_colon_rsq_12_14_24.pdf", width = 12, height = 9)
draw(colon_rsq)
dev.off()


colon_cor <- get_meta_associations(colon_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='pval')


# plot donor scores
colon_cor <- plot_donor_matrix(colon_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                               show_donor_ids = TRUE,
                               add_meta_associations='pval')

# show the donor scores heatmap
colon_pval <- colon_cor$plots$donor_matrix

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure3_12_14_24/DonorMatrixHeatmap_batch1_13_colon_pval_12_14_24.pdf", width = 12, height = 9)
draw(colon_pval)
dev.off()

#####Figure 3c#######


batch_colors <- c("batch1" = "#648FFF", "batch3" = "#785EF0", "batch4" = "#DC267F", "batch5" = "#FE6100", "batch6" = "#FFB000", "batch7" = "#BBCDF1", "batch8" = "#1CF815", "batch9" = "#9B97E1", "batch10" = "#51CFC6", "batch11" = "#2319A5", "batch12" = "#B34C90", "batch13" = "#CC6677")
macroIF_colors <- c("IF" = "#E66100", "NI" = "#5D3A9B")
microIF_colors <- c("IF" = "#F88379", "NI" = "#AA336A", "N/A" = "#a3bac3")
sex_colors <- c("F" = "#dab1da", "M" = "yellow")
ancestry_colors <- c("EA" = "#c9def4", "AA" = "#f5ccd4", "Asian" = "#b8a4c9")


# get donor scores-metadata associations
rectum_cor <- get_meta_associations(rectum_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='rsq')

# plot donor scores
rectum_cor <- plot_donor_matrix(rectum_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                                show_donor_ids = TRUE,
                                add_meta_associations='rsq')

# show the donor scores heatmap
rectum_rsq <- rectum_cor$plots$donor_matrix

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure3_12_14_24/DonorMatrixHeatmap_batch1_13_rectum_rsq_12_14_24.pdf", width = 12, height = 9)
draw(rectum_rsq)
dev.off()


rectum_cor <- get_meta_associations(rectum_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='pval')


# plot donor scores
rectum_cor <- plot_donor_matrix(rectum_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                                show_donor_ids = TRUE,
                                add_meta_associations='pval')

# show the donor scores heatmap
rectum_pval <- rectum_cor$plots$donor_matrix

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure3_12_14_24/DonorMatrixHeatmap_batch1_13_rectum_pval_12_14_24.pdf", width = 12, height = 9)
draw(rectum_pval)
dev.off()

#------------------Figure 4------------------#

#####Figure 4a#######

#pathway bias results 

#---read in counts----# 
ileum_counts <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_group_avg_log2_counts.csv")

rownames(ileum_counts) <- ileum_counts$X

colnames(ileum_counts) <- c("genes", "group_1", "group_2")

#calculate log2FC
ileum_counts$log2FC <- ileum_counts$group_1 - ileum_counts$group_2

#positive is up in group 1
#negative is up in group 2

# Assuming your dataframe is called df and the column you are checking is called column_name
#ileum_counts$sig <- ifelse(ileum_counts$log2FC > 0, 1, -1)

ileum_counts$log2FC <- as.numeric(ileum_counts$log2FC)
ileum_counts$group_1 <- as.numeric(ileum_counts$group_1)
ileum_counts$group_2 <- as.numeric(ileum_counts$group_2)


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
table(overlap$sig_strict) #4, 19

#response to TNF
go_term <- "GO:0034612"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0034612 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db)


intersect(ileum_counts$genes, genes_GO0034612$SYMBOL) #235

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0034612$SYMBOL,]
table(overlap$sig_strict) #18,15

#cytokine activity 

#GO:0005125

go_term <- "GO:0005125"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0005125 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db)


intersect(ileum_counts$genes, genes_GO0005125$SYMBOL) #203

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0005125$SYMBOL,]
table(overlap$sig_strict) #15,6

#response to LPS 

#GO:0032496

go_term <- "GO:0032496"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0032496 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db)


intersect(ileum_counts$genes, genes_GO0032496$SYMBOL) #315

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0032496$SYMBOL,]
table(overlap$sig_strict) #21,15 

#collagen binding 

#GO:0005518

go_term <- "GO:0005518"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0005518 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #339


intersect(ileum_counts$genes, genes_GO0005518$SYMBOL) #67

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0005518$SYMBOL,]
table(overlap$sig_strict) #5,2

#inflammatory response 

go_term <- "GO:0006954"

# Convert GO term to gene symbols
genes_GO0006954 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #339


intersect(ileum_counts$genes, genes_GO0006954$SYMBOL) #752

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0006954$SYMBOL,]
table(overlap$sig_strict) #50,32

#collagen catabolic process 

go_term <- "GO:0030574"

genes_GO0030574 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #339


intersect(ileum_counts$genes, genes_GO0030574$SYMBOL) #43

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0030574$SYMBOL,]
table(overlap$sig_strict) #4,4

#neutrophil chemotaxis 

#GO:0030593

go_term <- "GO:0030593"

genes_GO0030593 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #106


intersect(ileum_counts$genes, genes_GO0030593$SYMBOL) #102

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0030593$SYMBOL,]
table(overlap$sig_strict) #21,6

#neutrophil migration 

#GO:1990266

go_term <- "GO:1990266"

genes_GO1990266 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #130


intersect(ileum_counts$genes, genes_GO1990266$SYMBOL) #124

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO1990266$SYMBOL,]
table(overlap$sig_strict) #21,6

#chemokine receptor binding 

#GO:0042379

go_term <- "GO:0042379"

genes_GO0042379 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #71


intersect(ileum_counts$genes, genes_GO0042379$SYMBOL) #59

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0042379$SYMBOL,]
table(overlap$sig_strict) #8,1

#response to bacterium 

#GO:0009617

go_term <- "GO:0009617"

genes_GO0009617 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #786


intersect(ileum_counts$genes, genes_GO0009617$SYMBOL) #668

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0009617$SYMBOL,]
table(overlap$sig_strict) #63, 39

#response to wounding 

#GO:0009611

go_term <- "GO:0009611"

genes_GO0009611 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #577


intersect(ileum_counts$genes, genes_GO0009611$SYMBOL) #549

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0009611$SYMBOL,]
table(overlap$sig_strict) #23,40

#cytokine receptor binding 

#GO:0005126

go_term <- "GO:0005126"

genes_GO0005126 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #273


intersect(ileum_counts$genes, genes_GO0005126$SYMBOL) #235

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0005126$SYMBOL,]
table(overlap$sig_strict) #13,6

#ECM disassembly 

#GO:0022617

go_term <- "GO:0022617"

genes_GO0022617 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #62


intersect(ileum_counts$genes, genes_GO0022617$SYMBOL) #60

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0022617$SYMBOL,]
table(overlap$sig_strict) #5,3


#growth factor binding 

#GO:0019838

go_term <- "GO:0019838"

genes_GO0019838 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #132


intersect(ileum_counts$genes, genes_GO0019838$SYMBOL) #128

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0019838$SYMBOL,]
table(overlap$sig_strict) #7,3


#ecm organization 

#GO:0030198

go_term <- "GO:0030198"

genes_GO0030198 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #318


intersect(ileum_counts$genes, genes_GO0030198$SYMBOL) #303

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0030198$SYMBOL,]
table(overlap$sig_strict) #5,12


#fibronectin binding 

#GO:0001968

go_term <- "GO:0001968"

genes_GO0001968 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db) #31


intersect(ileum_counts$genes, genes_GO0001968$SYMBOL) #29

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0001968$SYMBOL,]
table(overlap$sig_strict) #1,2

#ECM structural constituent

#GO:0005201

go_term <- "GO:0005201"  # Replace with your GO term ID

# Convert GO term to gene symbols
genes_GO0005201 <- bitr(get(go_term, org.Hs.egGO2ALLEGS), fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Hs.eg.db)


intersect(ileum_counts$genes, genes_GO0005201$SYMBOL) #159

overlap <- ileum_counts[ileum_counts$genes %in% genes_GO0005201$SYMBOL,]
table(overlap$sig) #2,5

#-------read in pathways for bar chart--------#

#load in file that has pathways and number of genes in each group 
bar_plot_pathways <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_group_pathway_enrichment_log2.csv")

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


#number of genes were obtained on 12/16/2024

plot <- test %>%
  ggplot(aes(x = factor(pathways, levels = c("Chemokine Receptor Binding (71)", "Neutrophil Chemotaxis (106)", "Neutrophil Migration (130)", "Cytokine Activity (237)", "Collagen Binding (69)", "Growth Factor Binding (132)", "Cytokine Receptor Binding (273)", "ECM Disassembly (62)", "Response to Bacterium (786)", "Inflammatory Response (848)", "Response to LPS (339)", "Response to TNF (254)", "Collagen Catabolic Process (45)", "Response to Wounding (577)", "Fibronectin Binding (31)", "ECM Organization (318)", "ECM Structural Constituent (173)", "Metallopeptidase activity (184)")), y = PercentGenes, fill = group))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c('#F8766D', '#619CFF')) +
  coord_flip()+
  scale_y_continuous(breaks = breaks_values,
                     labels = abs(breaks_values)) + 
  geom_col(position = "dodge") +
  xlab("Pathways") + ylab("Percent of Genes in Pathway") +
  theme_minimal()+
  theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
  guides(fill=guide_legend(title="Groups")) 
#geom_text(aes(label=abs(number))) 
plot

#save as png 
png(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure4_12_16_24/BarPlot_helm_batch1_13_ileum_group_bias_log2.png", width = 3000, height = 2000, res = 300)
plot(plot)
dev.off()


#save as pdf 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure4_12_16_24/BarPlot_helm_batch1_13_ileum_group_bias_log2.pdf", width = 10, height = 8)
plot(plot)
dev.off()


#####Figure 4b#######

#--------metallopeptidase activity module score---------#

#load in seurat object with module score for metallopeptidase activity 
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_group_mod_score.rds")

plot1 <- VlnPlot(ileum, features = "signature_1met_pep_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") + theme(axis.title.x=element_blank())
png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure4_12_16_24/VlnPlot_helm_batch1_13_ileum_group_met_pep_module.png", width = 2800, height = 1800, res = 300)
plot(plot1 + ggtitle("Metallopeptidase activity") +
       ylab("Module score") + theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)))
dev.off()

#save as pdf 
pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure4_12_16_24/VlnPlot_helm_batch1_13_ileum_group_met_pep_module.pdf", width = 10, height = 8)
plot(plot1 + ggtitle("Metallopeptidase activity") +
       ylab("Module score") + theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)))
dev.off()


#####Figure 4c#######

#use ileumfrom figure 4b

#chemokine receptor binding 

chem_genes <- intersect(ileum_counts$genes, genes_GO0042379$SYMBOL) #59

chem_genes <- paste(shQuote(chem_genes), collapse=", ")

chem_genes <- c('JAK1', 'S100A14', 'NES', 'XCL2', 'XCL1', 'CNIH4', 'STAT1', 'CCL20', 'CX3CR1', 'CCR2', 'CCRL2', 'CXCL8', 'CXCL6', 'PF4V1', 'CXCL1', 'PF4', 'PPBP', 'CXCL5', 'CXCL3', 'CXCL2', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13', 'CCL28', 'CXCL14', 'CCL26', 'CCL24', 'DEFB1', 'DEFB4A', 'CCL27', 'CCL19', 'CCL21', 'MSMP', 'C5', 'CXCL12', 'CCL22', 'CX3CL1', 'CCL17', 'CKLF', 'CXCL16', 'CCL2', 'CCL7', 'CCL11', 'CCL8', 'CCL13', 'CCL1', 'CCL5', 'CCL16', 'CCL14', 'CCL15', 'CCL23', 'CCL18', 'CCL3', 'CCL4', 'CCL3L1', 'CCL25', 'ITCH', 'TFF2')

chem_genes_ileum <- intersect(chem_genes, ileum_gene) #59

ileum <- AddModuleScore_UCell(ileum, features = list(chem_genes), name = "chem_sig")

VlnPlot(ileum, features = "signature_1chem_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") 

plot3 <- VlnPlot(ileum, features = "signature_1chem_sig", group.by = "cell_typev2", split.by = "group", cols = c('#F8766D', '#619CFF'), pt.size = 0, sort = "decreasing") + theme(axis.title.x=element_blank())
png("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure4_12_16_24/VlnPlot_helm_batch1_13_ileum_group_chemokine_recept_module_12_16_24.png", width = 2800, height = 1800, res = 300)
plot(plot3 + ggtitle("Chemokine receptor binding module") +
       ylab("Module score") + theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)))
dev.off()


#------------------Figure 5------------------#

#####Figure 5a#######

#myeloid cell volcano plot of group 1 vs. group 2 donors 

#switch to R v. 4.3

#load in results from dreamlet 
res.dl <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/helm_batch1_13_colon_group_batch_cor_dreamlet_12_10_24.rds")

#extract genes from myeloid cells 

myeloid_cell <- topTable(res.dl[["Myeloid Cell"]], coef = "groupgroup2", number = Inf)

myeloid_cell$group <- 'NS'
myeloid_cell$ID <- rownames(myeloid_cell)
myeloid_cell$group[myeloid_cell$logFC < (0) & myeloid_cell$adj.P.Val<0.05] <-'DOWN'
myeloid_cell$group[myeloid_cell$logFC > (0) & myeloid_cell$adj.P.Val<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
myeloid_cell$label <- ifelse(myeloid_cell$group == "NS", NA, myeloid_cell$group) #might have to do this 

myeloid_cell$label[myeloid_cell$group !='NS'] <- myeloid_cell$ID[myeloid_cell$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=myeloid_cell,mapping=aes(x=logFC,y=-log10(adj.P.Val),color=group,label=label))+
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
  labs(color = myeloid_cell$group) + 
  scale_color_manual(values=c("DOWN"= "#D65F5F", "NS"="#999999", "UP"="#7986CB"), name = "Group", labels = c("DOWN" = "group 1", "UP" = "group 2"))

volcano_gaussian

#save the volcano plot - use the pdf 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure5_12_16_24/VolcanoPlot_helm_batch1_13_colon_group_myeloid_cell_12_16_24.pdf", width = 7, height = 7)
plot(volcano_gaussian)
dev.off()


#####Figure 5b#######

#pathway enrichment plot 

myeloid_cell <- topTable(res.dl[["Myeloid Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #263

myeloid_DEG_down <- myeloid_cell[which(myeloid_cell$logFC < (0)),] #10

myeloid_DEG_up <- myeloid_cell[which(myeloid_cell$logFC > (0)),] #253

myl_g2 <- rownames(myeloid_DEG_up)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = myl_g2, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_myl_g2 <- data.frame(ego) #172 pathways

#visualize results
dotplot(ego, showCategory = 20)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure5_12_16_24/BarPlot_helm_batch1_13_colon_myeloid_cell_pathway_enrichment_12_16_24.pdf", width = 12, height = 18)
barplot(ego, showCategory = 20, font.size = 16)
dev.off()

#####Figure 5c#######

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
  xlab("Mean macrophage activation score") +
  ylab("Mean regulation of inflammatory response score") +
  theme_bw() +
  scale_color_manual(labels = c("group 1", "group 2"), values = c('#D65F5F', '#7986CB'))

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

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure5_12_16_24/ScatterPlot_helm_batch1_13_colon_myeloid_cell_mac_activation_vs_inflam_response_score_3_6_25.pdf", width = 7, height = 7)
plot(plot)
dev.off()


#####Figure 5d#######

#switch back to R version 4.3.1

#load in sce object to extract metadata information 
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/helm_batch1_13_colon_scITD_group_sce_11_27_24.rds")

colon$id <- paste(colon$group, colon$donor, sep = "_")

#create pseudobulk data; specify cluster_id and sample_id
#use counts data stored in assay field
colon$batch <- as.factor(colon$batch)
colon$donor <- as.factor(colon$donor)
colon$ind <- as.factor(colon$ind)
pb <- aggregateToPseudoBulk(colon, assay = "counts", cluster_id = "ctypes", sample_id = "id", verbose = FALSE)

metadata <- pb@colData

metadata <- as.data.frame(metadata)

#PCA of macrophage activation genes colon group myeloid cells 

#load in myeloid genes 
myl <- read.csv(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/colon_pseuodbulk_myeloid_cell.csv")

rownames(myl) <- myl$X


#1/14/2025 - update with relaxed threshold DEGs
mac_gene <- c("SLC11A1", "TLR2", "LRRK2", "C5AR1", "ITGAM", "CCL3", "FCGR3A", "CD93", "PTPRC", "TYROBP", "JAK2")


myl$gene <- rownames(myl)
myl <- myl[,-1]

myl_filt <- myl[myl$gene %in% mac_gene,]

myl_filt <- myl_filt[,-33] #8 genes; update has 11 genes

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


#relaxed DEG threshold
plot <- pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(group))) +
  geom_point(size = 3) +
  theme_bw() +
  guides(color = guide_legend(title = "Group")) +
  theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  xlab("PC1 (61.7%)") +
  ylab("PC2 (11.1%)") +
  scale_color_manual(labels = c("Group 1", "Group 2"), values = c('#D65F5F', '#7986CB'))

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure5_12_16_24/PCA_helm_batch1_13_colon_myeloid_cell_mac_activation_3_6_25.pdf", width = 5.5, height = 5)
plot(plot)
dev.off()


#relaxed DEG threshold
pc_score_metadata$macro_IF <- factor(pc_score_metadata$macro_IF, levels = c("NI", "IF"))
plot <- pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(macro_IF))) +
  geom_point(size = 3) +
  theme_bw() +
  guides(color = guide_legend(title = "Macro IF Status")) +
  theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  xlab("PC1 (61.7%)") +
  ylab("PC2 (11.1%)") +
  scale_color_manual(labels = c("NI", "IF"), values = c('#0E4D92', '#A94064'))
plot

pdf("/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure5_12_16_24/PCA_helm_batch1_13_colon_myeloid_cell_mac_activation_macroIFstatus_3_6_25.pdf", width = 6, height = 5)
plot(plot)
dev.off()

#####Figure 5e#######

#cell chat ADGRE pathway 

cellchat_g1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/colon/data/helm_batch1_13_colon_group1_cellchat_obj_12_5_24.rds")

cellchat_g2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/colon/data/helm_batch1_13_colon_group2_cellchat_obj_12_5_24.rds")

object.list <- list(G1 = cellchat_g1, G2 = cellchat_g2)

#load in cellchat object 
cellchat <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/colon/data/helm_batch1_13_colon_group_comb_cellchat_obj_12_5_24.rds")

#vector of colors to use for figure - load in the colors from figure 2
colon_color <- c("Naive B Cell" = green[1], "Memory B Cell" = green[2], "Cycling B Cell" = green[4], "IgA Plasma Cell" = green[3], "IgG Plasma Cell" = green[5], "CD4 T Cell" = blue[1],"Regulatory T Cell" = blue[5], "NK/CD8 T Cell" = blue[2], "Macrophage" = orange[1], "Monocyte" = orange[2], "Mast Cell" = orange[4], "Stem/TA Cell" = pink[1], "Colonocyte" = pink[2], "Mature Colonocyte" = pink[10], "BEST4 Colonocyte" = pink [4], "Goblet Cell" = pink[3], "Tuft Cell" = pink[5], "Mesenchymal Cell" = purple[1], "Endothelial Cell" = purple[2])

pathways.show <- c("ADGRE") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = c(pathways.show)) # control the edge weights across different datasets
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure5_12_16_24/CirclePlot_batch1_13_colon_group_ADGRE_pathway_12_16_24.pdf", width = 15, height = 8)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], pt.title = 16, color.use = colon_color, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure5_12_16_24/CirclePlot_batch1_13_colon_group1_ADGRE_pathway_2_3_25.pdf", width = 11, height = 12)
netVisual_aggregate(cellchat_g1, pt.title = 16, color.use = colon_color, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1])
dev.off()


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure5_12_16_24/CirclePlot_batch1_13_colon_group2_ADGRE_pathway_2_3_25.pdf", width = 11, height = 12)
netVisual_aggregate(cellchat_g2, pt.title = 16, color.use = colon_color, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1])
dev.off()

#------------------Figure 6------------------#

#####Figure 6a#######

#colonocyte volcano plot 

col <- topTable(res.dl[["Colonocyte"]], coef = "groupgroup1", number = Inf)

col$group <- 'NS'
col$ID <- rownames(col)
col$logFC <- col$logFC * -1
col$group[col$logFC < (0) & col$adj.P.Val<0.05] <-'DOWN'
col$group[col$logFC > (0) & col$adj.P.Val<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
col$label <- ifelse(col$group == "NS", NA, col$group) #might have to do this 

col$label[col$group !='NS'] <- col$ID[col$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=col,mapping=aes(x=logFC,y=-log10(adj.P.Val),color=group,label=label))+
  geom_point()+
  #xlim(-0.228,0.228) +
  #labs(title='volcano plot of significance against DEG recurrence and non-recurrence') +
  xlab("log2FC") +
  geom_label_repel(na.rm = TRUE, show.legend = F, box.padding = 0.5, size = 5) +
  labs() +
  theme_minimal()+ 
  theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16), legend.position = "right") +
  #geom_vline(xintercept=c(-1, 1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  labs(color = col$group) + 
  scale_color_manual(values=c("DOWN"= "#F06292", "NS"="#999999", "UP"="#6A9FB5"), name = "Group", labels = c("DOWN" = "group 1", "UP" = "group 2"))

volcano_gaussian

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure6_12_16_24/VolcanoPlot_helm_batch1_13_rectum_group_colonocyte_12_16_24.pdf", width = 7, height = 7)
plot(volcano_gaussian)
dev.off()


#####Figure 6b#######

#colonocyte pathway enriched in group 1

colonocyte <- topTable(res.dl[["Colonocyte"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #30 DEG
colonocyte_DEG_up <- colonocyte[which(colonocyte$logFC > (0)),] #27
col_g1 <- rownames(colonocyte_DEG_up)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = col_g1, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_col_g1 <- data.frame(ego)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure6_12_16_24/BarPlot_helm_batch1_13_rectum_colonocyte_group1_pathways_12_16_24.pdf", width = 12, height = 16)
barplot(ego, showCategory = 20, font.size = 16)
dev.off()

#####Figure 6c#######

#stem cell volcano plot 

stem <- topTable(res.dl[["Stem Cell"]], coef = "groupgroup1", number = Inf)

stem$group <- 'NS'
stem$ID <- rownames(stem)
stem$logFC <- stem$logFC * -1
stem$group[stem$logFC < (0) & stem$adj.P.Val<0.05] <-'DOWN'
stem$group[stem$logFC > (0) & stem$adj.P.Val<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
stem$label <- ifelse(stem$group == "NS", NA, stem$group) #might have to do this 

stem$label[stem$group !='NS'] <- stem$ID[stem$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=stem,mapping=aes(x=logFC,y=-log10(adj.P.Val),color=group,label=label))+
  geom_point()+
  #xlim(-0.228,0.228) +
  #labs(title='volcano plot of significance against DEG recurrence and non-recurrence') +
  xlab("log2FC") +
  geom_label_repel(na.rm = TRUE, show.legend = F, box.padding = 0.5, size = 5) +
  theme_minimal()+ 
  labs() +
  theme(axis.text.x = element_text(color="black", size=16),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16), legend.position = "right") +
  #geom_vline(xintercept=c(-1, 1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  labs(color = stem$group) + 
  scale_color_manual(values=c("DOWN"= "#F06292", "NS"="#999999", "UP"="#6A9FB5"), name = "Group", labels = c("DOWN" = "group 1", "UP" = "group 2"))

volcano_gaussian

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure6_12_16_24/VolcanoPlot_helm_batch1_13_rectum_group_stem_cell_12_16_24.pdf", width = 7, height = 7)
plot(volcano_gaussian)
dev.off()

#####Figure 6d#######

#stem cell pathway enrichment plot 

stem_cell <- topTable(res.dl[["Stem Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #44 DEG

stem_DEG_up <- stem_cell[which(stem_cell$logFC > (0)),]

#stem group 1 
stem_g1 <- rownames(stem_DEG_up)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = stem_g1, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_stem_g1 <- data.frame(ego)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure6_12_16_24/BarPlot_helm_batch1_13_rectum_stem_cell_group1_pathways_12_16_24.pdf", width = 10, height = 8)
barplot(ego, showCategory = 7, font.size = 16)
dev.off()

#####Figure 6e#######

#switch to R 4.2.1

#lollipop figure of IFN gamma and TNF module score 

#load in rectum with IFN and TNF response 
rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_group_scITD_mod_score_12_11_24.rds")

#re-order donors by group 
rectum$donor <- factor(rectum$donor, levels = c("donor4", "donor10", "donor14", "donor16", "donor18", "donor21", "donor22", "donor25", "donor28", "donor30", "donor31", "donor33", "donor34", "donor1", "donor5", "donor6", "donor7", "donor11", "donor12", "donor13", "donor15", "donor17", "donor19", "donor20", "donor26", "donor27", "donor29", "donor32"))

rectum_metadata <- rectum@meta.data

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

#create

ifn_sig_name <- "IFN Signature"
tnf_sig_name <- "TNF Signature"


test <- ggplot(df_long_cell, aes(x = cell_typev2, y = score, color = group, group = cell_typev2)) +
  geom_line(size = 1.2) +  # Connect IFN and TNF scores with a line
  geom_point(size = 4) +   # Add points for scores
  facet_wrap(~score_type, scales = "free_y", labeller = labeller(mean_ifn_score = ifn_sig_name, mean_tnf_score = tnf_sig_name)) + # Optional: Separate by group
  theme_minimal(base_size = 16) +
  labs(x = "Donor", y = "Score", color = "Group") +
  scale_color_manual(labels = c("group 1", "group 2"), values = c("#F06292", "#6A9FB5")) +
  coord_flip()

test

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure6_12_16_24/Lollipop_helm_batch1_13_rectum_ifn_tnf_mod_score_12_16_24.pdf", width = 12, height = 8)
plot(test)
dev.off()


#------------Figure 7-------------#

#####Figure 7a#######

#GWAS enrichment in ileum data 

#load in marginal results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_conti_2025.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_marg = P) %>% 
  dplyr::arrange(P_marg)

#re order cell types to match UMAP
df_marg$VARIABLE <- factor(df_marg$VARIABLE, levels = c("Naive_B_Cell", "Memory_B_Cell", "Plasma_Cell", "Plasmablast", "CD4_T_Cell", "NK_CD8_T_Cell", "Ambig__T_Cell", "Macrophage", "Monocyte", "Pro_Inflammatory_Monocyte", "Mast_Cell", "Stem_Paneth_Cell", "Enterocyte", "Goblet_Cell", "TA_Cell", "Tuft_Cell", "Mesenchymal_Cell", "Endothelial_Cell"))

cols <- c(green[1], green[2], green[3], green[6], blue[1], blue[2], blue[3], orange[1], orange[2], orange[3], orange[4], pink[1], pink[2], pink[3], pink[6], pink[5], purple[1], purple[2])

bonferroni_threshold <- 0.05 / nrow(df_marg)

plot <- df_marg %>% 
  # create the plot
  ggplot(aes(x = VARIABLE, y = -log10(P_marg), colour = factor(VARIABLE))) +
  geom_point(aes(size = -log10(P_marg))) +
  geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed", color = "red", size = 0.8) +
  theme_classic() +
  scale_size_area(max_size = 9) +
  #guides(color = guide_legend(title = "Cell type")) +
  theme(legend.position = "none", axis.text.x = element_text(color="black", size=16, angle = 60, hjust = 1),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  scale_x_discrete(labels=c("Naive_B_Cell" = "Naive B Cell", "Memory_B_Cell" = "Memory B Cell", "Plasma_Cell" = "Plasma Cell", "Plasmablast" = "Plasmablast", "CD4_T_Cell" = "CD4 T Cell", "NK_CD8_T_Cell" = "NK/CD8 T Cell", "Ambig__T_Cell" = "Ambig. T Cell", "Macrophage" = "Macrophage", "Monocyte" = "Monocyte", "Pro_Inflammatory_Monocyte" = "Pro-Inflam. Monocyte", "Mast_Cell" = "Mast Cell", "Stem_Paneth_Cell" = "Stem/Paneth Cell", "Enterocyte" = "Enterocyte", "Goblet_Cell" = "Goblet Cell", "TA_Cell" = "TA Cell", "Tuft_Cell" = "Tuft Cell", "Mesenchymal_Cell" = "Mesenchymal Cell", "Endothelial_Cell" = "Endothelial Cell")) + 
  xlab("Cell types") +
  ylab("-log10(p value)") +
  scale_color_manual(values = cols)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure7_1_28_25/ScatterPlot_ileum_topcell_bonferroni.pdf", width = 8, height = 7)
plot
dev.off()

#####Figure 7b#######

#GWAS enrichment in colon data 

#read in marginal results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_colon_conti_spe_2025.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_marg = P) %>% 
  dplyr::arrange(P_marg)

df_marg$VARIABLE <- factor(df_marg$VARIABLE, levels = c("Naive_B_Cell", "Memory_B_Cell", "Cycling_B_Cell", "IgA_Plasma_Cell", "IgG_Plasma_Cell", "CD4_T_Cell", "Regulatory_T_Cell", "NK_CD8_T_Cell", "Macrophage", "Monocyte", "Mast_Cell", "Stem_TA_Cell", "Colonocyte", "Mature_Colonocyte", "BEST4_Colonocyte", "Goblet_Cell", "Tuft_Cell", "Mesenchymal_Cell", "Endothelial_Cell"))

cols <- c(green[1], green[2], green[4], green[3], green[5], blue[1], blue[5], blue[2], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink [4], pink[3], pink[5], purple[1], purple[2])

bonferroni_threshold <- 0.05/nrow(df_marg)

plot <- df_marg %>% 
  # create the plot
  ggplot(aes(x = VARIABLE, y = -log10(P_marg), colour = factor(VARIABLE))) +
  geom_point(aes(size = -log10(P_marg))) +
  geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed", color = "red", size = 0.8) +
  theme_classic() +
  scale_size_area(max_size = 9) +
  #guides(color = guide_legend(title = "Cell type")) +
  theme(legend.position = "none", axis.text.x = element_text(color="black", size=16, angle = 60, hjust = 1),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  scale_x_discrete(labels = c("Naive_B_Cell" = "Naive B Cell", "Memory_B_Cell" = "Memory B Cell", "Cycling_B_Cell" = "Cycling B Cell", "IgA_Plasma_Cell" = "IgA Plasma Cell", "IgG_Plasma_Cell" = "IgG Plasma Cell", "CD4_T_Cell" = "CD4 T Cell", "Regulatory_T_Cell" = "Regulatory T Cell", "NK_CD8_T_Cell" = "NK/CD8 T Cell", "Macrophage" = "Macrophage", "Monocyte" = "Monocyte", "Mast_Cell" = "Mast Cell", "Stem_TA_Cell" = "Stem/TA_Cell", "Colonocyte" = "Colonocyte", "Mature_Colonocyte" = "Mature Colonocyte", "BEST4_Colonocyte" = 'BEST4 Colonocyte', "Goblet_Cell" = "Goblet Cell", "Tuft_Cell" = "Tuft Cell", "Mesenchymal_Cell" = "Mesenchymal Cell", "Endothelial_Cell" = "Endothelial Cell")) +
  xlab("Cell types") +
  ylab("-log10(p value)") +
  scale_color_manual(values = cols)


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure7_1_28_25/ScatterPlot_colon_topcell_bonferroni_2_3_25.pdf", width = 8, height = 7)
plot
dev.off()

#####Figure 7c#######

#GWAS enrichment in rectum data 

#load in the marginal results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_rectum_conti_spe_2025.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_marg = P) %>% 
  dplyr::arrange(P_marg)

#re order cell types to match UMAP
df_marg$VARIABLE <- factor(df_marg$VARIABLE, levels = c("Naive_B_Cell", "Memory_B_Cell", "Cycling_B_Cell", "Plasma_Cell", "Ambig__B_Cell", "CD4_T_Cell", "CD8_T_Cell", "NK_CD8_T_Cell", "ILC3", "Ambig__T_Cell", "Macrophage", "Monocyte", "Mast_Cell", "Stem_Cell", "Colonocyte", "Mature_Colonocyte", "Goblet_Cell","TA_Cell", "Tuft_Cell", "Enteroendocrine_Cell", "Mesenchymal_Cell", "Endothelial_Cell"))

cols <- c(green[1], green[2], green[4], green[3], green[7], blue[1], blue[6], blue[2], blue[8], blue[3], orange[1], orange[2], orange[4], pink[1], pink[2], pink[10], pink[3], pink[6], pink[5], pink[12], purple[1], purple[2])


bonferroni_threshold <- 0.05/nrow(df_marg)
plot <- df_marg %>% 
  # create the plot
  ggplot(aes(x = VARIABLE, y = -log10(P_marg), colour = factor(VARIABLE))) +
  geom_point(aes(size = -log10(P_marg))) +
  geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed", color = "red", size = 0.8) +
  theme_classic() +
  scale_size_area(max_size = 9) +
  #guides(color = guide_legend(title = "Cell type")) +
  theme(legend.position = "none", axis.text.x = element_text(color="black", size=16, angle = 60, hjust = 1),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  scale_x_discrete(labels = c("Naive_B_Cell" = "Naive B Cell", "Memory_B_Cell" = "Memory B Cell", "Cycling_B_Cell" = "Cycling B Cell", "Plasma_Cell" = "Plasma Cell", "Ambig__B_Cell" = "Ambig. B Cell", "CD4_T_Cell" = "CD4 T Cell", "CD8_T_Cell" = 'CD8 T Cell', "NK_CD8_T_Cell" = "NK/CD8 T Cell", "ILC3" = "ILC3", "Ambig__T_Cell" = "Ambig. T Cell", "Ambig__T_Cell", "Macrophage" = "Macrophage", "Monocyte" = "Monocyte", "Mast_Cell" = "Mast Cell", "Stem_Cell" = "Stem Cell", "TA_Cell" = "TA Cell", "Colonocyte" = "Colonocyte", "Mature_Colonocyte" = "Mature Colonocyte", "Goblet_Cell" = "Goblet Cell", "Tuft_Cell" = "Tuft Cell", "Enteroendocrine_Cell" = "Enteroendocrine Cell", "Mesenchymal_Cell" = "Mesenchymal Cell", "Endothelial_Cell" = "Endothelial Cell")) +
  xlab("Cell types") +
  ylab("-log10(p value)") +
  scale_color_manual(values = cols)


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/figures/manuscript_figures/figure7_1_28_25/ScatterPlot_rectum_topcell_bonferroni_2_3_25.pdf", width = 8, height = 7)
plot
dev.off()
