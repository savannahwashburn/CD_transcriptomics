#Helmsley GCA Batch 1-13 RECTUM scITD analysis 

#load libraries 
library(scITD)
library(Seurat)
library(dplyr)
library(WGCNA)
library(pheatmap)
library(ComplexHeatmap)

#--------set up the container----------#

#load in seurat object 
#contains donors and ctypes to use for tensor construction 
rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_cell_annotation_12_8_24.rds")

#extract raw counts matrix 
rectum_counts <- rectum@assays$RNA@counts

#create metadata matrix - metadata for each cell 
rectum_metadata <- rectum@meta.data

rectum_metadata <- rectum_metadata[, -c(2:4, 7, 8, 13:15)]

colnames(rectum_metadata) <- c("sample", "batch", "tissue", "Macro_IF", "Micro_IF", "sex", "ancestry", 'donors', 'ctypes')

rectum_metadata$batch <- as.factor(rectum_metadata$batch)
rectum_metadata$donors <- as.factor(rectum_metadata$donors)

#set up project parameters - use same cell types that were previously used; including mature colonocyte separately 
param_list <- initialize_params(ctypes_use = c("B Cell", "Plasma Cell", "NK/T Cell", "Myeloid Cell", "Stem Cell", "TA Cell", "Colonocyte", "Mature Colonocyte", "Goblet Cell"),
                                ncores = 30, rand_seed = 10)

# create project container
rectum_cor <- make_new_container(count_data=rectum_counts, 
                                 meta_data=rectum_metadata,
                                 params=param_list,
                                 label_donor_sex = FALSE)

# form the tensor from the data, include batch effect correction
rectum_cor <- form_tensor(rectum_cor, donor_min_cells=5,
                          norm_method='regular', scale_factor=1000000, batch_var = "batch",
                          vargenes_method='norm_var', vargenes_thresh=500,
                          scale_var = TRUE, var_scale_power = 0.5)

#keeping 28 donors

#check for overdispersion
print(length(rectum_cor[["all_vargenes"]])) #2584

#determine ranks to use 
# get assistance with rank determination
rectum_cor <- determine_ranks_tucker(rectum_cor, max_ranks_test=c(10,15),
                                     shuffle_level='cells', 
                                     num_iter=100, 
                                     batch_var = 'batch',
                                     norm_method='regular',
                                     scale_factor=1000000,
                                     scale_var=TRUE,
                                     var_scale_power=0.5)

rectum_cor$plots$rank_determination_plot

#4 factors 
#10 gene sets 

#run the tucker decomposition
rectum_cor <- run_tucker_ica(rectum_cor, ranks=c(4,10),
                             tucker_type = 'regular', rotation_type = 'hybrid')

# get donor scores-metadata associations
rectum_cor <- get_meta_associations(rectum_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='rsq')

# plot donor scores
rectum_cor <- plot_donor_matrix(rectum_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                                show_donor_ids = TRUE,
                                add_meta_associations='rsq')

# show the donor scores heatmap
rectum_cor$plots$donor_matrix

rectum_cor <- get_meta_associations(rectum_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='pval')


# plot donor scores
rectum_cor <- plot_donor_matrix(rectum_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                                show_donor_ids = TRUE,
                                add_meta_associations='pval')

# show the donor scores heatmap
rectum_cor$plots$donor_matrix

#-------assess stability------------# 
rectum_cor <- run_stability_analysis(rectum_cor, ranks = c(4,10), n_iterations = 100, subset_type = 'subset', sub_prop = .80)
rectum_cor$plots$stability_plot_dsc
rectum_cor$plots$stability_plot_lds

#-------factor loadings------#
#get significant genes
rectum_cor <- get_lm_pvals(rectum_cor)

# generate the loadings plots
rectum_cor <- get_all_lds_factor_plots(rectum_cor, 
                                       use_sig_only=TRUE,
                                       nonsig_to_zero=TRUE,
                                       sig_thresh=.05,
                                       display_genes=FALSE,
                                       gene_callouts = TRUE,
                                       callout_n_gene_per_ctype=3,
                                       show_var_explained = TRUE)

# arrange the plots into a figure and show the figure
rectum_loading_cor <- render_multi_plots(rectum_cor,data_type='loadings')


plot1 <- rectum_loading_cor
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/rectum/figures/Heatmap_factor_loadings_batch1_13_rectum_4fact_10gene_batch_cor_12_8_24.pdf", width = 15, height = 20)
plot(plot1)
dev.off()


#run GSEA on each factor 
rectum_cor <- run_gsea_one_factor(rectum_cor, factor_select=1, method="fgsea", thresh=0.001, db_use=c("GO"), signed=TRUE)
rectum_cor <- run_gsea_one_factor(rectum_cor, factor_select=2, method="fgsea", thresh=0.001, db_use=c("GO"), signed=TRUE)
rectum_cor <- run_gsea_one_factor(rectum_cor, factor_select=3, method="fgsea", thresh=0.05, db_use=c("GO"), signed=TRUE)
rectum_cor <- run_gsea_one_factor(rectum_cor, factor_select=4, method="fgsea", thresh=0.01, db_use=c("GO"), signed=TRUE)


plot <- rectum_cor[["plots"]][["gsea"]][["1"]]
plot <- draw(plot, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 50), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/rectum/figures/factor1_batch_cor_gsea_12_8_24.pdf", width = 18, height = 20)
plot(plot)
dev.off()

plot2 <- rectum_cor[["plots"]][["gsea"]][["2"]]
plot2 <- draw(plot2, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 120), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/rectum/figures/factor2_batch_cor_gsea_12_8_24.pdf", width = 12, height = 16)
plot(plot2)
dev.off()

plot3 <- rectum_cor[["plots"]][["gsea"]][["3"]]
plot3 <- draw(plot3, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 120), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/rectum/figures/factor3_batch_cor_gsea_12_8_24.pdf", width = 14, height = 10)
plot(plot3)
dev.off()

plot4 <- rectum_cor[["plots"]][["gsea"]][["4"]]
plot4 <- draw(plot4, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 120), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/rectum/figures/factor4_batch_cor_gsea_12_8_24.pdf", width = 14, height = 20)
plot(plot4)
dev.off()

f1_pvals <- get_one_factor_gene_pvals(rectum_cor, factor_select = 1)

#number of significant genes 
sum(f1_pvals < 0.05) #2246

f2_pvals <- get_one_factor_gene_pvals(rectum_cor, factor_select = 2)

#number of significant genes 
sum(f2_pvals < 0.05) #2699

f3_pvals <- get_one_factor_gene_pvals(rectum_cor, factor_select = 3)

#number of significant genes 
sum(f3_pvals < 0.05) #1886

f4_pvals <- get_one_factor_gene_pvals(rectum_cor, factor_select = 4)

#number of significant genes 
sum(f4_pvals < 0.05) #1159


#------extract factor 1 scores------#

f1_data <- get_one_factor(rectum_cor, factor_select = 1)
f1_scores <- f1_data[[1]]
f1_loadings <- f1_data[[2]]

f1_scores <- as.data.frame(f1_scores)

#top 5% quantile values 
top_5_percentile_value <- quantile(f1_scores$V1, 0.95)
top_5_percent_df <- f1_scores[f1_scores$V1 >= top_5_percentile_value, ]
#donor 4 and 30 

#bottom 5% quantile values 
bottom_5_percentile_value <- quantile(f1_scores$V1, 0.05)

# Filter the dataframe to get the bottom 5% values
bottom_5_percent_df <- f1_scores[f1_scores$V1 <= bottom_5_percentile_value, ]
#donor 27 and 4

#sort the donor scores 
# sort by mpg
f1_score_sort <- f1_scores[order(f1_scores$V1), , drop = F]

#save the colon container 

saveRDS(rectum_cor, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/colon/data/helm_batch1_13_rectum_batch_cor_container_12_8_24.rds")
