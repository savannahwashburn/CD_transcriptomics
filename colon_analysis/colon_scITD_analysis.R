#Helmsley GCA Batch 1-13 COLON scITD analysis 

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
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_cell_annotation_11_26_24.rds")

#extract raw counts matrix 
colon_counts <- colon@assays$RNA@counts

#create metadata matrix - metadata for each cell 
colon_metadata <- colon@meta.data

colon_metadata <- colon_metadata[, -c(2:4, 7, 8, 13:14, 16)]

colnames(colon_metadata) <- c("sample", "batch", "tissue", "Macro_IF", "Micro_IF", "sex", "ancestry", 'donors', 'ctypes')

colon_metadata$batch <- as.factor(colon_metadata$batch)
colon_metadata$donors <- as.factor(colon_metadata$donors)

#set up project parameters - use same cell types that were previously used; including mature colonocyte separately 
param_list <- initialize_params(ctypes_use = c("B Cell", "Plasma Cell", "NK/T Cell", "Myeloid Cell", "Stem/TA Cell", "Colonocyte", "Mature Colonocyte", "Goblet Cell", "Mesenchymal Cell"),
                                ncores = 30, rand_seed = 10)

# create project container
colon_cor <- make_new_container(count_data=colon_counts, 
                                meta_data=colon_metadata,
                                params=param_list,
                                label_donor_sex = FALSE)

# form the tensor from the data, include batch effect correction
colon_cor <- form_tensor(colon_cor, donor_min_cells=5,
                         norm_method='regular', scale_factor=1000000, batch_var = "batch",
                         vargenes_method='norm_var', vargenes_thresh=500,
                         scale_var = TRUE, var_scale_power = 0.5)

#keeping 32 donors

#check for overdispersion
print(length(colon_cor[["all_vargenes"]])) #2547

#determine ranks to use 
# get assistance with rank determination
colon_cor <- determine_ranks_tucker(colon_cor, max_ranks_test=c(10,15),
                                    shuffle_level='cells', 
                                    num_iter=100, 
                                    batch_var = 'batch',
                                    norm_method='regular',
                                    scale_factor=1000000,
                                    scale_var=TRUE,
                                    var_scale_power=0.5)

colon_cor$plots$rank_determination_plot

#4 factors 
#10 gene sets 

#run the tucker decomposition
colon_cor <- run_tucker_ica(colon_cor, ranks=c(4,10),
                            tucker_type = 'regular', rotation_type = 'hybrid')
# get donor scores-metadata associations
colon_cor <- get_meta_associations(colon_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='rsq')

# plot donor scores
colon_cor <- plot_donor_matrix(colon_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                               show_donor_ids = TRUE,
                               add_meta_associations='rsq')

# show the donor scores heatmap
colon_cor$plots$donor_matrix

colon_cor <- get_meta_associations(colon_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='pval')


# plot donor scores
colon_cor <- plot_donor_matrix(colon_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                               show_donor_ids = TRUE,
                               add_meta_associations='pval')

# show the donor scores heatmap
colon_cor$plots$donor_matrix

#-------assess stability------------# 
colon_cor <- run_stability_analysis(colon_cor, ranks = c(4,10), n_iterations = 100, subset_type = 'subset', sub_prop = .80)
colon_cor$plots$stability_plot_dsc
colon_cor$plots$stability_plot_lds

#-------factor loadings------#
#get significant genes
colon_cor <- get_lm_pvals(colon_cor)

# generate the loadings plots
colon_cor <- get_all_lds_factor_plots(colon_cor, 
                                      use_sig_only=TRUE,
                                      nonsig_to_zero=TRUE,
                                      sig_thresh=.05,
                                      display_genes=FALSE,
                                      gene_callouts = TRUE,
                                      callout_n_gene_per_ctype=3,
                                      show_var_explained = TRUE)

# arrange the plots into a figure and show the figure
colon_loading_cor <- render_multi_plots(colon_cor,data_type='loadings')


plot1 <- colon_loading_cor
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/colon/figures/Heatmap_factor_loadings_batch1_13_colon_4fact_10gene_11_26_24.pdf", width = 15, height = 20)
plot(plot1)
dev.off()


#run GSEA on each factor 
colon_cor <- run_gsea_one_factor(colon_cor, factor_select=1, method="fgsea", thresh=0.05, db_use=c("GO"), signed=TRUE)
colon_cor <- run_gsea_one_factor(colon_cor, factor_select=2, method="fgsea", thresh=0.001, db_use=c("GO"), signed=TRUE)
colon_cor <- run_gsea_one_factor(colon_cor, factor_select=3, method="fgsea", thresh=0.05, db_use=c("GO"), signed=TRUE)
colon_cor <- run_gsea_one_factor(colon_cor, factor_select=4, method="fgsea", thresh=0.001, db_use=c("GO"), signed=TRUE)


plot <- colon_cor[["plots"]][["gsea"]][["1"]]
plot <- draw(plot, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 50), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/colon/figures/factor1_batch_cor_gsea_11_26_24.pdf", width = 9, height = 11)
plot(plot)
dev.off()

plot2 <- colon_cor[["plots"]][["gsea"]][["2"]]
plot2 <- draw(plot2, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 120), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/colon/figures/factor2_batch_cor_gsea_11_26_24.pdf", width = 14, height = 20)
plot(plot2)
dev.off()

plot3 <- colon_cor[["plots"]][["gsea"]][["3"]]
plot3 <- draw(plot3, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 120), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/colon/figures/factor3_batch_cor_gsea_11_26_24.pdf", width = 14, height = 10)
plot(plot3)
dev.off()

plot4 <- colon_cor[["plots"]][["gsea"]][["4"]]
plot4 <- draw(plot4, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 120), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/colon/figures/factor4_batch_cor_gsea_11_26_24.pdf", width = 14, height = 20)
plot(plot4)
dev.off()

f1_pvals <- get_one_factor_gene_pvals(colon_cor, factor_select = 1)

#number of significant genes 
sum(f1_pvals < 0.05) #1880

f2_pvals <- get_one_factor_gene_pvals(colon_cor, factor_select = 2)

#number of significant genes 
sum(f2_pvals < 0.05) #2788

f3_pvals <- get_one_factor_gene_pvals(colon_cor, factor_select = 3)

#number of significant genes 
sum(f3_pvals < 0.05) #2617

f4_pvals <- get_one_factor_gene_pvals(colon_cor, factor_select = 4)

#number of significant genes 
sum(f4_pvals < 0.05) #1965


#------extract factor 2 scores------#

f2_data <- get_one_factor(colon_cor, factor_select = 2)
f2_scores <- f2_data[[1]]
f2_loadings <- f2_data[[2]]

f2_scores <- as.data.frame(f2_scores)

#top 5% quantile values 
top_5_percentile_value <- quantile(f2_scores$V1, 0.95)
top_5_percent_df <- f2_scores[f2_scores$V1 >= top_5_percentile_value, ]
#donor 4 and 27 

#bottom 5% quantile values 
bottom_5_percentile_value <- quantile(f2_scores$V1, 0.05)

# Filter the dataframe to get the bottom 5% values
bottom_5_percent_df <- f2_scores[f2_scores$V1 <= bottom_5_percentile_value, ]
#donor 16_s1 and 30

#sort the donor scores 
# sort by mpg
f2_score_sort <- f2_scores[order(f2_scores$V1), , drop = F]

#save the colon container 

saveRDS(colon_cor, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/colon/data/helm_batch1_13_colon_batch_cor_container_11_27_24.rds")
