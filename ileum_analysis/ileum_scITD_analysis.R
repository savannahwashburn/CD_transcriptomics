#11/19/2024 Helmsley GCA Batch 1-13 ILEUM scITD analysis 

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
ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_11_19_24.rds")

#extract raw counts matrix 
ileum_counts <- ileum@assays$RNA@counts

#create metadata matrix - metadata for each cell 
ileum_metadata <- ileum@meta.data

ileum_metadata <- ileum_metadata[, -c(2:4, 7, 8, 13:15)]

colnames(ileum_metadata) <- c("sample", "batch", "tissue", "Macro_IF", "Micro_IF", "sex", "ancestry", 'donors', 'ctypes')

ileum_metadata$batch <- as.factor(ileum_metadata$batch)
ileum_metadata$donors <- as.factor(ileum_metadata$donors)

#set up project parameters 
param_list <- initialize_params(ctypes_use = c("B Cell", "Plasma Cell", "NK/T Cell", "Myeloid Cell", "Stem/Paneth Cell", "Enterocyte", "Goblet Cell", "TA Cell"),
                                ncores = 30, rand_seed = 10)

# create project container
ileum_cor <- make_new_container(count_data=ileum_counts, 
                                meta_data=ileum_metadata,
                                params=param_list,
                                label_donor_sex = FALSE)

# form the tensor from the data, include batch effect correction
ileum_cor <- form_tensor(ileum_cor, donor_min_cells=5,
                         norm_method='regular', scale_factor=1000000, batch_var = "batch",
                         vargenes_method='norm_var', vargenes_thresh=500,
                         scale_var = TRUE, var_scale_power = 0.5)

#check for overdispersion
print(length(ileum_cor[["all_vargenes"]])) #2513

#determine ranks to use 
# get assistance with rank determination
ileum_cor <- determine_ranks_tucker(ileum_cor, max_ranks_test=c(10,15),
                                    shuffle_level='cells', 
                                    num_iter=100, 
                                    batch_var = 'batch',
                                    norm_method='regular',
                                    scale_factor=1000000,
                                    scale_var=TRUE,
                                    var_scale_power=0.5)

ileum_cor$plots$rank_determination_plot

#2 factors 
#9 gene sets 

#run the tucker decomposition
ileum_cor <- run_tucker_ica(ileum_cor, ranks=c(2,9),
                            tucker_type = 'regular', rotation_type = 'hybrid')
# get donor scores-metadata associations
ileum_cor <- get_meta_associations(ileum_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='rsq')

# plot donor scores
ileum_cor <- plot_donor_matrix(ileum_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                               show_donor_ids = TRUE,
                               add_meta_associations='rsq')

# show the donor scores heatmap
ileum_cor$plots$donor_matrix

ileum_cor <- get_meta_associations(ileum_cor, vars_test=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'), stat_use='pval')


# plot donor scores
ileum_cor <- plot_donor_matrix(ileum_cor, meta_vars=c('batch', 'Macro_IF', 'Micro_IF', 'sex', 'ancestry'),
                               show_donor_ids = TRUE,
                               add_meta_associations='pval')

# show the donor scores heatmap
ileum_cor$plots$donor_matrix

#-------assess stability------------# 
ileum_cor <- run_stability_analysis(ileum_cor, ranks = c(2,9), n_iterations = 100, subset_type = 'subset', sub_prop = .80)
ileum_cor$plots$stability_plot_dsc
ileum_cor$plots$stability_plot_lds

#-------factor loadings------#
#get significant genes
ileum_cor <- get_lm_pvals(ileum_cor)

# generate the loadings plots
ileum_cor <- get_all_lds_factor_plots(ileum_cor, 
                                      use_sig_only=TRUE,
                                      nonsig_to_zero=TRUE,
                                      sig_thresh=.05,
                                      display_genes=FALSE,
                                      gene_callouts = TRUE,
                                      callout_n_gene_per_ctype=3,
                                      show_var_explained = TRUE)

# arrange the plots into a figure and show the figure
ileum_loading_cor <- render_multi_plots(ileum_cor,data_type='loadings')
ileum_loading_cor

#run GSEA on each factor 
ileum_cor <- run_gsea_one_factor(ileum_cor, factor_select=1, method="fgsea", thresh=0.05, db_use=c("GO"), signed=TRUE)
ileum_cor <- run_gsea_one_factor(ileum_cor, factor_select=2, method="fgsea", thresh=0.001, db_use=c("GO"), signed=TRUE)


plot <- ileum_cor[["plots"]][["gsea"]][["1"]]
plot <- draw(plot, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 50), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/ileum/figures/factor1_batch_cor_gsea_11_19_24.pdf", width = 9, height = 11)
plot(plot)
dev.off()

plot2 <- ileum_cor[["plots"]][["gsea"]][["2"]]
plot2 <- draw(plot2, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 120), "mm"))
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/ileum/figures/factor2_batch_cor_gsea_11_19_24.pdf", width = 14, height = 20)
plot(plot2)
dev.off()

f1_pvals <- get_one_factor_gene_pvals(ileum_cor, factor_select = 1)

#number of significant genes 
sum(f1_pvals < 0.05) #4770

f2_pvals <- get_one_factor_gene_pvals(ileum_cor, factor_select = 2)

#number of significant genes 
sum(f2_pvals < 0.05) #1893

#------extract factor 2 scores------#

f2_data <- get_one_factor(ileum_cor, factor_select = 2)
f2_scores <- f2_data[[1]]
f2_loadings <- f2_data[[2]]

#top 5% quantile values 
top_5_percentile_value <- quantile(f2_scores$V1, 0.95)
top_5_percent_df <- f2_scores[f2_scores$V1 >= top_5_percentile_value, ]
#donor 4 and 23 

#bottom 5% quantile values 
bottom_5_percentile_value <- quantile(f2_scores$V1, 0.05)

# Filter the dataframe to get the bottom 5% values
bottom_5_percent_df <- f2_scores[f2_scores$V1 <= bottom_5_percentile_value, ]
#donor 10 and 24

#sort the donor scores 
# sort by mpg
f2_score_sort <- f2_scores[order(f2_scores$V1), , drop = F]

#save the ileum container 

saveRDS(ileum_cor, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/scITD/ileum/helm_batch1_13_ileum_batch_cor_container_11_20_24.rds")

