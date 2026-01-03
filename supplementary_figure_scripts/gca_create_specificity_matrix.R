#script to perform merged analysis - MAGMA for all 3 sites 


#create new label names - to facilitate combined analysis
#some cell types are dataset specific

ileum$MAGMA_label <- ileum$cell_typev2
colon$MAGMA_label <- colon$cell_typev1
rectum$MAGMA_label <- rectum$cell_typev2

#merge the data
gca <- merge(ileum, y = c(colon, rectum))

#set RNA assay 
DefaultAssay(gca) <- "RNA"

#calculate average expression 
gca_counts <- AverageExpression(gca, group.by = "MAGMA_label", assay = "RNA")

#save 
write.csv(gca_counts, file = "/storage/home/swashburn30/Helmsley/Batch1_13/MAGMA_analysis/data/helm_batch1_13_gca_avg_gene_exp.csv")

#------------create specificity matrix----------------#

gene_coordinates <- 
  read_tsv("/storage/home/swashburn30/Helmsley/Batch1_13/MAGMA_analysis/NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  select(1:4) %>% 
  rename(end="X4", start="X3", chr="X2", ENTREZ="X1") %>% 
  mutate(chr=paste0("chr",chr))

entrez_ensembl <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL)

# Get Entrez to Gene Symbol mappings
entrez_symbol <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL)

entrez_ensembl_unique_genes_entrez <- entrez_ensembl %>% count(gene_id) %>% filter(n==1)
entrez_ensembl_unique_genes_ens <- entrez_ensembl %>% count(ensembl_id) %>% filter(n==1)
entrez_ensembl <- filter(entrez_ensembl,gene_id%in%entrez_ensembl_unique_genes_entrez$gene_id & ensembl_id %in% entrez_ensembl_unique_genes_ens$ensembl_id)
colnames(entrez_ensembl)[1] <- "ENTREZ"
colnames(entrez_ensembl)[2] <- "Gene"

# Merge with Gene Symbols
entrez_ensembl <- left_join(entrez_ensembl, entrez_symbol, by = c("ENTREZ" = "gene_id"))

# Rename Gene Symbol column
colnames(entrez_ensembl)[3] <- "Gene_Symbol"

gene_coordinates <- inner_join(entrez_ensembl,gene_coordinates) %>% as.tibble()

gca_count <- read.csv(file = "/storage/home/swashburn30/Helmsley/Batch1_13/MAGMA_analysis/data/helm_batch1_13_gca_avg_gene_exp.csv")

gca_count <- gca_count %>% add_count(X) %>% 
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-X) %>%
  as.tibble()

exp_agg <- gca_count %>% rename(ClusterID=column, Expr_sum_mean=Expr)

not_expressed <- exp_agg %>% 
  group_by(X) %>% 
  summarise(total_sum=sum(Expr_sum_mean)) %>% 
  filter(total_sum==0) %>% 
  select(X) %>% unique() 

exp_agg <- filter(exp_agg,!X%in%not_expressed$X)

exp_agg <- exp_agg %>% 
  group_by(ClusterID) %>% 
  mutate(Expr_sum_mean=Expr_sum_mean*1000/sum(Expr_sum_mean)) %>% 
  ungroup()

exp_agg <- exp_agg %>% 
  group_by(X) %>% 
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>% 
  ungroup()

colnames(exp_agg) <- c("Gene_Symbol", "ClusterID", "Expr_sum_mean", "specificity")
exp_agg <- inner_join(exp_agg,gene_coordinates,by="Gene_Symbol")

exp_conti_spe <- exp_agg %>% select(ENTREZ, ClusterID, specificity) %>% spread(ClusterID, specificity)
colnames(exp_conti_spe)[1] <- "GENE"

#save specificity matrix 
write_tsv(exp_conti_spe, file = "/storage/home/swashburn30/Helmsley/Batch1_13/MAGMA_analysis/data/gca_conti_specificity_matrix.txt")

#----------filter for top cell types------------#

#!/usr/bin/env Rscript

#script to create top results for conditional analysis 

#load libraries 
library(readr)
library(dplyr)


#---------------gca results--------------#

#read in specificity matrix 
df <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/gca_conti_specificity_matrix.txt")

#read in step 2 results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_gca_conti_spe_2025.gsa.out", skip=4) %>%
  arrange(P) %>%
  filter(P < 0.05/n()) #keep significant clusters after bonferroni correction

sig_list <- df_marg$VARIABLE

df_sig <- df[c("GENE", sig_list)] #regulatory t cell, monocyte, and pro-inflammatory monocyte

#save the output 
write_tsv(df_sig, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_gca_conti_matrix_sig.txt")

#---------------conditional analysis-----------------#

#read in significant results 
df_cond <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_gca_joint_spe_sig.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_cond = P)

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

df_marg %>% 
  # create the plot
  ggplot(aes(x = VARIABLE, y = -log10(P_marg), colour = factor(VARIABLE))) +
  geom_point(aes(size = -log10(P_marg))) +
  geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed", color = "red", size = 0.8) +
  theme_classic() +
  scale_size_area(max_size = 9) +
  #guides(color = guide_legend(title = "Cell type")) +
  theme(legend.position = "none", axis.text.x = element_text(color="black", size=16, angle = 60, hjust = 1),axis.text.y = element_text(color="black", size=16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) + 
  scale_x_discrete(labels = c("RNA.Naive.B.Cell" = "Naive B Cell", "RNA.Memory.B.Cell" = "Memory B Cell", "Cycling_B_Cell" = "Cycling B Cell", "Plasma_Cell" = "Plasma Cell", "Ambig__B_Cell" = "Ambig. B Cell", "CD4_T_Cell" = "CD4 T Cell", "CD8_T_Cell" = 'CD8 T Cell', "NK_CD8_T_Cell" = "NK/CD8 T Cell", "ILC3" = "ILC3", "Ambig__T_Cell" = "Ambig. T Cell", "Ambig__T_Cell", "Macrophage" = "Macrophage", "Monocyte" = "Monocyte", "Mast_Cell" = "Mast Cell", "Stem_Cell" = "Stem Cell", "TA_Cell" = "TA Cell", "Colonocyte" = "Colonocyte", "Mature_Colonocyte" = "Mature Colonocyte", "Goblet_Cell" = "Goblet Cell", "Tuft_Cell" = "Tuft Cell", "Enteroendocrine_Cell" = "Enteroendocrine Cell", "Mesenchymal_Cell" = "Mesenchymal Cell", "Endothelial_Cell" = "Endothelial Cell")) +
  xlab("Cell types") +
  ylab("-log10(p value)")
