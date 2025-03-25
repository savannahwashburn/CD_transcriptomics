#create specificity matrix 


#calculate average gene expression - use Seurat AverageGene

#------Ileum Low Resolution----------#
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

ileum_count <- read.csv(file = "/storage/home/swashburn30/Helmsley/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_avg_gene_exp.csv")

ileum_count_test <- ileum_count %>% add_count(X) %>% 
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-X) %>%
  as.tibble()

exp_agg <- ileum_count_test %>% rename(ClusterID=column, Expr_sum_mean=Expr)

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
write_tsv(exp_conti_spe, file = "ileum_low_res_conti_sepcificity_matrix.txt")

#-------ileum high resolution--------#

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

ileum_count <- read.csv(file = "/storage/home/swashburn30/Helmsley/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_avg_gene_exp_high_res.csv")

ileum_count_test <- ileum_count %>% add_count(X) %>% 
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-X) %>%
  as.tibble()

exp_agg <- ileum_count_test %>% rename(ClusterID=column, Expr_sum_mean=Expr)

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
write_tsv(exp_conti_spe, file = "ileum_high_res_conti_sepcificity_matrix.txt")


#-------ileum referenece-------#
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

ileum_count <- read.csv(file = "/storage/home/swashburn30/Helmsley/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_avg_gene_exp_ref.csv")

ileum_count_test <- ileum_count %>% add_count(X) %>% 
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-X) %>%
  as.tibble()

exp_agg <- ileum_count_test %>% rename(ClusterID=column, Expr_sum_mean=Expr)

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
write_tsv(exp_conti_spe, file = "ileum_ref_conti_sepcificity_matrix.txt")


