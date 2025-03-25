#!/usr/bin/env Rscript

#script to perform forward selection condition - MAGMA

#load libraries 
library(readr)
library(dplyr)
library(tidyr)

#read in significant results 
df_cond <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_joint_spe_sig.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_cond = P)

#read in marginal results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_conti_2025.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_marg = P) %>% 
  dplyr::arrange(P_marg)

n_sig <- ceiling(sqrt(nrow(df_cond)))
n_clusters <- nrow(df_marg)

#combine marginal and conditional results 
df_comb <- left_join(df_cond, df_marg, by = "VARIABLE")
df_comb <- df_comb %>% 
  mutate(VarCode = rep(c("a", "b"),times=nrow(df_comb)/2)) %>%
  pivot_wider(names_from = VarCode, values_from=c(VARIABLE, P_cond, P_marg)) %>%
  mutate(PS_a=log10(P_cond_a)/log10(P_marg_a), PS_b=log10(P_cond_b)/log10(P_marg_b)) %>%
  stopifnot(df_comb$P_marg_a <= df_comb$P_marg_b)
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
#pro inflammatory monocyte and CD4 t cell


#-----------high resolution stably assigned cell conditional analysis-----------#

#read in significant results 
df_cond <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_high_res_joint_spe_sig.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_cond = P)

#read in marginal results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_high_res_conti_spe_2025.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_marg = P) %>% 
  dplyr::arrange(P_marg)

n_sig <- ceiling(sqrt(nrow(df_cond)))
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
#inflammatory monocyte and regulatory T cell 

#----------Reference------------#

#read in significant results 
df_cond <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_ref_joint_spe_sig.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_cond = P)

#read in marginal results 
df_marg <- read_table(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/MAGMA_analysis/data/helm_batch1_13_ileum_ref_conti_spe_2025.gsa.out", skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_marg = P) %>% 
  dplyr::arrange(P_marg)

n_sig <- ceiling(sqrt(nrow(df_cond)))
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
#inflammatory monocyte and regulatory T cell 
