#Helmsley GCA Batch 1-13 Rectum Dreamlet Analysis 

#load in libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(DESeq2)
library(ggplot2)
library(sctransform)
library(SeuratDisk)
library(SeuratData)
library(RColorBrewer)
library(ArchR)
library(pheatmap)
library(Dune)
library(dreamlet)
library(muscat)
library(ExperimentHub)
library(zenith)
library(scater)
library(SingleCellExperiment)

#---------account for batch---------#

#--------load in single cell experiment object--------#


rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/helm_batch1_13_rectum_scITD_info_sce_12_8_24.rds")


#create unique identifier for each sample
rectum$id <- paste(rectum$group, rectum$donor, sep = "_")

#create pseudobulk data; specify cluster_id and sample_id
#use counts data stored in assay field
rectum$batch <- as.factor(rectum$batch)
rectum$donor <- as.factor(rectum$donor)

pb <- aggregateToPseudoBulk(rectum, assay = "counts", cluster_id = "ctypes", sample_id = "id", verbose = FALSE)
assayNames(pb)

#normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, ~batch + group, min.count = 5)

details(res.proc)

plotVoom(res.proc)


#variance partition
vp.lst <- fitVarPart(res.proc, ~batch + group + macro_IF + micro_IF + sex + ancestry)
genes <- vp.lst$gene[1:5]
plotPercentBars(vp.lst[vp.lst$gene %in% genes, ])
plotVarPart(vp.lst, label.angle = 60)

#-------perform differential expression----------#
# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~batch + group)

# names of estimated coefficients
coefNames(res.dl)

# the resulting object of class dreamletResult
# stores results and other information
res.dl

#volcano plot 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/figures/VolcanoPlot_batch1_13_rectum_batch_cor_group_12_9_24.pdf", width = 10, height = 15)
plotVolcano(res.dl, coef = "groupgroup1")
dev.off()

#extract results 
b_cell <- topTable(res.dl[["B Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #1

write.csv(b_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_b_cell_batch_cor_group1_12_9_24.csv")

plasma_cell <- topTable(res.dl[["Plasma Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #2 DEGs 

write.csv(plasma_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_plasma_cell_batch_cor_group1_12_9_24.csv")

nkt_cell <- topTable(res.dl[["NK/T Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #0 DEGs

myeloid_cell <- topTable(res.dl[["Myeloid Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #73


myeloid_DEG_down <- myeloid_cell[which(myeloid_cell$logFC < (0)),] #7

write.csv(myeloid_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_myeloid_cell_batch_cor_group2_relax_12_10_24.csv")

myeloid_DEG_up <- myeloid_cell[which(myeloid_cell$logFC > (0)),] #66

write.csv(myeloid_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_myeloid_cell_batch_cor_group1_relax_12_10_24.csv")

stem_cell <- topTable(res.dl[["Stem Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #44 DEG

stem_DEG_down <- stem_cell[which(stem_cell$logFC < (0)),] #2

write.csv(stem_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_stem_cell_batch_cor_group2_relax_12_10_24.csv")

stem_DEG_up <- stem_cell[which(stem_cell$logFC > (0)),] #42

write.csv(stem_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_stemcell_batch_cor_group1_relax_12_10_24.csv")

ta_cell <- topTable(res.dl[["TA Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #6 DEG

write.csv(ta_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_ta_cell_batch_cor_group1_12_9_24.csv")

colonocyte <- topTable(res.dl[["Colonocyte"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #30 DEG


colonocyte_DEG_down <- colonocyte[which(colonocyte$logFC < (0)),] #3

write.csv(colonocyte_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_colonocyte_batch_cor_group2_relax_12_10_24.csv")


colonocyte_DEG_up <- colonocyte[which(colonocyte$logFC > (0)),] #27

write.csv(colonocyte_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_colonocyte_batch_cor_group1_relax_12_10_24.csv")

mat_col <- topTable(res.dl[["Mature Colonocyte"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #49 DEGs


mat_col_DEG_down <- mat_col[which(mat_col$logFC < (0)),] #8

#write.csv(mat_col_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_mat_colonocyte_batch_cor_group2_12_9_24.csv")

mat_col_DEG_up <- mat_col[which(mat_col$logFC > (0)),] #41

write.csv(mat_col_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_mat_colonocyte_batch_cor_group1_relax_12_10_24.csv")

goblet_cell <- topTable(res.dl[["Goblet Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #1

write.csv(goblet_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_gob_batch_cor_group1_12_9_24.csv")

tuft_cell <- topTable(res.dl[["Tuft Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #11

tuft_DEG_down <- tuft_cell[which(tuft_cell$logFC < (0)),] #1

write.csv(tuft_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_tuft_cell_batch_cor_group2_relax_12_10_24.csv")

entero <- topTable(res.dl[["Enteroendocrine Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #none

mes_cell <- topTable(res.dl[["Mesenchymal Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #1

write.csv(mes_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/data/DEG_mesenchymal_cell_batch_cor_group1_12_9_24.csv")

endo_cell <- topTable(res.dl[["Endothelial Cell"]], coef = "groupgroup1", p.value = 0.05, number = Inf) #none

#check direction 
# get data
df <- extractData(res.proc, "Colonocyte")


# make plot
ggplot(df, aes(group, GBP1)) +
  geom_boxplot() +
  #thm +
  ylab(bquote(Expression ~ (log[2] ~ CPM))) +
  ggtitle("GBP1")

#positive log2FC is up in group 1
#negative log2FC is up in group 2

#-------pathway analysis---------#


library(DOSE)
library(pathview)
library(clusterProfiler)

#colonocyte group 1

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

#visualize results
dotplot(ego, showCategory = 27)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/figures/pathway_col_group1_relax_12_10_24.pdf", width = 10, height = 14)
barplot(ego, showCategory = 20)
dev.off()


#mature colonocyte group 1
mat_col_g1 <- rownames(mat_col_DEG_up)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = mat_col_g1, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_mat_col_g1 <- data.frame(ego)

#visualize results
dotplot(ego, showCategory = 2)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/figures/pathway_mat_col_group1_relax_12_10_24.pdf", width = 8, height = 10)
barplot(ego, showCategory = 2)
dev.off()

#myeloid cell group 1
myeloid_g1 <- rownames(myeloid_DEG_up)

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

#visualize results
dotplot(ego, showCategory = 4)


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/figures/pathway_myl_col_group1_relax_12_10_24.pdf", width = 8, height = 10)
barplot(ego, showCategory = 8)
dev.off()

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

#visualize results
dotplot(ego, showCategory = 7)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/rectum/figures/pathway_stem_group1_relax_12_10_24.pdf", width = 8, height = 10)
barplot(ego, showCategory = 7)
dev.off()

#-------volcano plots of DEG results---------#

#-----flipped sign so red is group 1 and blue is group 2-------------#

#colonocyte
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
  geom_label_repel(na.rm = TRUE, show.legend = T, box.padding = 0.5) +
  theme_minimal()+ 
  theme(legend.position="none")+
  #geom_vline(xintercept=c(-1, 1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  labs(color = col$group) + 
  scale_color_manual(values=c("DOWN"= "#800000", "NS"="#999999", "UP"="#57A0D2"))

volcano_gaussian


#mature colonocyte
mat_col$group <- 'NS'
mat_col$ID <- rownames(mat_col)
mat_col$logFC <- mat_col$logFC * -1
mat_col$group[mat_col$logFC < (0) & mat_col$adj.P.Val<0.05] <-'DOWN'
mat_col$group[mat_col$logFC > (0) & mat_col$adj.P.Val<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
mat_col$label <- ifelse(mat_col$group == "NS", NA, mat_col$group) #might have to do this 

mat_col$label[mat_col$group !='NS'] <- mat_col$ID[mat_col$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=mat_col,mapping=aes(x=logFC,y=-log10(adj.P.Val),color=group,label=label))+
  geom_point()+
  #xlim(-0.228,0.228) +
  #labs(title='volcano plot of significance against DEG recurrence and non-recurrence') +
  xlab("log2FC") +
  geom_label_repel(na.rm = TRUE, show.legend = T, box.padding = 0.5) +
  theme_minimal()+ 
  theme(legend.position="none")+
  #geom_vline(xintercept=c(-1, 1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  labs(color = mat_col$group) + 
  scale_color_manual(values=c("DOWN"= "#800000", "NS"="#999999", "UP"="#57A0D2"))

volcano_gaussian

#myeloid cells
myeloid_cell$group <- 'NS'
myeloid_cell$ID <- rownames(myeloid_cell)
myeloid_cell$logFC <- myeloid_cell$logFC * -1
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
  geom_label_repel(na.rm = TRUE, show.legend = T, box.padding = 0.5) +
  theme_minimal()+ 
  theme(legend.position="none")+
  #geom_vline(xintercept=c(-1, 1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  labs(color = myeloid_cell$group) + 
  scale_color_manual(values=c("DOWN"= "#800000", "NS"="#999999", "UP"="#57A0D2"))

volcano_gaussian

#stem cell 
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
  geom_label_repel(na.rm = TRUE, show.legend = T, box.padding = 0.5) +
  theme_minimal()+ 
  theme(legend.position="none")+
  #geom_vline(xintercept=c(-1, 1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  labs(color = stem$group) + 
  scale_color_manual(values=c("DOWN"= "#800000", "NS"="#999999", "UP"="#57A0D2"))

volcano_gaussian
