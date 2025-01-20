#11/20/2024 Helmsley GCA Batch 1-13 Ileum Dreamlet Analysis 

#switch to R v 4.3.1

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

ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/helm_batch1_13_ileum_scITD_group_sce_11_20_24.rds")


#create unique identifier for each sample
ileum$id <- paste(ileum$group, ileum$donor, sep = "_")

#create pseudobulk data; specify cluster_id and sample_id
#use counts data stored in assay field
ileum$batch <- as.factor(ileum$batch)
ileum$donor <- as.factor(ileum$donor)
pb <- aggregateToPseudoBulk(ileum, assay = "counts", cluster_id = "ctypes", sample_id = "id", verbose = FALSE)
assayNames(pb)


#----------test DEG without correcting for batch (which is confounded by individual)---------#

#normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, ~group, min.count = 5)

details(res.proc)

plotVoom(res.proc)

#-------perform differential expression----------#
# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group)

# names of estimated coefficients
coefNames(res.dl)

# the resulting object of class dreamletResult
# stores results and other information
res.dl

#volcano plot 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/figures/VolcanoPlot_group_11_20_24.pdf", width = 10, height = 15)
plotVolcano(res.dl, coef = "groupgroup2")
dev.off()

#extract results 
b_cell <- topTable(res.dl[["B Cell"]], coef = "groupgroup2", p.value = 0.10, number = Inf) #110 DEGs

plasma_cell <- topTable(res.dl[["Plasma Cell"]], coef = "groupgroup2", number = Inf) #none

nkt_cell <- topTable(res.dl[["NK/T Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #63 DEGs

ambigt_cell <- topTable(res.dl[["Ambig. T Cell"]], coef = "groupgroup2", number = Inf) #none

myeloid_cell <- topTable(res.dl[["Myeloid Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #3 DEGs

stempan_cell <- topTable(res.dl[["Stem/Paneth Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #5 DEG

ent <- topTable(res.dl[["Enterocyte"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #286 DEG

goblet_cell <- topTable(res.dl[["Goblet Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #185 DEGs

ta_cell <- topTable(res.dl[["TA Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #13 DEGs

tuft_cell <- topTable(res.dl[["Tuft Cell"]], coef = "groupgroup2", number = Inf) #none

mes_cell <- topTable(res.dl[["Mesenchymal Cell"]], coef = "groupgroup2", number = Inf) #none

endo_cell <- topTable(res.dl[["Endothelial Cell"]], coef = "groupgroup2", number = Inf) #none

#--save the DEG output---#

#b cell 
write.csv(b_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_b_cell.csv")

b_DEG_up <- b_cell[which(b_cell$logFC >= (1)),]#25

b_DEG_up$gene <- rownames(b_DEG_up)

#save
write.csv(b_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_b_cell_group2.csv")

b_DEG_down <- b_cell[which(b_cell$logFC <= (-1)),] #6

#save
write.csv(b_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_b_cell_group1.csv")

#nkt cell 
write.csv(nkt_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_nk_t_cell.csv")

nkt_DEG_up <- nkt_cell[which(nkt_cell$logFC >= (1)),]#6

#save
write.csv(nkt_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_nkt_cell_group2.csv")

nkt_DEG_down <- nkt_cell[which(nkt_cell$logFC <= (-1)),] #24

#save
write.csv(nkt_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_nkt_cell_group1.csv")


#myeloid cell 
write.csv(myeloid_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_myeloid_cell.csv")

myl_DEG_up <- myeloid_cell[which(myeloid_cell$logFC >= (1)),]#3

#save
write.csv(myl_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_myeloid_cell_group2.csv")

myl_DEG_down <- myeloid_cell[which(myeloid_cell$logFC <= (-1)),] #0

#save
#write.csv(myl_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_nkt_cell_group1.csv")

#stem/paneth cell 
write.csv(stempan_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_stem_paneth_cell.csv")

stempan_DEG_up <- stempan_cell[which(stempan_cell$logFC >= (1)),]#3

#save
write.csv(stempan_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_stem_pan_cell_group2.csv")

stempan_DEG_down <- stempan_cell[which(stempan_cell$logFC <= (-1)),] #1

#save
write.csv(stempan_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_stem_pan_cell_group1.csv")


#enterocyte
write.csv(ent, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_enterocyte.csv")

ent_DEG_up <- ent[which(ent$logFC >= (1)),]#68

#save enterocyte with relaxed thresholds 1/2/2024
ent_DEG_up <- ent[which(ent$logFC > (0)),] #143
write.csv(ent_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_enterocyte_group2_1_2_25.csv")

ent_DEG_down <- ent[which(ent$logFC < (0)),] #143
write.csv(ent_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_enterocyte_group1_1_2_25.csv")

#save
write.csv(ent_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_ent_cell_group1.csv")



#save
write.csv(ent_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_ent_cell_group2.csv")

ent_DEG_down <- ent[which(ent$logFC <= (-1)),] #63

#save
write.csv(ent_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_ent_cell_group1.csv")

#goblet cell 
write.csv(goblet_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_goblet_cell.csv")

gob_DEG_up <- goblet_cell[which(goblet_cell$logFC >= (1)),]#50

#save
write.csv(gob_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_goblet_cell_group2.csv")

gob_DEG_down <- goblet_cell[which(goblet_cell$logFC <= (-1)),] #23

#save
write.csv(gob_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_goblet_cell_group1.csv")

#ta cell 
write.csv(ta_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_ta_cell.csv")

ta_DEG_up <- ta_cell[which(ta_cell$logFC >= (1)),]#13

#save
write.csv(ta_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/data/DEG_ta_cell_group2.csv")



#check direction 

# get data
df <- extractData(res.proc, "NK/T Cell")


# make plot
ggplot(df, aes(group, CUBN)) +
  geom_boxplot() +
  #thm +
  ylab(bquote(Expression ~ (log[2] ~ CPM))) +
  ggtitle("CUBN")

#group 2 is positive log2FC
#group 1 is negative log2FC

#------use the non-batch corrected data---------#

library(DOSE)
library(pathview)
library(clusterProfiler)

#enterocytes group 2

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
go_summary_b_g2 <- data.frame(ego)

#visualize results
dotplot(ego, showCategory = 16)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/figures/pathway_enterocyte_group2_11_20_24.pdf", width = 8, height = 10)
barplot(ego, showCategory = 16)
dev.off()

#enterocytes group 1
ent_g1 <- rownames(ent_DEG_down)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = ent_g1, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "MF", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_ent_g1 <- data.frame(ego)

#visualize results
dotplot(ego, showCategory = 16)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/figures/pathway_enterocyte_group1_11_20_24.pdf", width = 8, height = 10)
barplot(ego, showCategory = 16)
dev.off()

ent_test <- topTable(res.dl[["Enterocyte"]], coef = "groupgroup2", number = Inf)
ent_DEG_up <- ent$logFC > 1 & ent$adj.P.Val<0.05#77
ent_DEG_up <- ent_DEG_up[which(ent_DEG_up$adj.P.Val <= 0.05),]


ent_DEG_down <- filtered_b_cell_tests[which(filtered_b_cell_tests$logFC <= (-1)),]#77

ent_test$group <- 'NS'
ent_test$ID <- rownames(ent_test)
ent_test$group[ent_test$logFC < (0) & ent_test$adj.P.Val<0.05] <-'DOWN'
ent_test$group[ent_test$logFC > (0) & ent_test$adj.P.Val<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
ent_test$label <- ifelse(ent_test$group == "NS", NA, ent_test$group) #might have to do this 

ent_test$label[ent_test$group !='NS'] <- ent_test$ID[ent_test$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=ent_test,mapping=aes(x=logFC,y=-log10(adj.P.Val),color=group,label=label))+
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
  labs(color = ent_test$group) + 
  scale_color_manual(values=c("DOWN"= "#800000", "NS"="#999999", "UP"="#57A0D2"))

volcano_gaussian

gob_g2 <- rownames(gob_DEG_up)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = gob_g2, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_gob_g2 <- data.frame(ego)

#visualize results
dotplot(ego, showCategory = 16)

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/ileum/figures/pathway_gob_group2_11_20_24.pdf", width = 8, height = 10)
barplot(ego, showCategory = 16)
dev.off()

gob_test <- topTable(res.dl[["Goblet Cell"]], coef = "groupgroup2", number = Inf)

gob_test$group <- 'NS'
gob_test$ID <- rownames(gob_test)
gob_test$group[gob_test$logFC < (0) & gob_test$adj.P.Val<0.05] <-'DOWN'
gob_test$group[gob_test$logFC > (0) & gob_test$adj.P.Val<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
gob_test$label <- ifelse(gob_test$group == "NS", NA, gob_test$group) #might have to do this 

gob_test$label[gob_test$group !='NS'] <- gob_test$ID[gob_test$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=gob_test,mapping=aes(x=logFC,y=-log10(adj.P.Val),color=group,label=label))+
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
  labs(color = gob_test$group) + 
  scale_color_manual(values=c("DOWN"= "#800000", "NS"="#999999", "UP"="#57A0D2"))

volcano_gaussian

#12/9/2024 NK/T cells

#NK/T cells 

nkt_g2 <- rownames(nkt_DEG_up)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = nkt_g2, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

nkt_g1 <- rownames(nkt_DEG_down)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = nkt_g1, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_gob_g2 <- data.frame(ego) #no significant pathways


nkt_cell <- topTable(res.dl[["NK/T Cell"]], coef = "groupgroup2", number = Inf)

nkt_cell$group <- 'NS'
nkt_cell$ID <- rownames(nkt_cell)
nkt_cell$group[nkt_cell$logFC < (0) & nkt_cell$adj.P.Val<0.05] <-'DOWN'
nkt_cell$group[nkt_cell$logFC > (0) & nkt_cell$adj.P.Val<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
nkt_cell$label <- ifelse(nkt_cell$group == "NS", NA, nkt_cell$group) #might have to do this 

nkt_cell$label[nkt_cell$group !='NS'] <- nkt_cell$ID[nkt_cell$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=nkt_cell,mapping=aes(x=logFC,y=-log10(adj.P.Val),color=group,label=label))+
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
  labs(color = nkt_cell$group) + 
  scale_color_manual(values=c("DOWN"= "#800000", "NS"="#999999", "UP"="#57A0D2"))

volcano_gaussian

