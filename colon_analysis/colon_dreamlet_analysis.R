#Helmsley GCA Batch 1-13 Colon Dreamlet Analysis 

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


colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/helm_batch1_13_colon_scITD_group_sce_11_27_24.rds")

#include indvidual as random effect 

#create unique identifier for each sample
colon$id <- paste(colon$group, colon$donor, sep = "_")

#create pseudobulk data; specify cluster_id and sample_id
#use counts data stored in assay field
colon$batch <- as.factor(colon$batch)
colon$donor <- as.factor(colon$donor)
colon$ind <- as.factor(colon$ind)
pb <- aggregateToPseudoBulk(colon, assay = "counts", cluster_id = "ctypes", sample_id = "id", verbose = FALSE)
assayNames(pb)

#normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, ~batch + group + (1|ind), min.count = 5)

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
res.dl <- dreamlet(res.proc, ~batch + group + (1|ind))

# names of estimated coefficients
coefNames(res.dl)

# the resulting object of class dreamletResult
# stores results and other information
res.dl

#volcano plot 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/figures/VolcanoPlot_batch1_13_colon_batch_cor_group_re_ind_11_27_24.pdf", width = 10, height = 15)
plotVolcano(res.dl, coef = "groupgroup2")
dev.off()

#extract results 
b_cell <- topTable(res.dl[["B Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #1

write.csv(b_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_b_cell_batch_cor_re_ind_group2_11_27_24.csv")

plasma_cell <- topTable(res.dl[["Plasma Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #1 DEGs 

write.csv(plasma_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_plasma_cell_batch_cor_re_ind_group2_11_27_24.csv")

nkt_cell <- topTable(res.dl[["NK/T Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #2 DEGs

write.csv(nkt_cell, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_nkt_cell_batch_cor_re_ind_group2_11_27_24.csv")

myeloid_cell <- topTable(res.dl[["Myeloid Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #263

myeloid_DEG_down <- myeloid_cell[which(myeloid_cell$logFC < (0)),] #10
write.csv(myeloid_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_myeloid_cell_batch_cor_re_ind_group1_all_12_10_24.csv")

myeloid_DEG_up <- myeloid_cell[which(myeloid_cell$logFC > (0)),] #253
write.csv(myeloid_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_myeloid_cell_batch_cor_re_ind_group2_all_12_10_24.csv")

stempan_cell <- topTable(res.dl[["Stem/TA Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #3 DEG

stempan_DEG_down <- stempan_cell[which(stempan_cell$logFC < (0)),] #1
write.csv(stempan_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_stem_ta_cell_batch_cor_re_ind_group1_all_12_10_24.csv")

stempan_DEG_up <- stempan_cell[which(stempan_cell$logFC >= (1)),] #2

write.csv(stempan_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_stem_ta_cell_batch_cor_re_ind_group2_11_27_24.csv")

colonocyte <- topTable(res.dl[["Colonocyte"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #4 DEG

colonocyte_DEG_down <- colonocyte[which(colonocyte$logFC <= (-1)),] #2

write.csv(colonocyte_DEG_down, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_colonocyte_batch_cor_re_ind_group1_11_27_24.csv")

colonocyte_DEG_up <- colonocyte[which(colonocyte$logFC >= (1)),] #2

write.csv(colonocyte_DEG_up, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_colonocyte_batch_cor_re_ind_group2_11_27_24.csv")

mat_col <- topTable(res.dl[["Mature Colonocyte"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #2 DEGs

write.csv(mat_col, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/DEG_mat_colonocyte_batch_cor_re_ind_group2_11_27_24.csv")


goblet_cell <- topTable(res.dl[["Goblet Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #0


tuft_cell <- topTable(res.dl[["Tuft Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #none

mes_cell <- topTable(res.dl[["Mesenchymal Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #none

endo_cell <- topTable(res.dl[["Endothelial Cell"]], coef = "groupgroup2", p.value = 0.05, number = Inf) #none


#check direction 

# get data
df <- extractData(res.proc, "B Cell")


# make plot
ggplot(df, aes(group, SPINK4)) +
  geom_boxplot() +
  #thm +
  ylab(bquote(Expression ~ (log[2] ~ CPM))) +
  ggtitle("SPINK4")

#group 2 is positive log2FC
#group 1 is negative log2FC

#save the dreamlet object 
saveRDS(res.dl, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/data/helm_batch1_13_colon_group_batch_cor_dreamlet_12_10_24.rds")

#----------pathway analysis--------#
library(DOSE)
library(pathview)
library(clusterProfiler)

#pathway analysis and volcano plot 
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

pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/pseudobulk_deg/colon/figures/pathway_myeloid_group2_relax_12_10_24.pdf", width = 10, height = 14)
barplot(ego, showCategory = 20)
dev.off()

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
  geom_label_repel(na.rm = TRUE, show.legend = T, box.padding = 0.5) +
  theme_minimal()+ 
  theme(legend.position="none")+
  #geom_vline(xintercept=c(-1, 1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  labs(color = myeloid_cell$group) + 
  scale_color_manual(values=c("DOWN"= "#800000", "NS"="#999999", "UP"="#57A0D2"))

volcano_gaussian


#---------save gene expression results-----------#
#load in dreamlet from group comparison 
# Save each element to a separate CSV file
for (i in seq_along(res.proc)) {
  write.csv(res.proc[[i]][["E"]], paste0("colon_pseuodbulk_", i, ".csv"), row.names = TRUE)
}


#overlap between macrophage activation and myeloid cell differentiation pathway 
go_genes <- go_summary_myl_g2$geneID
myl_dif_gene <- go_genes[6]

myl_dif_gene <- strsplit(myl_dif_gene, "/")

myl_dif_gene <- unlist(myl_dif_gene, use.names=FALSE)

mac_gene <- go_genes[91]

mac_gene <- strsplit(mac_gene, "/")

mac_gene <- unlist(mac_gene, use.names=FALSE)

intersect(myl_dif_gene, mac_gene)

mac_gene <- c("SLC11A1", "TLR2", "LRRK2", "C5AR1", "ITGAM", "CCL3", "FCGR3A", "CD93", "PTPRC", "TYROBP", "JAK2")

