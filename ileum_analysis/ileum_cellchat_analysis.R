#Helmsley GCA Batch 1-13 ILEUM Cell Chat

#11/22/2024

#load libraries
library(CellChat)
library(patchwork)
library(dplyr)
library(NMF)
library(ggalluvial)


#---------Assess group 1 data-----------#

#load in seurat object - 34411 cells
ileum_g1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_group1_11_22_24.rds")

Idents(ileum_g1) <- ileum_g1$cell_typev2

ileum_g1$donor <- factor(ileum_g1$donor, levels = c("donor1", "donor2", "donor3", "donor4", "donor5", "donor8", "donor9", "donor11", "donor14", "donor18", "donor20", "donor23", "donor27"))

#prepare data for cellchat obj
data.input <- ileum_g1[["RNA"]]@data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(ileum_g1)

Idents(ileum_g1) <- ileum_g1$donor

samples <- Idents(ileum_g1)
meta <- data.frame(labels = labels, samples = samples, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_g1 <- createCellChat(object = data.input, meta = meta, group.by = "labels")

#set the ligand-receptor database 
CellChatDB <- CellChatDB.human 

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# set the used database in the object
cellchat_g1@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat_g1 <- subsetData(cellchat_g1) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat_g1 <- identifyOverExpressedGenes(cellchat_g1)
cellchat_g1 <- identifyOverExpressedInteractions(cellchat_g1)

#The number of highly variable ligand-receptor pairs used for signaling inference is 1900

#compute communication probability and infer cellular communication network
cellchat_g1 <- computeCommunProb(cellchat_g1, type = "triMean")

#filter
cellchat_g1 <- filterCommunication(cellchat_g1, min.cells = 10)

#infer cell-cell communication at signaling pathway level
cellchat_g1 <- computeCommunProbPathway(cellchat_g1)

#calculate aggregated cell-cell communication network 
cellchat_g1 <- aggregateNet(cellchat_g1)

groupSize <- as.numeric(table(cellchat_g1@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_g1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_g1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#examine signaling within each cell group 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/figures/helm_batch1_13_ileum_group1_agg_celltype_11_22_24.pdf", width = 10, height = 20)
mat <- cellchat_g1@net$weight
#par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  print(netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]))
}
dev.off()

#pathways in group 1
g1_pathways <- cellchat_g1@netP$pathways

#[1] "MIF"       "MHC-I"     "COLLAGEN"  "APP"       "MHC-II"    "CD45"     
#[7] "LAMININ"   "GALECTIN"  "CypA"      "MK"        "ApoA"      "CLEC"     
#[13] "ANNEXIN"   "ADGRG"     "CDH"       "CD99"      "ApoB"      "PECAM1"   
#[19] "JAM"       "CTSG"      "PECAM2"    "CEACAM"    "ADGRE"     "ICAM"     
#[25] "FN1"       "DESMOSOME" "BAFF"      "SEMA4"     "GRN"       "SELL"     
#[31] "IL1"       "CXCL"      "CLDN"      "PARs"      "CD160"     "CDH1"     
#[37] "IL16"      "SIRP"      "CCL"       "VEGF"      "GUCA"      "CD23"     
#[43] "LCK"       "SEMA7"     "NRG"       "CD96"      "CD40"      "TGFb"     
#[49] "CALCR"     "HSPG"      "CSF"       "TENASCIN"  "CDH5"      "LIFR"     
#[55] "CD86"      "PLAU"      "PVR"       "ESAM"      "LAIR1"     "EPHB"     
#[61] "PDGF"      "PTPRM"     "NOTCH"     "CADM" 

#visualize each pathway 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/figures/helm_batch1_13_ileum_group1_pathways_chord_11_23_24.pdf", width = 10, height = 10)
for (i in 1:length(g1_pathways)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  print(netVisual_aggregate(cellchat_g1, signaling = g1_pathways[i], layout = "chord"))
  #ggsave(filename=paste0(ni_pathways[i], "circle_plot_7_8_24.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
dev.off()


# Compute the network centrality scores
cellchat_g1 <- netAnalysis_computeCentrality(cellchat_g1, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network signaling pathways

netAnalysis_signalingRole_scatter(cellchat_g1)

#save cell chat object group 1
saveRDS(cellchat_g1, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/data/helm_batch1_13_ileum_group1_cellchat_obj_11_23_24.rds")

#-------------cell chat group 2------------#

#11/23/2024

#load in seurat object - 34411 cells
ileum_g2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_group2_11_22_24.rds")

Idents(ileum_g2) <- ileum_g2$cell_typev2

ileum_g2$donor <- factor(ileum_g2$donor, levels = c("donor6", "donor7", "donor10", "donor12", "donor17", "donor21", "donor24", "donor26", "donor30", "donor34"))

#prepare data for cellchat obj
data.input <- ileum_g2[["RNA"]]@data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(ileum_g2)

Idents(ileum_g2) <- ileum_g2$donor

samples <- Idents(ileum_g2)

meta <- data.frame(labels = labels, samples = samples, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_g2 <- createCellChat(object = data.input, meta = meta, group.by = "labels")

#set the ligand-receptor database 
CellChatDB <- CellChatDB.human 

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# set the used database in the object
cellchat_g2@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat_g2 <- subsetData(cellchat_g2) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat_g2 <- identifyOverExpressedGenes(cellchat_g2)
cellchat_g2 <- identifyOverExpressedInteractions(cellchat_g2)

#The number of highly variable ligand-receptor pairs used for signaling inference is 1896

#compute communication probability and infer cellular communication network
cellchat_g2 <- computeCommunProb(cellchat_g2, type = "triMean")

#filter
cellchat_g2 <- filterCommunication(cellchat_g2, min.cells = 10)

#infer cell-cell communication at signaling pathway level
cellchat_g2 <- computeCommunProbPathway(cellchat_g2)

#calculate aggregated cell-cell communication network 
cellchat_g2 <- aggregateNet(cellchat_g2)

groupSize <- as.numeric(table(cellchat_g2@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_g2@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_g2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#examine signaling within each cell group 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/figures/helm_batch1_13_ileum_group2_agg_celltype_11_23_24.pdf", width = 10, height = 20)
mat <- cellchat_g2@net$weight
#par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  print(netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]))
}
dev.off()

#pathways in group 2
g2_pathways <- cellchat_g2@netP$pathways

#"MIF"       "MHC-I"     "APP"       "COLLAGEN"  "MHC-II"    "CD45"      "CypA"     
#[8] "MK"        "GALECTIN"  "LAMININ"   "CLEC"      "ANNEXIN"   "BAFF"      "CD99"     
#[15] "ADGRE"     "PECAM2"    "CEACAM"    "PECAM1"    "CCL"       "CTSG"      "JAM"      
#[22] "ADGRG"     "ICAM"      "FN1"       "CXCL"      "CDH"       "SEMA4"     "DESMOSOME"
#[29] "CLDN"      "ApoA"      "ApoB"      "THBS"      "GRN"       "PARs"      "SIRP"     
#[36] "CD23"      "CDH1"      "LCK"       "IL1"       "IL16"      "SELL"      "CD160"    
#[43] "APRIL"     "CD40"      "VCAM"      "CHEMERIN"  "VEGF"      "CD96"      "CD86"     
#[50] "NECTIN"    "CDH5"      "LAIR1"     "LIFR"      "CADM"      "TGFb"      "CD34"     
#[57] "PVR"       "KLK"       "PDGF"   

overlap <- intersect(g1_pathways, g2_pathways) #52 pathways 

g1_only <- setdiff(g1_pathways, g2_pathways) #12 pathways 

#"GUCA"     "SEMA7"    "NRG"      "CALCR"    "HSPG"     "CSF"      "TENASCIN"
#[8] "PLAU"     "ESAM"     "EPHB"     "PTPRM"    "NOTCH"   

g2_only <- setdiff(g2_pathways, g1_pathways) #7 pathways 

#"THBS"     "APRIL"    "VCAM"     "CHEMERIN" "NECTIN"   "CD34"     "KLK"  

#visualize each pathway 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/figures/helm_batch1_13_ileum_group2_pathways_chord_11_23_24.pdf", width = 10, height = 10)
for (i in 1:length(g2_pathways)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  print(netVisual_aggregate(cellchat_g2, signaling = g2_pathways[i], layout = "chord"))
  #ggsave(filename=paste0(ni_pathways[i], "circle_plot_7_8_24.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
dev.off()


# Compute the network centrality scores
cellchat_g2 <- netAnalysis_computeCentrality(cellchat_g2, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network signaling pathways

netAnalysis_signalingRole_scatter(cellchat_g2)

#save cell chat object group 1
saveRDS(cellchat_g2, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/data/helm_batch1_13_ileum_group2_cellchat_obj_11_23_24.rds")


#--------------joint analysis-------------#

#11/23/2024

#load in cell chat g1 and cell chat g2

cellchat_g1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/data/helm_batch1_13_ileum_group1_cellchat_obj_11_23_24.rds")

cellchat_g2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/data/helm_batch1_13_ileum_group2_cellchat_obj_11_23_24.rds")

object.list <- list(G1 = cellchat_g1, G2 = cellchat_g2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#An object of class CellChat created from a merged object with multiple datasets 
#1104 signaling genes.
#67645 cells. 
#CellChat analysis of single cell RNA-seq data!

#compare interactions and strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#G1 vs G2 (red is increased in G2 blue is increased in g1)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weight - ", names(object.list)[i]))
}

#pathway figure
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)

#DEA 

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "G1"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "G1",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "G2",ligand.logFC = -0.05, receptor.logFC = NULL)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat) #111
gene.down <- extractGeneSubsetFromPair(net.down, cellchat) #88

#CDH (cadherin) signaling difference 

pathways.show <- c("CDH") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#ApoA

pathways.show <- c("ApoA") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#ApoB

pathways.show <- c("ApoB") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#icam 

pathways.show <- c("ICAM") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#BAFF

pathways.show <- c("BAFF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#save the merged cellchat object 
saveRDS(cellchat, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/ileum/data/helm_batch1_13_ileum_group_cellchat_obj_merge_11_25_24.rds")


