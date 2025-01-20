#Helmsley GCA Batch 1-13 RECTUM Cell Chat

#load libraries
library(CellChat)
library(patchwork)
library(dplyr)
library(NMF)
library(ggalluvial)


#---------Assess group 1 data-----------#

#load in seurat object - 49344 cells
rectum_g1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_g1_12_11_24.rds")

Idents(rectum_g1) <- rectum_g1$cell_typev2

rectum_g1$donor <- factor(rectum_g1$donor, levels = c("donor4", "donor10", "donor14", "donor16", "donor18", "donor21", "donor22", "donor25", "donor28", "donor30", "donor31", "donor33", "donor34"))

#prepare data for cellchat obj
data.input <- rectum_g1[["RNA"]]@data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(rectum_g1)

Idents(rectum_g1) <- rectum_g1$donor

samples <- Idents(rectum_g1)
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

#The number of highly variable ligand-receptor pairs used for signaling inference is 1976

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
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/figures/helm_batch1_13_rectum_group1_agg_celltype_12_11_24.pdf", width = 10, height = 20)
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


#visualize each pathway 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/figures/helm_batch1_13_rectum_group1_pathways_chord_12_11_24.pdf", width = 10, height = 10)
for (i in 1:length(g1_pathways)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  print(netVisual_aggregate(cellchat_g1, signaling = g1_pathways[i], layout = "chord"))
  #ggsave(filename=paste0(ni_pathways[i], "circle_plot_7_8_24.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
dev.off()

#[1] "MIF"       "MHC-I"     "APP"       "COLLAGEN"  "CypA"      "CEACAM"   
#[7] "CD45"      "CLEC"      "MK"        "LAMININ"   "CLDN"      "GALECTIN" 
#[13] "MHC-II"    "DESMOSOME" "CD99"      "CDH"       "ANNEXIN"   "ADGRE"    
#[19] "BAFF"      "THBS"      "PECAM1"    "ADGRG"     "JAM"       "PTN"      
#[25] "CCL"       "PARs"      "FN1"       "SIRP"      "SEMA4"     "SELL"     
#[31] "CTSG"      "IGFBP"     "VISFATIN"  "GUCA"      "PECAM2"    "CDH1"     
#[37] "GRN"       "NRG"       "CD23"      "CD160"     "CD96"      "SEMA7"    
#[43] "LCK"       "VEGF"      "SELPLG"    "NECTIN"    "IL16"      "CD40"     
#[49] "EPHB"      "BTLA"      "OX40"      "CXCL"      "KLK"       "ICAM"     
#[55] "LIFR"      "HSPG"      "CDH5"      "ncWNT"     "MPZ"       "PVR"      
#[61] "CALCR"     "SEMA3"     "PLAU"      "ESAM"      "MADCAM"    "CD34"     
#[67] "BAG"       "EPHA"      "CADM"      "PTPRM"     "GAP"       "OCLN"  

# Compute the network centrality scores
cellchat_g1 <- netAnalysis_computeCentrality(cellchat_g1, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network signaling pathways

netAnalysis_signalingRole_scatter(cellchat_g1)

netVisual_chord_cell(cellchat_g1, signaling = "BTLA", title.name = paste0("BTLA signaling network"))

netVisual_chord_cell(cellchat_g1, signaling = "ncWNT", title.name = paste0("ncWNT signaling network"))

netVisual_chord_cell(cellchat_g1, signaling = "CADM", title.name = paste0("CADM signaling network"))

#save cell chat object group 1
saveRDS(cellchat_g1, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/data/helm_batch1_13_rectum_group1_cellchat_obj_12_11_24.rds")

#------------Assess group 2 data------------#

#load in seurat object - 44914 cells
rectum_g2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_scITD_info_g2_12_11_24.rds")

Idents(rectum_g2) <- rectum_g2$cell_typev2

rectum_g2$donor <- factor(rectum_g2$donor, levels = c("donor1", "donor5", "donor6", "donor7", "donor11", "donor12", 'donor13', "donor15", "donor17", "donor19", "donor20", "donor26", "donor27", "donor29", "donor32"))

#prepare data for cellchat obj
data.input <- rectum_g2[["RNA"]]@data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(rectum_g2)

Idents(rectum_g2) <- rectum_g2$donor

samples <- Idents(rectum_g2)

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

#The number of highly variable ligand-receptor pairs used for signaling inference is 1985

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
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/figures/helm_batch1_13_rectum_group2_agg_celltype_12_11_24.pdf", width = 10, height = 20)
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

#[1] "MHC-I"     "MIF"       "APP"       "COLLAGEN"  "CypA"      "CD45"      "LAMININ"  
#[8] "CLEC"      "CEACAM"    "CLDN"      "CD99"      "DESMOSOME" "CDH"       "MK"       
#[15] "MHC-II"    "GALECTIN"  "FN1"       "ANNEXIN"   "PTN"       "JAM"       "GUCA"     
#[22] "CCL"       "CXCL"      "PECAM1"    "ADGRG"     "CDH1"      "BAFF"      "PARs"     
#[29] "CTSG"      "SIRP"      "THBS"      "ADGRE"     "SELL"      "VISFATIN"  "IGFBP"    
#[36] "SEMA4"     "PECAM2"    "LCK"       "GRN"       "CD46"      "CD160"     "IL1"      
#[43] "CDH5"      "LIFR"      "CD96"      "LT"        "CD23"      "EPHB"      "CD40"     
#[50] "VEGF"      "LIGHT"     "IL16"      "ICAM"      "EGF"       "TNF"       "TENASCIN" 
#[57] "OCLN"      "OX40"      "NECTIN"    "MADCAM"    "PVR"       "SEMA3"     "PLAU"     
#[64] "NOTCH"     "KLK"       "ESAM"      "GAP"       "PDGF"      "HSPG"      "EDN"      
#[71] "GIPR"      "EPHA"      "MPZ"       "WNT"       "PTPRM"     "NEGR"  

overlap <- intersect(g1_pathways, g2_pathways) #63 pathways 

g1_only <- setdiff(g1_pathways, g2_pathways) #6 pathways 

#"NRG"    "SEMA7"  "SELPLG" "BTLA"   "ncWNT"  "CALCR"  "CD34"   "BAG"    "CADM"  

g2_only <- setdiff(g2_pathways, g1_pathways) #15 pathways 

#[1] "CD46"     "IL1"      "LT"       "LIGHT"    "EGF"      "TNF"      "TENASCIN" "NOTCH"   
#[9] "PDGF"     "EDN"      "GIPR"     "WNT"      "NEGR"     

#visualize each pathway 
pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/figures/helm_batch1_13_rectum_group2_pathways_chord_12_11_24.pdf", width = 10, height = 10)
for (i in 1:length(g2_pathways)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  print(netVisual_aggregate(cellchat_g2, signaling = g2_pathways[i], layout = "chord"))
  #ggsave(filename=paste0(ni_pathways[i], "circle_plot_7_8_24.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
dev.off()

netVisual_chord_cell(cellchat_g2, signaling = "TNF", title.name = paste0("TNF signaling network"))

netVisual_chord_cell(cellchat_g2, signaling = "LT", title.name = paste0("LT signaling network"))

netVisual_chord_cell(cellchat_g2, signaling = "LIGHT", title.name = paste0("LIGHT signaling network"))

netVisual_chord_cell(cellchat_g2, signaling = "APP", title.name = paste0("APP signaling network"))

netVisual_chord_cell(cellchat_g2, signaling = "CD160", title.name = paste0("C160 signaling network"))


# Compute the network centrality scores
cellchat_g2 <- netAnalysis_computeCentrality(cellchat_g2, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network signaling pathways

netAnalysis_signalingRole_scatter(cellchat_g2)

#save cell chat object group 2
saveRDS(cellchat_g2, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/data/helm_batch1_13_rectum_group2_cellchat_obj_12_11_24.rds")



#-------group comparison-------#


#load in cell chat g1 and cell chat g2

cellchat_g1 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/data/helm_batch1_13_rectum_group1_cellchat_obj_12_11_24.rds")

cellchat_g2 <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/data/helm_batch1_13_rectum_group2_cellchat_obj_12_11_24.rds")

object.list <- list(G1 = cellchat_g1, G2 = cellchat_g2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat

#An object of class CellChat created from a merged object with multiple datasets 
#1104 signaling genes.
#72760 cells. 
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


gene.up <- extractGeneSubsetFromPair(net.up, cellchat) #112
gene.down <- extractGeneSubsetFromPair(net.down, cellchat) #108


#OX40 signaling difference 

pathways.show <- c("OX40") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#MHC-CII signaling difference 
pathways.show <- c("MHC-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


pdf(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/figures/ChordPlot_batch1_13_rectum_group_mhcii_pathway_12_11_24.pdf", width = 15, height = 15)
pathways.show <- c("MHC-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

#MHC-I signaling differences 
pathways.show <- c("MHC-I") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


pathways.show <- c("BAFF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CEACAM") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("KLK") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("ADGRE") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("GALECTIN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

saveRDS(cellchat, file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/cell_chat/rectum/data/helm_batch1_13_rectum_group_comb_cellchat_obj_12_11_24.rds")

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("G1", "G2")) # set factor level
plotGeneExpression(cellchat, signaling = c("TNF"), split.by = "datasets", type = "violin")

