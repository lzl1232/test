library(Seurat)
library(tidyverse)
library(tidyverse)
library(Seurat)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(magrittr)
library(Seurat)
library(CellChat)
library(tidyverse)
library(clustree)
library(SeuratData)
library(patchwork)
load("EYE.Rdata")


####NC####
stim.object <- subset(EYE , stim == "NC")
stim.data.input <- GetAssayData(stim.object , assay = "RNA" , slot = "data")
stim.meta <- stim.object@meta.data[ , c("zhushi" , "stim")]
stim.meta$zhushi %<>%as.vector(.)

stim.cellchat <- createCellChat(object = stim.data.input)
stim.cellchat <- addMeta(stim.cellchat , meta = stim.meta)
stim.cellchat <- setIdent(stim.cellchat , ident.use = "zhushi")

levels(stim.cellchat@idents)
groupsize <- as.numeric(table(stim.cellchat@idents))

stim.cellchat@DB <- CellChatDB.human

dplyr::glimpse(CellChatDB.human$interaction)


stim.cellchat <- subsetData(stim.cellchat , features = NULL)
future::plan("multisession" , workers = 10)
stim.cellchat <- identifyOverExpressedGenes(stim.cellchat)
stim.cellchat <- identifyOverExpressedInteractions(stim.cellchat)


stim.cellchat <- projectData(stim.cellchat , PPI.human)

stim.cellchat <- computeCommunProb(stim.cellchat , raw.use = T)
stim.cellchat <- filterCommunication(stim.cellchat , min.cells = 10)
stim.cellchat <- computeCommunProbPathway(stim.cellchat)
stim.cellchat <- aggregateNet(stim.cellchat)
stim.cellchat <- netAnalysis_computeCentrality(stim.cellchat , slot.name = "netP")

group1.net <- subsetCommunication(stim.cellchat)
write.csv(group1.net , file = "group1_net_inter.csv")
saveRDS(stim.cellchat , file = "NC_cellchat.Rdata")

groupsize <- as.numeric(table(stim.cellchat@idents))

#整体细胞互作count
?netVisual_circle
netVisual_circle(stim.cellchat@net$count ,
                 vertex.weight = groupsize ,
                 weight.scale = T ,
                 label.edge = F ,
                 title.name = "Number of interactions" ,
                 color.use = c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b"))
#整体细胞互作weight
netVisual_circle(stim.cellchat@net$weight ,
                 vertex.weight = groupsize ,
                 weight.scale = T ,
                 label.edge = F ,
                 title.name = "Interaction weights/strength" ,
                 color.use = c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b"))
#单个细胞互作count
dev.off()
mat <- stim.cellchat@net$count
par(mfrow = c(3 , 4) , xpd = T)
for (i in 1 : nrow(mat)) {
  mat2 <- matrix(0 , nrow = nrow(mat) , ncol = ncol(mat) , dimnames = dimnames(mat))
  mat2[i , ] <- mat[i , ]
  netVisual_circle(mat2 , vertex.weight = groupsize , weight.scale = T , arrow.width = 0.2 ,
                   arrow.size = 0.1 , edge.weight.max = max(mat) , title.name = rownames(mat)[i])

}  ## 2000 * 1000
#单个细胞互作weight
dev.off()
mat <- stim.cellchat@net$weight
par(mfrow = c(3 , 4) , xpd = T)
for (i in 1 : nrow(mat)) {
  mat2 <- matrix(0 , nrow = nrow(mat) , ncol = ncol(mat) , dimnames = dimnames(mat))
  mat2[i , ] <- mat[i , ]
  netVisual_circle(mat2 , vertex.weight = groupsize , weight.scale = T , arrow.width = 0.2 ,
                   arrow.size = 0.1 , edge.weight.max = max(mat) , title.name = rownames(mat)[i])

}  ## 2000 * 1000


levels(stim.cellchat@idents)
stim.cellchat@netP$pathways
##
pathshow <- c("CNTN")
vertex.receiver = c(1 , 2 , 4, 8)
netVisual_aggregate(stim.cellchat , signaling = pathshow , vertex.receiver = vertex.receiver , layout = "hierarchy")
dev.off()
netVisual_aggregate(stim.cellchat , signaling = pathshow , layout = "circle")

netVisual_heatmap(stim.cellchat , signaling = pathshow , color.heatmap = "Reds")

netAnalysis_contribution(stim.cellchat , signaling = pathshow)

pairLR.PTN <- extractEnrichedLR(stim.cellchat , signaling = pathshow , geneLR.return = F)

LR.show <- pairLR.PTN[1 , ]
vertex.receiver = c(1 , 2 , 4, 8)
netVisual_individual(stim.cellchat , signaling = pathshow , pairLR.use = LR.show , vertex.receiver = vertex.receiver , layout = "hierarchy")
netVisual_individual(stim.cellchat , signaling = pathshow , pairLR.use = LR.show , vertex.receiver = vertex.receiver , layout = "circle")
##


netVisual_bubble(stim.cellchat ,
                 sources.use = "Muller" ,   ###指定来源
                 #targets.use = "" ,           ###指定去向
                 remove.isolate = F ,
                 #signaling = ""               ###指定路线
                 )
#pairLR.PTN <- extractEnrichedLR(stim.cellchat , signaling = pathshow , geneLR.return = F)   ###指定小通路


p <- plotGeneExpression(stim.cellchat , signaling = "CNTN1")
p
####LIM####
stim.object_LIM <- subset(EYE , stim == "LIM")
stim.data.input_LIM <- GetAssayData(stim.object_LIM , assay = "RNA" , slot = "data")
stim.meta_LIM <- stim.object_LIM@meta.data[ , c("zhushi" , "stim")]
stim.meta_LIM$zhushi %<>%as.vector(.)

stim.cellchat_LIM <- createCellChat(object = stim.data.input_LIM)
stim.cellchat_LIM <- addMeta(stim.cellchat_LIM , meta = stim.meta_LIM)
stim.cellchat_LIM <- setIdent(stim.cellchat_LIM , ident.use = "zhushi")

levels(stim.cellchat_LIM@idents)
groupsize_LIM <- as.numeric(table(stim.cellchat_LIM@idents))

stim.cellchat_LIM@DB <- CellChatDB.human

dplyr::glimpse(CellChatDB.human$interaction)


stim.cellchat_LIM <- subsetData(stim.cellchat_LIM , features = NULL)
future::plan("multisession" , workers = 10)
stim.cellchat_LIM <- identifyOverExpressedGenes(stim.cellchat_LIM)
stim.cellchat_LIM <- identifyOverExpressedInteractions(stim.cellchat_LIM)


stim.cellchat_LIM <- projectData(stim.cellchat_LIM , PPI.human)

stim.cellchat_LIM <- computeCommunProb(stim.cellchat_LIM , raw.use = T)
stim.cellchat_LIM <- filterCommunication(stim.cellchat_LIM , min.cells = 10)
stim.cellchat_LIM <- computeCommunProbPathway(stim.cellchat_LIM)
stim.cellchat_LIM <- aggregateNet(stim.cellchat_LIM)
stim.cellchat_LIM <- netAnalysis_computeCentrality(stim.cellchat_LIM , slot.name = "netP")

group2.net <- subsetCommunication(stim.cellchat_LIM)
write.csv(group1.net , file = "group1_net_inter.csv")
saveRDS(stim.cellchat_LIM , file = "LIM.cellchat.Rdata")

groupsize_LIM <- as.numeric(table(stim.cellchat_LIM@idents))

#整体细胞互作count
netVisual_circle(stim.cellchat_LIM@net$count ,
                 vertex.weight = groupsize ,
                 weight.scale = T ,
                 label.edge = F ,
                 title.name = "Number of interactions" ,
                 color.use = c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b"))
#整体细胞互作weight
netVisual_circle(stim.cellchat_LIM@net$weight ,
                 vertex.weight = groupsize ,
                 weight.scale = T ,
                 label.edge = F ,
                 title.name = "Interaction weights/strength" ,
                 color.use = c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b"))
#单个细胞互作count
dev.off()
mat <- stim.cellchat_LIM@net$count
par(mfrow = c(3 , 4) , xpd = T)
for (i in 1 : nrow(mat)) {
  mat2 <- matrix(0 , nrow = nrow(mat) , ncol = ncol(mat) , dimnames = dimnames(mat))
  mat2[i , ] <- mat[i , ]
  netVisual_circle(mat2 , vertex.weight = groupsize , weight.scale = T , arrow.width = 0.2 ,
                   arrow.size = 0.1 , edge.weight.max = max(mat) , title.name = rownames(mat)[i])

}
#单个细胞互作weight
dev.off()
mat <- stim.cellchat_LIM@net$weight
par(mfrow = c(3 , 4) , xpd = T)
for (i in 1 : nrow(mat)) {
  mat2 <- matrix(0 , nrow = nrow(mat) , ncol = ncol(mat) , dimnames = dimnames(mat))
  mat2[i , ] <- mat[i , ]
  netVisual_circle(mat2 , vertex.weight = groupsize , weight.scale = T , arrow.width = 0.2 ,
                   arrow.size = 0.1 , edge.weight.max = max(mat) , title.name = rownames(mat)[i])

}


levels(stim.cellchat_LIM@idents)
stim.cellchat_LIM@netP$pathways
##
pathshow <- c("CNTN")
vertex.receiver = c(1 , 2 , 4, 8)
netVisual_aggregate(stim.cellchat_LIM , signaling = pathshow , vertex.receiver = vertex.receiver , layout = "hierarchy")
dev.off()
netVisual_aggregate(stim.cellchat_LIM , signaling = pathshow , layout = "circle")

netVisual_heatmap(stim.cellchat_LIM , signaling = pathshow , color.heatmap = "Reds")

netAnalysis_contribution(stim.cellchat_LIM , signaling = pathshow)

pairLR.PTN <- extractEnrichedLR(stim.cellchat_LIM , signaling = pathshow , geneLR.return = F)

LR.show <- pairLR.PTN[1 , ]
vertex.receiver = c(1 , 2 , 4, 8)
netVisual_individual(stim.cellchat_LIM , signaling = pathshow , pairLR.use = LR.show , vertex.receiver = vertex.receiver , layout = "hierarchy")
netVisual_individual(stim.cellchat_LIM , signaling = pathshow , pairLR.use = LR.show , vertex.receiver = vertex.receiver , layout = "circle")
##


netVisual_bubble(stim.cellchat_LIM ,
                 sources.use = "Muller" ,   ###指定来源
                 #targets.use = "" ,           ###指定去向
                 remove.isolate = F ,
                 #signaling = ""               ###指定路线
)
#pairLR.PTN <- extractEnrichedLR(stim.cellchat , signaling = pathshow , geneLR.return = F)   ###指定小通路


p <- plotGeneExpression(stim.cellchat_LIM , signaling = "PTN")
p

####EA####
stim.object_EA <- subset(EYE , stim == "EA")
stim.data.input_EA <- GetAssayData(stim.object_EA , assay = "RNA" , slot = "data")
stim.meta_EA <- stim.object_EA@meta.data[ , c("zhushi" , "stim")]
stim.meta_EA$zhushi %<>%as.vector(.)

stim.cellchat_EA <- createCellChat(object = stim.data.input_EA)
stim.cellchat_EA <- addMeta(stim.cellchat_EA , meta = stim.meta_EA)
stim.cellchat_EA <- setIdent(stim.cellchat_EA , ident.use = "zhushi")

levels(stim.cellchat_EA@idents)
groupsize_EA <- as.numeric(table(stim.cellchat_EA@idents))

stim.cellchat_EA@DB <- CellChatDB.human

dplyr::glimpse(CellChatDB.human$interaction)


stim.cellchat_EA <- subsetData(stim.cellchat_EA , features = NULL)
future::plan("multisession" , workers = 10)
stim.cellchat_EA <- identifyOverExpressedGenes(stim.cellchat_EA)
stim.cellchat_EA <- identifyOverExpressedInteractions(stim.cellchat_EA)


stim.cellchat_EA <- projectData(stim.cellchat_EA , PPI.human)

stim.cellchat_EA <- computeCommunProb(stim.cellchat_EA , raw.use = T)
stim.cellchat_EA <- filterCommunication(stim.cellchat_EA , min.cells = 10)
stim.cellchat_EA <- computeCommunProbPathway(stim.cellchat_EA)
stim.cellchat_EA <- aggregateNet(stim.cellchat_EA)
stim.cellchat_EA <- netAnalysis_computeCentrality(stim.cellchat_EA , slot.name = "netP")

group3.net <- subsetCommunication(stim.cellchat_EA)
write.csv(group1.net , file = "group1_net_inter.csv")
saveRDS(stim.cellchat , file = "stim.cellchat.Rdata")

groupsize_EA <- as.numeric(table(stim.cellchat_EA@idents))

#整体细胞互作count
netVisual_circle(stim.cellchat_EA@net$count ,
                 vertex.weight = groupsize ,
                 weight.scale = T ,
                 label.edge = F ,
                 title.name = "Number of interactions" ,
                 color.use = c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b"))
#整体细胞互作weight
netVisual_circle(stim.cellchat_EA@net$weight ,
                 vertex.weight = groupsize ,
                 weight.scale = T ,
                 label.edge = F ,
                 title.name = "Interaction weights/strength" ,
                 color.use = c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b"))
#单个细胞互作count
dev.off()
mat <- stim.cellchat_EA@net$count
par(mfrow = c(3 , 4) , xpd = T)
for (i in 1 : nrow(mat)) {
  mat2 <- matrix(0 , nrow = nrow(mat) , ncol = ncol(mat) , dimnames = dimnames(mat))
  mat2[i , ] <- mat[i , ]
  netVisual_circle(mat2 , vertex.weight = groupsize , weight.scale = T , arrow.width = 0.2 ,
                   arrow.size = 0.1 , edge.weight.max = max(mat) , title.name = rownames(mat)[i])

}
#单个细胞互作weight
dev.off()
mat <- stim.cellchat_EA@net$weight
par(mfrow = c(3 , 4) , xpd = T)
for (i in 1 : nrow(mat)) {
  mat2 <- matrix(0 , nrow = nrow(mat) , ncol = ncol(mat) , dimnames = dimnames(mat))
  mat2[i , ] <- mat[i , ]
  netVisual_circle(mat2 , vertex.weight = groupsize , weight.scale = T , arrow.width = 0.2 ,
                   arrow.size = 0.1 , edge.weight.max = max(mat) , title.name = rownames(mat)[i])

}


levels(stim.cellchat_EA@idents)
stim.cellchat_EA@netP$pathways
##
pathshow <- c("CNTN")
vertex.receiver = c(1 , 2 , 8)
netVisual_aggregate(stim.cellchat_EA , signaling = pathshow , vertex.receiver = vertex.receiver , layout = "hierarchy")
dev.off()
netVisual_aggregate(stim.cellchat_EA , signaling = pathshow , layout = "circle")

netVisual_heatmap(stim.cellchat_EA , signaling = pathshow , color.heatmap = "Reds")

netAnalysis_contribution(stim.cellchat_EA , signaling = pathshow)

pairLR.PTN <- extractEnrichedLR(stim.cellchat_EA , signaling = pathshow , geneLR.return = F)

LR.show <- pairLR.PTN[1 , ]
vertex.receiver = c(1 , 2 , 8)
netVisual_individual(stim.cellchat_EA , signaling = pathshow , pairLR.use = LR.show , vertex.receiver = vertex.receiver , layout = "hierarchy")
netVisual_individual(stim.cellchat_EA , signaling = pathshow , pairLR.use = LR.show , vertex.receiver = vertex.receiver , layout = "circle")
##


netVisual_bubble(stim.cellchat_EA ,
                 sources.use = "Muller" ,   ###指定来源
                 #targets.use = "" ,           ###指定去向
                 remove.isolate = F ,
                 #signaling = ""               ###指定路线
)
#pairLR.PTN <- extractEnrichedLR(stim.cellchat , signaling = pathshow , geneLR.return = F)   ###指定小通路


p <- plotGeneExpression(stim.cellchat_EA , signaling = "PTN")
p

####合并-NC-LIM####
cco.list <- list(NC = stim.cellchat , LIM = stim.cellchat_LIM)
cellchat <- mergeCellChat(cco.list , add.names = names(cco.list) , cell.prefix = T)
?compareInteractions
compareInteractions(cellchat , show.legend = F , group = c(1 : 2) , measure = "count" , color.use = c("#a2d2e7" , "#ff9d9f"))
compareInteractions(cellchat , show.legend = F , group = c(1 : 2) , measure = "weight", color.use = c("#a2d2e7" , "#ff9d9f"))

par(mfrow = c(1 : 2))
?netVisual_diffInteraction
netVisual_diffInteraction(cellchat , weight.scale = T , color.use =   c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b") ,
                          color.edge = c( "#a2d2e7", "#ff9d9f"))
netVisual_diffInteraction(cellchat , weight.scale = T , measure = "weight" , color.use =   c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b") ,
                          color.edge = c( "#a2d2e7" ,"#ff9d9f" ))

par(mfrow = c(1 , 1))
?netVisual_heatmap
h1 <- netVisual_heatmap(cellchat , color.heatmap =  c( "#a2d2e7" , "#ff9d9f"))
h2 <- netVisual_heatmap(cellchat , measure = "weight" ,color.heatmap =  c(  "#a2d2e7" , "#ff9d9f"))
h1 + h2

gg1 <- rankNet(cellchat , mode = "comparison" , stacked = T , do.stat = T , color.use = c( "#a2d2e7" , "#ff9d9f"))
#gg2 <- rankNet(cellchat , mode = "comparison" , stacked = F , do.stat = T , color.use = c("#a2d2e7" , "#ff9d9f"))
gg1

levels(cellchat@idents$joint)
netVisual_bubble(cellchat ,
                 sources.use = c(1 , 2) ,
                 targets.use = c(3 : 11) ,
                 comparison = c(1 , 2) ,
                 angle.x = 45)


####合并-EA-LIM####
cco.list1 <- list(LIM = stim.cellchat_LIM , EA = stim.cellchat_EA)
cellchat1 <- mergeCellChat(cco.list1 , add.names = names(cco.list1) , cell.prefix = T)

compareInteractions(cellchat1 , show.legend = F , group = c(1 : 2) , measure = "count" , c("#ff9d9f" , "#6fb3a3"))
compareInteractions(cellchat1 , show.legend = F , group = c(1 : 2) , measure = "weight" , c("#ff9d9f", "#6fb3a3"))

par(mfrow = c(1 : 2))
netVisual_diffInteraction(cellchat1 , weight.scale = T ,color.use =  c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b") ,
                          color.edge = c("#ff9d9f" , "#6fb3a3"))
netVisual_diffInteraction(cellchat1 , weight.scale = T , measure = "weight" , color.use =  c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b") ,
                          color.edge = c( "#ff9d9f" , "#6fb3a3"))

par(mfrow = c(1 , 1))
h1 <- netVisual_heatmap(cellchat1 , color.heatmap = c( "#ff9d9f" , "#6fb3a3"))
h2 <- netVisual_heatmap(cellchat1 , measure = "weight" , color.heatmap = c("#ff9d9f" , "#6fb3a3"))
h1 + h2

gg1 <- rankNet(cellchat1 , mode = "comparison" , stacked = T , do.stat = T , color.use = c("#ff9d9f" , "#6fb3a3"))
#gg2 <- rankNet(cellchat1 , mode = "comparison" , stacked = F , do.stat = T)
gg1 #+ gg2

levels(cellchat1@idents$joint)
netVisual_bubble(cellchat1 ,
                 sources.use = c(1:9) ,
                 targets.use = c(6) ,
                 comparison = c(1 , 2) ,
                 angle.x = 45)
