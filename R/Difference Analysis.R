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
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(scRNAtoolVis)
library(tidyverse)
library(ggVolcano)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(data.table)

load("EYE.Rdata")
Idents(EYE) <- "zhushi"
Idents(EYE)
Bipolar <- subset(EYE  , idents = "Bipolar")
Cone <- subset(EYE  , idents = "Cone")
Muller <- subset(EYE , idents = "Muller")
Rod <- subset(EYE , idents = "Rod")
Microglia <- subset(EYE , idents = "Microglia")
Ganglion <- subset(EYE , idents = "Ganglion")
Neuron <- subset(EYE , idents = "Neuron")
Rpe <- subset(EYE , idents = "Rpe")
Amacrine <- subset(EYE , idents = "Amacrine")


table(EYE$zhushi)

####标准化####
Muller <- NormalizeData(Muller , normalization.method = "LogNormalize" , scale.factor = 10000)
Muller <- FindVariableFeatures(Muller , selection.method = "vst" , nfeatures = 2000)
all.genes <- rownames(Muller)
Muller <- ScaleData(Muller , features = all.genes)
Muller <- RunPCA(Muller , features = VariableFeatures(object = Muller))
Muller <- FindNeighbors(Muller , dims = 1:30 )
Muller <- FindClusters(Muller , resolution = 0.5)
Muller <- RunUMAP(Muller , dims = 1:30)


Bipolar <- NormalizeData(Bipolar , normalization.method = "LogNormalize" , scale.factor = 10000)
Bipolar <- FindVariableFeatures(Bipolar , selection.method = "vst" , nfeatures = 2000)
all.genes <- rownames(Bipolar)
Bipolar <- ScaleData(Bipolar , features = all.genes)
Bipolar <- RunPCA(Bipolar , features = VariableFeatures(object = Bipolar))
Bipolar <- FindNeighbors(Bipolar , dims = 1:30 )
Bipolar <- FindClusters(Bipolar , resolution = 0.5)
Bipolar <- RunUMAP(Bipolar , dims = 1:30)

Cone <- NormalizeData(Cone , normalization.method = "LogNormalize" , scale.factor = 10000)
Cone <- FindVariableFeatures(Cone , selection.method = "vst" , nfeatures = 2000)
all.genes <- rownames(Cone)
Cone <- ScaleData(Cone , features = all.genes)
Cone <- RunPCA(Cone , features = VariableFeatures(object = Cone))
Cone <- FindNeighbors(Cone , dims = 1:30 )
Cone <- FindClusters(Cone , resolution = 0.5)
Cone <- RunUMAP(Cone , dims = 1:30)

Rod <- NormalizeData(Rod , normalization.method = "LogNormalize" , scale.factor = 10000)
Rod <- FindVariableFeatures(Rod , selection.method = "vst" , nfeatures = 2000)
all.genes <- rownames(Rod)
Rod <- ScaleData(Rod , features = all.genes)
Rod <- RunPCA(Rod , features = VariableFeatures(object = Rod))
Rod <- FindNeighbors(Rod , dims = 1:30 )
Rod <- FindClusters(Rod , resolution = 0.5)
Rod <- RunUMAP(Rod , dims = 1:30)

Microglia <- NormalizeData(Microglia , normalization.method = "LogNormalize" , scale.factor = 10000)
Microglia <- FindVariableFeatures(Microglia , selection.method = "vst" , nfeatures = 2000)
all.genes <- rownames(Microglia)
Microglia <- ScaleData(Microglia , features = all.genes)
Microglia <- RunPCA(Microglia , features = VariableFeatures(object = Microglia))
Microglia <- FindNeighbors(Microglia , dims = 1:30 )
Microglia <- FindClusters(Microglia , resolution = 0.5)
Microglia <- RunUMAP(Microglia , dims = 1:30)

Ganglion <- NormalizeData(Ganglion , normalization.method = "LogNormalize" , scale.factor = 10000)
Ganglion <- FindVariableFeatures(Ganglion , selection.method = "vst" , nfeatures = 2000)
all.genes <- rownames(Ganglion)
Ganglion <- ScaleData(Ganglion , features = all.genes)
Ganglion <- RunPCA(Ganglion , features = VariableFeatures(object = Ganglion))
Ganglion <- FindNeighbors(Ganglion , dims = 1:30 )
Ganglion <- FindClusters(Ganglion , resolution = 0.5)
Ganglion <- RunUMAP(Ganglion , dims = 1:30)

Neuron <- NormalizeData(Neuron , normalization.method = "LogNormalize" , scale.factor = 10000)
Neuron <- FindVariableFeatures(Neuron , selection.method = "vst" , nfeatures = 2000)
all.genes <- rownames(Neuron)
Neuron <- ScaleData(Neuron , features = all.genes)
Neuron <- RunPCA(Neuron , features = VariableFeatures(object = Neuron))
Neuron <- FindNeighbors(Neuron , dims = 1:30 )
Neuron <- FindClusters(Neuron , resolution = 0.5)
Neuron <- RunUMAP(Neuron , dims = 1:30)

Rpe <- NormalizeData(Rpe , normalization.method = "LogNormalize" , scale.factor = 10000)
Rpe <- FindVariableFeatures(Rpe , selection.method = "vst" , nfeatures = 2000)
all.genes <- rownames(Rpe)
Rpe <- ScaleData(Rpe , features = all.genes)
Rpe <- RunPCA(Rpe , features = VariableFeatures(object = Rpe))
Rpe <- FindNeighbors(Rpe , dims = 1:30 )
Rpe <- FindClusters(Rpe , resolution = 0.5)
Rpe <- RunUMAP(Rpe , dims = 1:30)

Amacrine <- NormalizeData(Amacrine , normalization.method = "LogNormalize" , scale.factor = 10000)
Amacrine <- FindVariableFeatures(Amacrine , selection.method = "vst" , nfeatures = 2000)
all.genes <- rownames(Amacrine)
Amacrine <- ScaleData(Amacrine , features = all.genes)
Amacrine <- RunPCA(Amacrine , features = VariableFeatures(object = Amacrine))
Amacrine <- FindNeighbors(Amacrine , dims = 1:30 )
Amacrine <- FindClusters(Amacrine , resolution = 0.5)
Amacrine <- RunUMAP(Amacrine , dims = 1:30)


####差异基因####
Idents(Rod) = "stim"
Idents(Bipolar) = "stim"
Idents(Cone) = "stim"
Idents(Muller) = "stim"
Idents(Microglia) = "stim"
Idents(Ganglion) = "stim"
Idents(Neuron) = "stim"
Idents(Rpe) = "stim"
Idents(Amacrine) = "stim"


marker_Rod1 <- FindMarkers(Rod , ident.1 = "LIM" , ident.2 = "NC")
marker_Rod2 <- FindMarkers(Rod , ident.1 = "EA" , ident.2 = "LIM")

marker_Bipolar1 <- FindMarkers(Bipolar , ident.1 = "LIM" , ident.2 = "NC")
marker_Bipolar2 <- FindMarkers(Bipolar , ident.1 = "EA" , ident.2 = "LIM")

marker_Cone1 <- FindMarkers(Cone , ident.1 = "LIM" , ident.2 = "NC")
marker_Cone2 <- FindMarkers(Cone , ident.1 = "EA" , ident.2 = "LIM")

marker_Microglia1 <- FindMarkers(Microglia , ident.1 = "LIM" , ident.2 = "NC")
marker_Microglia2 <- FindMarkers(Microglia , ident.1 = "EA" , ident.2 = "LIM")

marker_Neuron1 <- FindMarkers(Neuron , ident.1 = "LIM" , ident.2 = "NC")
marker_Neuron2 <- FindMarkers(Neuron , ident.1 = "EA" , ident.2 = "LIM")

marker_Ganglion1 <- FindMarkers(Ganglion , ident.1 = "LIM" , ident.2 = "NC")
marker_Ganglion2 <- FindMarkers(Ganglion , ident.1 = "EA" , ident.2 = "LIM")

marker_Rpe1 <- FindMarkers(Rpe , ident.1 = "LIM" , ident.2 = "NC")
marker_Rpe2 <- FindMarkers(Rpe , ident.1 = "EA" , ident.2 = "LIM")

marker_Amacrine1 <- FindMarkers(Amacrine , ident.1 = "LIM" , ident.2 = "NC")
marker_Amacrine2 <- FindMarkers(Amacrine , ident.1 = "EA" , ident.2 = "LIM")

marker_Muller1 <- FindMarkers(Muller , ident.1 = "LIM" , ident.2 = "NC")
marker_Muller2 <- FindMarkers(Muller , ident.1 = "EA" , ident.2 = "LIM")

#
marker_Rod1$gene = rownames(marker_Rod1)
marker_Rod2$gene = rownames(marker_Rod2)
mer_Rod <- merge(marker_Rod1 , marker_Rod2 , by = "gene")

marker_Bipolar1$gene = rownames(marker_Bipolar1)
marker_Bipolar2$gene = rownames(marker_Bipolar2)
mer_Bipolar <- merge(marker_Bipolar1 , marker_Bipolar2 , by = "gene")

marker_Cone1$gene = rownames(marker_Cone1)
marker_Cone2$gene = rownames(marker_Cone2)
mer_Cone <- merge(marker_Cone1 , marker_Cone2 , by = "gene")

marker_Microglia1$gene = rownames(marker_Microglia1)
marker_Microglia2$gene = rownames(marker_Microglia2)
mer_Microglia <- merge(marker_Microglia1 , marker_Microglia2 , by = "gene")

marker_Neuron1$gene = rownames(marker_Neuron1)
marker_Neuron2$gene = rownames(marker_Neuron2)
mer_Neuron <- merge(marker_Neuron1 , marker_Neuron2 , by = "gene")

marker_Ganglion1$gene = rownames(marker_Ganglion1)
marker_Ganglion2$gene = rownames(marker_Ganglion2)
mer_Ganglion <- merge(marker_Ganglion1 , marker_Ganglion2 , by = "gene")

marker_Rpe1$gene = rownames(marker_Rpe1)
marker_Rpe2$gene = rownames(marker_Rpe2)
mer_Rpe <- merge(marker_Rpe1 , marker_Rpe2 , by = "gene")

marker_Amacrine1$gene = rownames(marker_Amacrine1)
marker_Amacrine2$gene = rownames(marker_Amacrine2)
mer_Amacrine <- merge(marker_Amacrine1 , marker_Amacrine2 , by = "gene")

marker_Muller1$gene = rownames(marker_Muller1)
marker_Muller2$gene = rownames(marker_Muller2)
mer_Muller <- merge(marker_Muller1 , marker_Muller2 , by = "gene")

write.csv(marker_Amacrine2 , file = "marker_Amarcine2.csv")

####定义up ， down####
index <- grep("^ENSCPOG" , marker_Amacrine1$gene)
marker_Amacrine1 <- marker_Amacrine1[-index ,]
marker_Amacrine1 $ logp <- -log10(marker_Amacrine1$p_val)
marker_Amacrine1$group <- as.factor(ifelse(
  marker_Amacrine1$p_val <0.05 & abs(marker_Amacrine1$avg_log2FC)>=0.5,
  ifelse(marker_Amacrine1$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Amacrine1$group)


index <- grep("^ENSCPOG" , marker_Amacrine2$gene)
marker_Amacrine2 <- marker_Amacrine2[-index ,]
marker_Amacrine2 $ logp <- -log10(marker_Amacrine2$p_val)
marker_Amacrine2$group <- as.factor(ifelse(
  marker_Amacrine2$p_val <0.05 & abs(marker_Amacrine2$avg_log2FC)>=0.5,
  ifelse(marker_Amacrine2$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Amacrine2$group)


index <- grep("^ENSCPOG" , marker_Bipolar1$gene)
marker_Bipolar1 <- marker_Bipolar1[-index ,]
marker_Bipolar1 $ logp <- -log10(marker_Bipolar1$p_val)
marker_Bipolar1$group <- as.factor(ifelse(
  marker_Bipolar1$p_val <0.05 & abs(marker_Bipolar1$avg_log2FC)>=0.5,
  ifelse(marker_Bipolar1$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Bipolar1$group)


index <- grep("^ENSCPOG" , marker_Bipolar2$gene)
marker_Bipolar2 <- marker_Bipolar2[-index ,]
marker_Bipolar2$ logp <- -log10(marker_Bipolar2$p_val)
marker_Bipolar2$group <- as.factor(ifelse(
  marker_Bipolar2$p_val <0.05 & abs(marker_Bipolar2$avg_log2FC)>=0.5,
  ifelse(marker_Bipolar2$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Bipolar2$group)


index <- grep("^ENSCPOG" , marker_Cone1$gene)
marker_Cone1 <- marker_Cone1[-index ,]
marker_Cone1$ logp <- -log10(marker_Cone1$p_val)
marker_Cone1$group <- as.factor(ifelse(
  marker_Cone1$p_val <0.05 & abs(marker_Cone1$avg_log2FC)>=0.5,
  ifelse(marker_Cone1$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Cone1$group)


index <- grep("^ENSCPOG" , marker_Cone2$gene)
marker_Cone2 <- marker_Cone2[-index ,]
marker_Cone2$ logp <- -log10(marker_Cone2$p_val)
marker_Cone2$group <- as.factor(ifelse(
  marker_Cone2$p_val <0.05 & abs(marker_Cone2$avg_log2FC)>=0.5,
  ifelse(marker_Cone2$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Cone2$group)


index <- grep("^ENSCPOG" , marker_Ganglion1$gene)
marker_Ganglion1 <- marker_Ganglion1[-index ,]
marker_Ganglion1$ logp <- -log10(marker_Ganglion1$p_val)
marker_Ganglion1$group <- as.factor(ifelse(
  marker_Ganglion1$p_val <0.05 & abs(marker_Ganglion1$avg_log2FC)>=0.5,
  ifelse(marker_Ganglion1$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Ganglion1$group)


index <- grep("^ENSCPOG" , marker_Ganglion2$gene)
marker_Ganglion2 <- marker_Ganglion2[-index ,]
marker_Ganglion2$ logp <- -log10(marker_Ganglion2$p_val)
marker_Ganglion2$group <- as.factor(ifelse(
  marker_Ganglion2$p_val <0.05 & abs(marker_Ganglion2$avg_log2FC)>=0.5,
  ifelse(marker_Ganglion2$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Ganglion2$group)


index <- grep("^ENSCPOG" , marker_Microglia1$gene)
marker_Microglia1 <- marker_Microglia1[-index ,]
marker_Microglia1$ logp <- -log10(marker_Microglia1$p_val)
marker_Microglia1$group <- as.factor(ifelse(
  marker_Microglia1$p_val <0.05 & abs(marker_Microglia1$avg_log2FC)>=0.5,
  ifelse(marker_Microglia1$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Microglia1$group)


index <- grep("^ENSCPOG" , marker_Microglia2$gene)
marker_Microglia2 <- marker_Microglia2[-index ,]
marker_Microglia2$ logp <- -log10(marker_Microglia2$p_val)
marker_Microglia2$group <- as.factor(ifelse(
  marker_Microglia2$p_val <0.05 & abs(marker_Microglia2$avg_log2FC)>=0.5,
  ifelse(marker_Microglia2$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Microglia2$group)


index <- grep("^ENSCPOG" , marker_Neuron1$gene)
marker_Neuron1 <- marker_Neuron1[-index ,]
marker_Neuron1$ logp <- -log10(marker_Neuron1$p_val)
marker_Neuron1$group <- as.factor(ifelse(
  marker_Neuron1$p_val <0.05 & abs(marker_Neuron1$avg_log2FC)>=0.5,
  ifelse(marker_Neuron1$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Neuron1$group)


index <- grep("^ENSCPOG" , marker_Neuron2$gene)
marker_Neuron2 <- marker_Neuron2[-index ,]
marker_Neuron2$ logp <- -log10(marker_Neuron2$p_val)
marker_Neuron2$group <- as.factor(ifelse(
  marker_Neuron2$p_val <0.05 & abs(marker_Neuron2$avg_log2FC)>=0.5,
  ifelse(marker_Neuron2$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Neuron2$group)


index <- grep("^ENSCPOG" , marker_Muller1$gene)
marker_Muller1 <- marker_Muller1[-index ,]
marker_Muller1$ logp <- -log10(marker_Muller1$p_val)
marker_Muller1$group <- as.factor(ifelse(
  marker_Muller1$p_val <0.05 & abs(marker_Muller1$avg_log2FC)>=0.5,
  ifelse(marker_Muller1$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Muller1$group)


index <- grep("^ENSCPOG" , marker_Muller2$gene)
marker_Muller2 <- marker_Muller2[-index ,]
marker_Muller2$ logp <- -log10(marker_Muller2$p_val)
marker_Muller2$group <- as.factor(ifelse(
  marker_Muller2$p_val <0.05 & abs(marker_Muller2$avg_log2FC)>=0.5,
  ifelse(marker_Muller2$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Muller2$group)


index <- grep("^ENSCPOG" , marker_Rod1$gene)
marker_Rod1 <- marker_Rod1[-index ,]
marker_Rod1$ logp <- -log10(marker_Rod1$p_val)
marker_Rod1$group <- as.factor(ifelse(
  marker_Rod1$p_val <0.05 & abs(marker_Rod1$avg_log2FC)>=0.5,
  ifelse(marker_Rod1$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Rod1$group)


index <- grep("^ENSCPOG" , marker_Rod2$gene)
marker_Rod2 <- marker_Rod2[-index ,]
marker_Rod2$ logp <- -log10(marker_Rod2$p_val)
marker_Rod2$group <- as.factor(ifelse(
  marker_Rod2$p_val <0.05 & abs(marker_Rod2$avg_log2FC)>=0.5,
  ifelse(marker_Rod2$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Rod2$group)


index <- grep("^ENSCPOG" , marker_Rpe1$gene)
marker_Rpe1 <- marker_Rpe1[-index ,]
marker_Rpe1$ logp <- -log10(marker_Rpe1$p_val)
marker_Rpe1$group <- as.factor(ifelse(
  marker_Rpe1$p_val <0.05 & abs(marker_Rpe1$avg_log2FC)>=0.5,
  ifelse(marker_Rpe1$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Rpe1$group)


index <- grep("^ENSCPOG" , marker_Rpe2$gene)
marker_Rpe2 <- marker_Rpe2[-index ,]
marker_Rpe2$ logp <- -log10(marker_Rpe2$p_val)
marker_Rpe2$group <- as.factor(ifelse(
  marker_Rpe2$p_val <0.05 & abs(marker_Rpe2$avg_log2FC)>=0.5,
  ifelse(marker_Rpe2$avg_log2FC>=0.5,"up","down"),"NS"
))
table(marker_Rpe2$group)


write.csv(marker_Amacrine2 , file = "marker_Amarcine2.csv")


####差异基因绘图####

up_data <- filter(marker_Bipolar2, group == 'up') %>%
  distinct(gene , .keep_all = TRUE) %>%
  top_n(10, logp)
down_data <- filter(marker_Bipolar2, group == 'down') %>%
  distinct(gene , .keep_all = TRUE) %>%
  top_n(10, logp)

ggplot(marker_Bipolar2,aes(avg_log2FC,logp),)+
  geom_point(aes(color=group),size=5)+
  scale_color_manual(values = c("#2a3663","#d9c5a8","#c1402b"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10,color = "black"))+
  theme(axis.text.y = element_text(size = 10,color = "black"))+
  theme(legend.text = element_text(size = 24,color = "black"),
        #legend.key.height = unit(80,"pt"),
        legend.position = "top",legend.box = "horizontal")+
  guides(fill=guide_legend(nrow = 1, byrow=T))+
  theme(legend.title = element_blank())+

  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
  geom_point(data = up_data,
             aes(x = avg_log2FC, y = logp),
             color = 'red3', size = 4.5, alpha = 0.2) +
  geom_label_repel(data = up_data,
                   aes(x = avg_log2FC, y = logp, label = gene),
                   seed = 233,
                   size = 3.5,
                   color = 'black',
                   min.segment.length = 0,
                   force = 2,
                   force_pull = 2,
                   box.padding = 0.4,
                   max.overlaps = Inf) +
  geom_point(data = down_data,
             aes(x = avg_log2FC, y = logp),
             color = 'blue4', size = 4.5, alpha = 0.2) +
  geom_label_repel(data = down_data,
                   aes(x = avg_log2FC, y = logp, label = gene),
                   seed = 233,
                   size = 3.5,
                   color = 'black',
                   min.segment.length = 0,
                   force = 2,
                   force_pull = 2,
                   box.padding = 0.4,
                   max.overlaps = Inf) +
  labs(title = "Bipolar_EA") +
  theme(plot.title = element_text(face = "bold" ,
                                  colour = "black" ,
                                  size = 24)) +
  theme(axis.title.x = element_text(face = "bold" ,
                                    colour = "black" ,
                                    size = 18)) +
  theme(axis.title.y = element_text(face = "bold" ,
                                    colour = "black" ,
                                    size = 18))

ggsave("Bipolar_LIM",width = 24,height = 24,units = "cm",dpi = 1000)

####差异基因富集####
  ###LIM###
#Amarcine#
UP <- subset(marker_Amacrine1 , marker_Amacrine1$group == "up")
DOWN <- subset(marker_Amacrine1 , marker_Amacrine1$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
              OrgDb = "org.Hs.eg.db",
              keyType = "ENTREZID",
              pvalueCutoff = 1 ,
              maxGSSize = 5000,
              minGSSize = 5,
              ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                 organism = "hsa",
                 keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-LIM-UP-Amarcine.Rdata")
save(EGO_DOWN , file = "EGO-LIM-DOWN-Amarcine.Rdata")
save(kEGG_UP , file = "KEGG-LIM-UP-Amarcine.Rdata")
save(kEGG_DOWN , file = "KEGG-LIM-DOWN-Amarcine.Rdata")


#Bipolar#
UP <- subset(marker_Bipolar2 , marker_Bipolar2$group == "up")
DOWN <- subset(marker_Bipolar2 , marker_Bipolar2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-LIM-UP-Bipolar.Rdata")
save(EGO_DOWN , file = "EGO-LIM-DOWN-Bipolar.Rdata")
save(kEGG_UP , file = "KEGG-LIM-UP-Bipolar.Rdata")
save(kEGG_DOWN , file = "KEGG-LIM-DOWN-Bipolar.Rdata")

#Cone#
UP <- subset(marker_Cone1 , marker_Cone1$group == "up")
DOWN <- subset(marker_Cone1 , marker_Cone1$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-LIM-UP-Cone.Rdata")
save(EGO_DOWN , file = "EGO-LIM-DOWN-Cone.Rdata")
save(kEGG_UP , file = "KEGG-LIM-UP-Cone.Rdata")
save(kEGG_DOWN , file = "KEGG-LIM-DOWN-Cone.Rdata")

#Ganglion#
UP <- subset(marker_Ganglion1 , marker_Ganglion1$group == "up")
DOWN <- subset(marker_Ganglion1 , marker_Ganglion1$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-LIM-UP-Ganglion.Rdata")
save(EGO_DOWN , file = "EGO-LIM-DOWN-Ganglion.Rdata")
save(kEGG_UP , file = "KEGG-LIM-UP-Ganglion.Rdata")
save(kEGG_DOWN , file = "KEGG-LIM-DOWN-Ganglion.Rdata")

#Microglia#

UP <- subset(marker_Microglia1 , marker_Microglia1$group == "up")
DOWN <- subset(marker_Microglia1 , marker_Microglia1$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-LIM-UP-Microglia.Rdata")
save(EGO_DOWN , file = "EGO-LIM-DOWN-Microglia.Rdata")
save(kEGG_UP , file = "KEGG-LIM-UP-Microglia.Rdata")
save(kEGG_DOWN , file = "KEGG-LIM-DOWN-Microglia.Rdata")


#Muller#

UP <- subset(marker_Muller2 , marker_Muller2$group == "up")
DOWN <- subset(marker_Muller2 , marker_Muller2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Muller.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Muller.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Muller.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Muller.Rdata")



#Neuron#
UP <- subset(marker_Neuron1 , marker_Neuron1$group == "up")
DOWN <- subset(marker_Neuron1 , marker_Neuron1$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-LIM-UP-Neuron.Rdata")
save(EGO_DOWN , file = "EGO-LIM-DOWN-Neuron.Rdata")
save(kEGG_UP , file = "KEGG-LIM-UP-Neuron.Rdata")
save(kEGG_DOWN , file = "KEGG-LIM-DOWN-Neuron.Rdata")


#Rod#
UP <- subset(marker_Rod1 , marker_Rod1$group == "up")
DOWN <- subset(marker_Rod1 , marker_Rod1$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-LIM-UP-Rod.Rdata")
save(EGO_DOWN , file = "EGO-LIM-DOWN-Rod.Rdata")
save(kEGG_UP , file = "KEGG-LIM-UP-Rod.Rdata")
save(kEGG_DOWN , file = "KEGG-LIM-DOWN-Rod.Rdata")


#Rpe#

UP <- subset(marker_Rpe1 , marker_Rpe1$group == "up")
DOWN <- subset(marker_Rpe1 , marker_Rpe1$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-LIM-UP-Rpe.Rdata")
save(EGO_DOWN , file = "EGO-LIM-DOWN-Rpe.Rdata")
save(kEGG_UP , file = "KEGG-LIM-UP-Rpe.Rdata")
save(kEGG_DOWN , file = "KEGG-LIM-DOWN-Rpe.Rdata")


###EA###

#Amarcine#
UP <- subset(marker_Amacrine2 , marker_Amacrine2$group == "up")
DOWN <- subset(marker_Amacrine2 , marker_Amacrine2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Amarcine.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Amarcine.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Amarcine.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Amarcine.Rdata")


#Bipolar#
UP <- subset(marker_Bipolar2 , marker_Bipolar2$group == "up")
DOWN <- subset(marker_Bipolar2 , marker_Bipolar2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Bipolar.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Bipolar.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Bipolar.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Bipolar.Rdata")

#Cone#
UP <- subset(marker_Cone2 , marker_Cone2$group == "up")
DOWN <- subset(marker_Cone2 , marker_Cone2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Cone.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Cone.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Cone.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Cone.Rdata")

#Ganglion#
UP <- subset(marker_Ganglion2 , marker_Ganglion2$group == "up")
DOWN <- subset(marker_Ganglion2 , marker_Ganglion2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Ganglion.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Ganglion.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Ganglion.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Ganglion.Rdata")

#Microglia#

UP <- subset(marker_Microglia2 , marker_Microglia2$group == "up")
DOWN <- subset(marker_Microglia2 , marker_Microglia2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Microglia.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Microglia.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Microglia.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Microglia.Rdata")


#Muller#

UP <- subset(marker_Muller2 , marker_Muller2$group == "up")
DOWN <- subset(marker_Muller2 , marker_Muller2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Muller.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Muller.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Muller.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Muller.Rdata")



#Neuron#
UP <- subset(marker_Neuron2 , marker_Neuron2$group == "up")
DOWN <- subset(marker_Neuron2 , marker_Neuron2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Neuron.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Neuron.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Neuron.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Neuron.Rdata")


#Rod#
UP <- subset(marker_Rod2 , marker_Rod2$group == "up")
DOWN <- subset(marker_Rod2 , marker_Rod2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Rod.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Rod.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Rod.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Rod.Rdata")


#Rpe#

UP <- subset(marker_Rpe2 , marker_Rpe2$group == "up")
DOWN <- subset(marker_Rpe2 , marker_Rpe2$group == "down")

ID<-bitr(UP$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_UP<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_UP<-setReadable(EGO_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_UP, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_UP<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_UP@result,x = "GeneRatio" , showCategory =20)
kEGG_UP<-setReadable(kEGG_UP,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")



ID2<-bitr(DOWN$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID2)
ID2 <- ID2[!duplicated_col , ]
gene_diff2<-ID2$ENTREZID
EGO_DOWN<-enrichGO(gene_diff2,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_DOWN<-setReadable(EGO_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_DOWN, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_DOWN<-enrichKEGG(ID2$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_DOWN , x = "GeneRatio" , color = "p.adjust" ,
        showCategory = 20 )
kEGG_DOWN<-setReadable(kEGG_DOWN,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


save(EGO_UP , file = "EGO-EA-UP-Rpe.Rdata")
save(EGO_DOWN , file = "EGO-EA-DOWN-Rpe.Rdata")
save(kEGG_UP , file = "KEGG-EA-UP-Rpe.Rdata")
save(kEGG_DOWN , file = "KEGG-EA-DOWN-Rpe.Rdata")

