library(Seurat)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(magrittr)
library(CellChat)
library(clustree)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(reshape2)
library(scRNAtoolVis)
library(ggVolcano)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(data.table)
library(ggridges)
library(ggplot2)
library(enrichplot)



#Amarcine#
Amarcine_LIM <- marker_Amacrine1$gene
Amarcine_LIM <- as.data.frame(Amarcine_LIM)
entrez<- bitr(Amarcine_LIM$Amarcine_LIM ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Amacrine1$avg_log2FC
names(genelist) <- marker_Amacrine1$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Amarcine-LIM.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Amarcine-LIM')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(1,4,10) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Amarcine-LIM")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


Amarcine_EA <- marker_Amacrine2$gene
Amarcine_EA <- as.data.frame(Amarcine_EA)
entrez<- bitr(Amarcine_EA$Amarcine_EA ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Amacrine2$avg_log2FC
names(genelist) <- marker_Amacrine2$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Amarcine-EA.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Amarcine-EA')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(2,4,8) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Amarcine-EA")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P



#Bipolar#
Bipolar_LIM <- marker_Bipolar1$gene
Bipolar_LIM <- as.data.frame(Bipolar_LIM)
entrez<- bitr(Bipolar_LIM$Bipolar_LIM ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Bipolar1$avg_log2FC
names(genelist) <- marker_Bipolar1$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Bipolar-LIM.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Bipolar-LIM')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(1,2,5) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Bipolar-LIM")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


Bipolar_EA <- marker_Bipolar2$gene
Bipolar_EA <- as.data.frame(Bipolar_EA)
entrez<- bitr(Bipolar_EA$Bipolar_EA ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Bipolar2$avg_log2FC
names(genelist) <- marker_Bipolar2$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Bipolar-EA.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Bipolar-EA')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(4,2,3) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Bipolar-EA")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


#Cone#
Cone_LIM <- marker_Cone1$gene
Cone_LIM <- as.data.frame(Cone_LIM)
entrez<- bitr(Cone_LIM$Cone_LIM ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Cone1$avg_log2FC
names(genelist) <- marker_Cone1$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Cone-LIM.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Cone-LIM')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(11,2,23) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Cone-LIM")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


Cone_EA <- marker_Cone2$gene
Cone_EA <- as.data.frame(Cone_EA)
entrez<- bitr(Cone_EA$Cone_EA ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Cone2$avg_log2FC
names(genelist) <- marker_Cone2$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Cone-EA.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Cone-EA')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(1,5,8) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Cone-EA")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


#Ganglion#
Ganglion_LIM <- marker_Ganglion1$gene
Ganglion_LIM <- as.data.frame(Ganglion_LIM)
entrez<- bitr(Ganglion_LIM$Ganglion_LIM ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Ganglion1$avg_log2FC
names(genelist) <- marker_Ganglion1$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Ganglion-LIM.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Ganglion-LIM')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(14,2,25) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Ganglion-LIM")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


Ganglion_EA <- marker_Ganglion2$gene
Ganglion_EA <- as.data.frame(Ganglion_EA)
entrez<- bitr(Ganglion_EA$Ganglion_EA ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Ganglion2$avg_log2FC
names(genelist) <- marker_Ganglion2$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Ganglion-EA.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Ganglion-EA')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(2,3,21) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Ganglion-EA")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


#Microglia#

Microglia_LIM <- marker_Microglia1$gene
Microglia_LIM <- as.data.frame(Microglia_LIM)
entrez<- bitr(Microglia_LIM$Microglia_LIM ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Microglia1$avg_log2FC
names(genelist) <- marker_Microglia1$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Microglia-LIM.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Microglia-LIM')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(2,1,5) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Microglia-LIM")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P



Microglia_EA <- marker_Microglia2$gene
Microglia_EA <- as.data.frame(Microglia_EA)
entrez<- bitr(Microglia_EA$Microglia_EA ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Microglia2$avg_log2FC
names(genelist) <- marker_Microglia2$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Microglia-EA.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Microglia-EA')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(5,8,14) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Microglia-EA")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


#Muller#

Muller_LIM <- marker_Muller1$gene
Muller_LIM <- as.data.frame(Muller_LIM)
entrez<- bitr(Muller_LIM$Muller_LIM ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Muller1$avg_log2FC
names(genelist) <- marker_Muller1$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Muller-LIM.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Muller-LIM')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(2,3,7) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Muller-LIM")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


Muller_EA <- marker_Muller2$gene
Muller_EA <- as.data.frame(Muller_EA)
entrez<- bitr(Muller_EA$Muller_EA ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Muller2$avg_log2FC
names(genelist) <- marker_Muller2$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Muller-EA.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Muller-EA')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(2,3,7) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Muller-EA")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


#Neuron#

Neuron_LIM <- marker_Neuron1$gene
Neuron_LIM <- as.data.frame(Neuron_LIM)
entrez<- bitr(Neuron_LIM$Neuron_LIM ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Neuron1$avg_log2FC
names(genelist) <- marker_Neuron1$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Neuron-LIM.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Neuron-LIM')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(8,19,11) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Neuron-LIM")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


Neuron_EA <- marker_Neuron2$gene
Neuron_EA <- as.data.frame(Neuron_EA)
entrez<- bitr(Neuron_EA$Neuron_EA ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Neuron2$avg_log2FC
names(genelist) <- marker_Neuron2$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Neuron-EA.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Neuron-EA')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(9,10,16) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Neuron-EA")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


#Rod#

Rod_LIM <- marker_Rod1$gene
Rod_LIM <- as.data.frame(Rod_LIM)
entrez<- bitr(Rod_LIM$Rod_LIM ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Rod1$avg_log2FC
names(genelist) <- marker_Rod1$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Rod-LIM.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Rod-LIM')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(2,13,18) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Rod-LIM")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P

Rod_EA <- marker_Rod2$gene
Rod_EA <- as.data.frame(Rod_EA)
entrez<- bitr(Rod_EA$Rod_EA ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Rod2$avg_log2FC
names(genelist) <- marker_Rod2$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Rod-EA.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Rod-EA')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(1,4,5) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Rod-EA")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P

#RPE#

Rpe_LIM <- marker_Rpe1$gene
Rpe_LIM <- as.data.frame(Rpe_LIM)
entrez<- bitr(Rpe_LIM$Rpe_LIM ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Rpe1$avg_log2FC
names(genelist) <- marker_Rpe1$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Rpe-LIM.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Rpe-LIM')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(1,3,4) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Rpe-LIM")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P


Rpe_EA <- marker_Rpe2$gene
Rpe_EA <- as.data.frame(Rpe_EA)
entrez<- bitr(Rpe_EA$Rpe_EA ,
              fromType= "SYMBOL" ,
              toType= "ENTREZID" ,
              OrgDb= "org.Hs.eg.db" )
genelist <- marker_Rpe2$avg_log2FC
names(genelist) <- marker_Rpe2$gene
genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
genelist<- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges<- gseKEGG(
  geneList= genelist ,
  organism= "hsa" ,
  minGSSize= 10 ,
  maxGSSize= 500 ,
  pvalueCutoff= 0.99 ,
  pAdjustMethod= "BH" ,
  verbose= FALSE ,
  eps= 0
)
KEGG_ges<- setReadable(KEGG_ges,
                       OrgDb= org.Hs.eg.db,
                       keyType= "ENTREZID")
KEGG_ges_result<- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('GSEA-Rpe-EA.csv'))


ridgeplot(KEGG_ges,
          showCategory= 15,
          fill= "pvalue",
          decreasing= T) +
  theme_minimal() +
  labs(title = 'GSEA Rpe-EA')

P <- gseaplot2(KEGG_ges ,
               geneSetID = c(2,6,8) ,
               color = c("#a2d2e7" , "#ffc17f" , "#ff9d9f") ,
               pvalue_table = F ,
               title = "Rpe-EA")
P[[1]] <- P[[1]] + geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)

P
P[[3]] <- P[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
P
