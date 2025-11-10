grey <- subset(modules, modules$color == 'green')
brown <- subset(modules, modules$color == 'brown')
blue <- subset(modules, modules$color == 'blue')
turquoise <- subset(modules, modules$color == 'turquoise')
yellow <- subset(modules, modules$color == 'yellow')

brown2 <- subset(hub_df, hub_df$module == 'brown')
write.csv(brown2, file = 'brown_gene.csv')
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


#### brown ####
ID<-bitr(brown$gene_name,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_brown<-enrichGO(gene_diff,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 pvalueCutoff = 1 ,
                 maxGSSize = 5000,
                 minGSSize = 5,
                 ont = "all")
EGO_brown<-setReadable(EGO_brown,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_brown, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_brown<-enrichKEGG(ID$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg")
barplot(kEGG_brown@result,x = "GeneRatio" , showCategory =20)
kEGG_brown<-setReadable(kEGG_brown,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

write.csv(kEGG_brown@result, file = 'brown_GO.csv')
write.csv(kEGG_brown@result, file = 'brown_KEGG.csv')

#### turquoise ####
ID<-bitr(turquoise$gene_name,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_turquoise<-enrichGO(gene_diff,
                    OrgDb = "org.Hs.eg.db",
                    keyType = "ENTREZID",
                    pvalueCutoff = 1 ,
                    maxGSSize = 5000,
                    minGSSize = 5,
                    ont = "all")
EGO_turquoise<-setReadable(EGO_turquoise,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_turquoise, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_turquoise<-enrichKEGG(ID$ENTREZID,
                       organism = "hsa",
                       keyType = "kegg")
barplot(kEGG_turquoise@result,x = "GeneRatio" , showCategory =20)
kEGG_turquoise<-setReadable(kEGG_turquoise,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

write.csv(EGO_turquoise@result, file = 'turquoise_GO.csv')
write.csv(kEGG_turquoise@result, file = 'turquoisen_KEGG.csv')


#### blue ####
ID<-bitr(blue$gene_name,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_blue<-enrichGO(gene_diff,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENTREZID",
                        pvalueCutoff = 1 ,
                        maxGSSize = 5000,
                        minGSSize = 5,
                        ont = "all")
EGO_blue<-setReadable(EGO_blue,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_blue, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_blue<-enrichKEGG(ID$ENTREZID,
                           organism = "hsa",
                           keyType = "kegg")
barplot(kEGG_blue@result,x = "GeneRatio" , showCategory =20)
kEGG_blue<-setReadable(kEGG_blue,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

write.csv(EGO_blue@result, file = 'blue_GO.csv')
write.csv(kEGG_blue@result, file = 'blue_KEGG.csv')

#### yellow ####
ID<-bitr(yellow$gene_name,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_yellow<-enrichGO(gene_diff,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   pvalueCutoff = 1 ,
                   maxGSSize = 5000,
                   minGSSize = 5,
                   ont = "all")
EGO_yellow<-setReadable(EGO_yellow,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_yellow, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_yellow<-enrichKEGG(ID$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg")
barplot(kEGG_yellow@result,x = "GeneRatio" , showCategory =20)
kEGG_yellow<-setReadable(kEGG_yellow,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

write.csv(EGO_yellow@result, file = 'yellow_GO.csv')
write.csv(kEGG_yellow@result, file = 'yellow_KEGG.csv')


#### yellow ####
ID<-bitr(grey$gene_name,"SYMBOL","ENTREZID","org.Hs.eg.db")
duplicated_col <-  duplicated(ID)
ID <- ID[!duplicated_col , ]
gene_diff<-ID$ENTREZID
EGO_yellow<-enrichGO(gene_diff,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     pvalueCutoff = 1 ,
                     maxGSSize = 5000,
                     minGSSize = 5,
                     ont = "all")
EGO_yellow<-setReadable(EGO_yellow,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(EGO_yellow, x = "GeneRatio", color = "p.adjust",
        showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
kEGG_yellow<-enrichKEGG(ID$ENTREZID,
                        organism = "hsa",
                        keyType = "kegg")
barplot(kEGG_yellow@result,x = "GeneRatio" , showCategory =20)
kEGG_yellow<-setReadable(kEGG_yellow,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

write.csv(EGO_yellow@result, file = 'yellow_GO.csv')
write.csv(kEGG_yellow@result, file = 'yellow_KEGG.csv')


library(MASS)
data("birthwt")
str(birthwt)
write.csv(birthwt, file = 'birthwt.csv')
