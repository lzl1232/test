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

load("sum.Rdata")
pbmc <- FindClusters(object = pbmc , resolution = c(seq(from = .1 , to = 1.6 , by = .2)))
clustree(pbmc@meta.data , prefix = "RNA_snn_res.")
DimPlot(pbmc , reduction = "umap" , label = T , group.by = "RNA_snn_res.0.5")
DimPlot(pbmc , reduction = "tSNE" , label = T)
ALLmarker <- FindAllMarkers(pbmc , logfc.threshold = 0.25 , only.pos = T)


ALLmarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() ->top10
DoHeatmap(pbmc , features = top10$gene)

pbmc$celltype <- recode(pbmc@meta.data$RNA_snn_res.0.5 , 
                        "0" = "0" , 
                        "6" = "3" , 
                        "8" = "3" , 
                        "19" = "0" , 
                        "1" = "1" , 
                        "4" = "1" , 
                        "21" = "1" , 
                        "2" = "2" , 
                        "3" = "2" , 
                        "7" = "2" , 
                        "5" = "5" , 
                        "9" = "9" , 
                        "10" = "10" , 
                        "11" = "11" , 
                        "12" = "12" , 
                        "13" = "13" , 
                        "14" = "14" , 
                        "15" = "15" , 
                        "16" = "16" , 
                        "17" = "17" , 
                        "18" = "18" , 
                        "20" = "20" , 
                        "22" = "22" , 
                        "23" = "23" , 
                        "24" = "24")
DimPlot(pbmc , reduction = "umap" , group.by = "celltype" , label = T)

Idents(pbmc) <- "celltype"
ALLLmarker <- FindAllMarkers(pbmc , logfc.threshold = 0.25 , only.pos = T)
ALLLmarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() ->top10
DoHeatmap(pbmc , features = top10$gene)

DotPlot(pbmc , features = unique(top10$gene) , assay = "RNA" ) + 
  coord_flip()

pbmc$zhushi <- recode(pbmc$celltype , 
                      "1" = "Rod" , 
                      "2" = "Cone" , 
                      "20" = "Microglia" , 
                      "22" = "Rpe" , 
                      "9" = "Amacrine" , 
                      "11" = "Bipolar" , 
                      "15" = "Bipolar" ,
                      "16" = "Bipolar" , 
                      "18" = "Bipolar" , 
                      "5" = "Bipolar" , 
                      "12" = "Ganglion" , 
                      "13" = "Ganglion" , 
                      "17" = "Ganglion" , 
                      "0" = "Muller" , 
                      "14" = "Neuron" , 
                      "10" = "Ganglion" , 
                      "3" = "Muller")     

FeaturePlot(pbmc , features = "CNTN2" , reduction = "umap")
Idents(pbmc) <- "zhushi"
DimPlot(pbmc , reduction = "umap" , label = T)

allll.marker <- FindAllMarkers(pbmc , only.pos = T , logfc.threshold = 0.25)
allll.marker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() ->top10
DoHeatmap(pbmc , features = top10$gene)
DimPlot(pbmc , label = T)
EYE <- subset(pbmc , idents = c("Rod" , "Cone" , "Microglia" , "Rpe"  , "Amacrine" , 
                                "Bipolar" , "Ganglion" , "Muller" , "Neuron"  ))
EYE@meta.data <- EYE@meta.data[ , -7]
save(EYE , file = "EYE.Rdata")
