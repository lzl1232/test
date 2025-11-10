# install.packages("Seurat")
# 安装seurat   dplyr在tidyverse包里

library(dplyr)
library(Seurat)
library(ggplot2)


NC_2 <- Read10X(data.dir = "single\\NC-2\\" , gene.column = 1)

NC_2_SEU <- CreateSeuratObject(counts = NC_2 ,
                               project = "NC_sample_2" , 
                               min.cells = 3 ,
                               min.features = 200)
NC_2_SEU[["percent.mt"]] <- PercentageFeatureSet(NC_2_SEU , pattern = "^MT-" )
#图1--C_F_P
pic1 <- VlnPlot(NC_2_SEU ,features = c("nCount_RNA" , "nFeature_RNA" , "percent.mt" , ncol = 3))
ggsave("N2_pic1.pdf" , pic1 , width = 30 , height = 10 ,dpi = 300)
#图2--C_F
pic2 <- FeatureScatter(NC_2_SEU , feature1 = "nCount_RNA" , feature2 = "nFeature_RNA")
ggsave("N2_pic2.pdf" , pic2 , width = 10 , height = 10 ,dpi = 300)
#图3--C_P
pic3 <- FeatureScatter(NC_2_SEU , feature1 = "nCount_RNA" , feature2 = "percent.mt")
ggsave("N2_pic3.pdf" , pic3 , width = 10 , height = 10 ,dpi = 300)
#质控
NC_2_SEU <- subset(NC_2_SEU , subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

#归一水平化
NC_2_SEU <- NormalizeData(NC_2_SEU , normalization.method = "LogNormalize" , scale.factor = 10000)
#寻找显著表达基因
NC_2_SEU <- FindVariableFeatures(NC_2_SEU , selection.method = "vst" , nfeatures = 2000)
#改为基因名
all.genes <- rownames(NC_2_SEU)
NC_2_SEU <- ScaleData(NC_2_SEU , features = all.genes)
#PCA降维
NC_2_SEU <- RunPCA(NC_2_SEU , features = VariableFeatures(object = NC_2_SEU))
top10 <- head(VariableFeatures(NC_2_SEU) , 10)
#图4--高表达基因展示
plot1 <- VariableFeaturePlot(NC_2_SEU)+
  ggtitle("NC_sample_2")
pic4 <- LabelPoints(plot = plot1 , points = top10 , repel = T)
ggsave("N2_pic4.pdf" , pic4 , width = 10 , height = 10 ,dpi = 300)
#图5--分组差异
pic5 <- ElbowPlot(NC_2_SEU , ndims = 50)
ggsave("N2_pic5.pdf" , pic5 , width = 10 , height = 10 ,dpi = 300)
#找邻居
NC_2_SEU <- FindNeighbors(NC_2_SEU , dims = 1:30 )
NC_2_SEU <- FindClusters(NC_2_SEU , resolution = 0.5)
NC_2_SEU <- RunUMAP(NC_2_SEU , dims = 1:30)
#图6--umap图
pic6 <- DimPlot(NC_2_SEU , reduction = "umap")
pic6
ggsave("N2_pic6.pdf" , pic6 , width = 10 , height = 10 ,dpi = 300)
NC_2_SEU <- RunTSNE(NC_2_SEU , dims = 1:30)
#图7--tsne图
pic7 <- DimPlot(NC_2_SEU , reduction = "tsne")
ggsave("N2_pic7.pdf" , pic7 , width = 10 , height = 10 ,dpi = 300)

save(NC_2_SEU , file = "NC_2.Rdata")
