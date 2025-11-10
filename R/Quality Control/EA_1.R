# install.packages("Seurat")
# 安装seurat   dplyr在tidyverse包里

library(dplyr)
library(Seurat)
library(ggplot2)


EA_1 <- Read10X(data.dir = "single\\EA-1\\" , gene.column = 1)

EA_1_SEU <- CreateSeuratObject(counts = EA_1 ,
                               project = "EA_sample_1" , 
                               min.cells = 3 ,
                               min.features = 200)
EA_1_SEU[["percent.mt"]] <- PercentageFeatureSet(EA_1_SEU , pattern = "^MT-" )
#图1--C_F_P
pic1 <- VlnPlot(EA_1_SEU ,features = c("nCount_RNA" , "nFeature_RNA" , "percent.mt" , ncol = 3))
ggsave("E1_pic1.pdf" , pic1 , width = 30 , height = 10 ,dpi = 300)
#图2--C_F
pic2 <- FeatureScatter(EA_1_SEU , feature1 = "nCount_RNA" , feature2 = "nFeature_RNA")
ggsave("E1_pic2.pdf" , pic2 , width = 10 , height = 10 ,dpi = 300)
#图3--C_P
pic3 <- FeatureScatter(EA_1_SEU , feature1 = "nCount_RNA" , feature2 = "percent.mt")
ggsave("E1_pic3.pdf" , pic3 , width = 10 , height = 10 ,dpi = 300)
#质控
EA_1_SEU <- subset(EA_1_SEU , subset = nFeature_RNA > 200 & nFeature_RNA < 7000)

#归一水平化
EA_1_SEU <- NormalizeData(EA_1_SEU , normalization.method = "LogNormalize" , scale.factor = 10000)
#寻找显著表达基因
EA_1_SEU <- FindVariableFeatures(EA_1_SEU , selection.method = "vst" , nfeatures = 2000)
#改为基因名
all.genes <- rownames(EA_1_SEU)
EA_1_SEU <- ScaleData(EA_1_SEU , features = all.genes)
#PCA降维
EA_1_SEU <- RunPCA(EA_1_SEU , features = VariableFeatures(object = EA_1_SEU))
top10 <- head(VariableFeatures(EA_1_SEU) , 10)
#图4--高表达基因展示
plot1 <- VariableFeaturePlot(EA_1_SEU)+
  ggtitle("EA_sample_1")
pic4 <- LabelPoints(plot = plot1 , points = top10 , repel = T)
ggsave("E1_pic4.pdf" , pic4 , width = 10 , height = 10 ,dpi = 300)
#图5--分组差异
pic5 <- ElbowPlot(EA_1_SEU , ndims = 50)
ggsave("E1_pic5.pdf" , pic5 , width = 10 , height = 10 ,dpi = 300)
#找邻居
EA_1_SEU <- FindNeighbors(EA_1_SEU , dims = 1:30 )
EA_1_SEU <- FindClusters(EA_1_SEU , resolution = 0.5)


EA_1_SEU <- RunUMAP(EA_1_SEU , dims = 1:30)
#图6--umap图
pic6 <- DimPlot(EA_1_SEU , reduction = "umap")
ggsave("E1_pic6.pdf" , pic6 , width = 10 , height = 10 ,dpi = 300)
EA_1_SEU <- RunTSNE(EA_1_SEU , dims = 1:30)
#图7--tsne图
pic7 <- DimPlot(EA_1_SEU , reduction = "tsne")
ggsave("E1_pic7.pdf" , pic7 , width = 10 , height = 10 ,dpi = 300)

save(EA_1_SEU , file = "EA_1.Rdata")

