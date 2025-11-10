# install.packages("Seurat")
# 安装seurat   dplyr在tidyverse包里

library(dplyr)
library(Seurat)
library(ggplot2)


LIM_1 <- Read10X(data.dir = "single\\LIM-1\\" , gene.column = 1)

LIM_1_SEU <- CreateSeuratObject(counts = LIM_1 ,
                               project = "LIM_sample_1" , 
                               min.cells = 3 ,
                               min.features = 200)
LIM_1_SEU[["percent.mt"]] <- PercentageFeatureSet(LIM_1_SEU , pattern = "^MT-" )
#图1--C_F_P
pic1 <- VlnPlot(LIM_1_SEU ,features = c("nCount_RNA" , "nFeature_RNA" , "percent.mt" , ncol = 3))
ggsave("L1_pic1.pdf" , pic1 , width = 30 , height = 10 ,dpi = 300)
#图2--C_F
pic2 <- FeatureScatter(LIM_1_SEU , feature1 = "nCount_RNA" , feature2 = "nFeature_RNA")
pic2
ggsave("L1_pic2.pdf" , pic2 , width = 10 , height = 10 ,dpi = 300)
#图3--C_P
pic3 <- FeatureScatter(LIM_1_SEU , feature1 = "nCount_RNA" , feature2 = "percent.mt")
ggsave("L1_pic3.pdf" , pic3 , width = 10 , height = 10 ,dpi = 300)
#质控
LIM_1_SEU <- subset(LIM_1_SEU , subset = nFeature_RNA > 200 & nFeature_RNA < 7000)

#归一水平化
LIM_1_SEU <- NormalizeData(LIM_1_SEU , normalization.method = "LogNormalize" , scale.factor = 10000)
#寻找显著表达基因
LIM_1_SEU <- FindVariableFeatures(LIM_1_SEU , selection.method = "vst" , nfeatures = 2000)
#改为基因名
all.genes <- rownames(LIM_1_SEU)
LIM_1_SEU <- ScaleData(LIM_1_SEU , features = all.genes)
#PCA降维
LIM_1_SEU <- RunPCA(LIM_1_SEU , features = VariableFeatures(object = LIM_1_SEU))
top10 <- head(VariableFeatures(LIM_1_SEU) , 10)
#图4--高表达基因展示
plot1 <- VariableFeaturePlot(LIM_1_SEU)+
  ggtitle("LIM_sample_1")
pic4 <- LabelPoints(plot = plot1 , points = top10 , repel = T)
ggsave("L1_pic4.pdf" , pic4 , width = 10 , height = 10 ,dpi = 300)
#图5--分组差异
pic5 <- ElbowPlot(LIM_1_SEU , ndims = 50)
ggsave("L1_pic5.pdf" , pic5 , width = 10 , height = 10 ,dpi = 300)
#找邻居
LIM_1_SEU <- FindNeighbors(LIM_1_SEU , dims = 1:30 )
LIM_1_SEU <- FindClusters(LIM_1_SEU , resolution = 0.5)


LIM_1_SEU <- RunUMAP(LIM_1_SEU , dims = 1:30)
#图6--umap图
DimPlot(LIM_1_SEU , reduction = "umap")
ggsave("L1_pic6.pdf" , pic6 , width = 10 , height = 10 ,dpi = 300)
LIM_1_SEU <- RunTSNE(LIM_1_SEU , dims = 1:30)
#图7--tsne图
DimPlot(LIM_1_SEU , reduction = "tsne")
ggsave("L1_pic7.pdf" , pic7 , width = 10 , height = 10 ,dpi = 300)

save(LIM_1_SEU , file = "LIM_1.Rdata")


