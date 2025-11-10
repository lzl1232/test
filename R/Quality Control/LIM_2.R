# install.packages("Seurat")
# 安装seurat   dplyr在tidyverse包里

library(dplyr)
library(Seurat)
library(ggplot2)


LIM_2 <- Read10X(data.dir = "single\\LIM-2\\" , gene.column = 1)

LIM_2_SEU <- CreateSeuratObject(counts = LIM_2 ,
                                project = "LIM_sample_2" , 
                                min.cells = 3 ,
                                min.features = 200)
LIM_2_SEU[["percent.mt"]] <- PercentageFeatureSet(LIM_2_SEU , pattern = "^MT-" )
#图1--C_F_P
pic1 <- VlnPlot(LIM_2_SEU ,features = c("nCount_RNA" , "nFeature_RNA" , "percent.mt" , ncol = 3))
ggsave("L2_pic1.pdf" , pic1 , width = 30 , height = 10 ,dpi = 300)
#图2--C_F
pic2 <- FeatureScatter(LIM_2_SEU , feature1 = "nCount_RNA" , feature2 = "nFeature_RNA")
ggsave("L2_pic2.pdf" , pic2 , width = 10 , height = 10 ,dpi = 300)
#图3--C_P
pic3 <- FeatureScatter(LIM_2_SEU , feature1 = "nCount_RNA" , feature2 = "percent.mt")
ggsave("L2_pic3.pdf" , pic3 , width = 10 , height = 10 ,dpi = 300)
#质控
LIM_2_SEU <- subset(LIM_2_SEU , subset = nFeature_RNA > 200 & nFeature_RNA < 7000)

#归一水平化
LIM_2_SEU <- NormalizeData(LIM_2_SEU , normalization.method = "LogNormalize" , scale.factor = 10000)
#寻找显著表达基因
LIM_2_SEU <- FindVariableFeatures(LIM_2_SEU , selection.method = "vst" , nfeatures = 2000)
#改为基因名
all.genes <- rownames(LIM_2_SEU)
LIM_2_SEU <- ScaleData(LIM_2_SEU , features = all.genes)
#PCA降维
LIM_2_SEU <- RunPCA(LIM_2_SEU , features = VariableFeatures(object = LIM_2_SEU))
top10 <- head(VariableFeatures(LIM_2_SEU) , 10)
#图4--高表达基因展示
plot1 <- VariableFeaturePlot(LIM_2_SEU)+
  ggtitle("LIM_sample_2")
pic4 <- LabelPoints(plot = plot1 , points = top10 , repel = T)
ggsave("L2_pic4.pdf" , pic4 , width = 10 , height = 10 ,dpi = 300)
#图5--分组差异
pic5 <- ElbowPlot(LIM_2_SEU , ndims = 50)
ggsave("L2_pic5.pdf" , pic5 , width = 10 , height = 10 ,dpi = 300)
#找邻居
LIM_2_SEU <- FindNeighbors(LIM_2_SEU , dims = 1:30 )
LIM_2_SEU <- FindClusters(LIM_2_SEU , resolution = 0.5)
LIM_2_SEU <- RunUMAP(LIM_2_SEU , dims = 1:30)
#图6--umap图
pic6 <- DimPlot(LIM_2_SEU , reduction = "umap")
ggsave("L2_pic6.pdf" , pic6 , width = 10 , height = 10 ,dpi = 300)
LIM_2_SEU <- RunTSNE(LIM_2_SEU , dims = 1:30)
#图7--tsne图
pic7 <- DimPlot(LIM_2_SEU , reduction = "tsne")
ggsave("L2_pic7.pdf" , pic7 , width = 10 , height = 10 ,dpi = 300)

save(LIM_2_SEU , file = "LIM_2.Rdata")


