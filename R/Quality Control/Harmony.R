library(Seurat)
library(usethis)
library(devtools)
# devtools::install_github("satijalab/seurat-data")
library(SeuratData)
library(patchwork)

load("NC_1.Rdata")
load("NC_2.Rdata")
load("EA_1.Rdata")
load("EA_2.Rdata")
load("LIM_1.Rdata")
load("LIM_2.Rdata")


NC_1_COUNTS <- GetAssayData(object = NC_1_SEU , layer = "counts")
NC_2_COUNTS <- GetAssayData(object = NC_2_SEU , layer = "counts")
EA_1_COUNTS <- GetAssayData(object = EA_1_SEU , layer = "counts")
EA_2_COUNTS <- GetAssayData(object = EA_2_SEU , layer = "counts")
LIM_1_COUNTS <- GetAssayData(object = LIM_1_SEU , layer = "counts")
LIM_2_COUNTS <- GetAssayData(object = LIM_2_SEU , layer = "counts")


rm(NC_1_SEU , NC_2_SEU , EA_1_SEU , EA_2_SEU , LIM_1_SEU , LIM_2_SEU)



Counts_jj <- intersect(row.names(NC_1_COUNTS) , row.names(EA_1_COUNTS))
Counts_jj <- intersect(Counts_jj , row.names(EA_1_COUNTS))
Counts_jj <- intersect(Counts_jj , row.names(EA_2_COUNTS))
Counts_jj <- intersect(Counts_jj , row.names(LIM_1_COUNTS))
Counts_jj <- intersect(Counts_jj , row.names(LIM_2_COUNTS))

NC_1_jjp <- match(Counts_jj , row.names(NC_1_COUNTS))
NC_2_jjp <- match(Counts_jj , row.names(NC_2_COUNTS))
EA_1_jjp <- match(Counts_jj , row.names(EA_1_COUNTS))
EA_2_jjp <- match(Counts_jj , row.names(EA_2_COUNTS))
LIM_1_jjp <- match(Counts_jj , row.names(LIM_1_COUNTS))
LIM_2_jjp <- match(Counts_jj , row.names(LIM_2_COUNTS))

NC_1_jjp <- NC_1_COUNTS[NC_1_jjp ,]
NC_2_jjp <- NC_2_COUNTS[NC_2_jjp ,]
EA_1_jjp <- NC_1_COUNTS[EA_1_jjp ,]
EA_2_jjp <- EA_2_COUNTS[EA_2_jjp ,]
LIM_1_jjp <- NC_1_COUNTS[LIM_1_jjp ,]
LIM_2_jjp <- LIM_2_COUNTS[LIM_2_jjp ,]


colnames(NC_1_jjp) <- paste("NC1_cell" , C(1 : ncol(NC_1_jjp)) , sep = "_")
colnames(NC_2_jjp) <- paste("NC2_cell" , C(1 : ncol(NC_2_jjp)) , sep = "_")
colnames(EA_1_jjp) <- paste("EA1_cell" , C(1 : ncol(EA_1_jjp)) , sep = "_")
colnames(EA_2_jjp) <- paste("EA2_cell" , C(1 : ncol(EA_2_jjp)) , sep = "_")
colnames(LIM_1_jjp) <- paste("LIM1_cell" , C(1 : ncol(LIM_1_jjp)) , sep = "_")
colnames(LIM_2_jjp) <- paste("LIM2_cell" , C(1 : ncol(LIM_2_jjp)) , sep = "_")


count_matrix <- list(NC_1_jjp,
                     NC_2_jjp,
                     EA_1_jjp,
                     EA_2_jjp,
                     LIM_1_jjp,
                     LIM_2_jjp,
                     )


count_matrix <- do.call("cbind" , count_matrix)
table(duplicated(colnames(count_matrix)))

g_p <- data.frame(sampleID = c(colnames(NC_1_jjp) ,
                               colnames(NC_2_jjp) , 
                               colnames(EA_1_jjp) ,
                               colnames(EA_2_jjp) , 
                               colnames(LIM_1_jjp),
                               colnames(LIM_2_jjp)),
                  group = c(rep("NC_sample1" , length(colnames(NC_1_jjp))) ,
                            rep("NC_sample2" , length(colnames(NC_2_jjp))) ,
                            rep("EA_sample1" , length(colnames(EA_1_jjp))) ,
                            rep("EA_sample2" , length(colnames(EA_2_jjp))) , 
                            rep("LIM_sample1" , length(colnames(LIM_1_jjp))) ,
                            rep("LIM_sample2" , length(colnames(LIM_2_jjp)))))


pbmc <- CreateSeuratObject(counts = as.matrix(count_matrix) , project = "Myopic" , min.cells = 0 , min.features = 0)

rm(NC_1_jjp , NC_2_jjp , EA_1_jjp , EA_2_jjp , LIM_1_jjp , LIM_2_jjp)
rm(NC_1_COUNTS , NC_2_COUNTS , EA_1_COUNTS , EA_2_COUNTS , LIM_1_COUNTS , LIM_2_COUNTS)


pp <- match(row.names(pbmc@meta.data) , g_p$sampleID)
pbmc@meta.data$group <- g_p$group[pp]


pbmc[["RNA"]] <- split(pbmc[["RNA"]] , f = pbmc$group)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc , dims = 1:30 , reduction = "pca")
pbmc <- FindClusters(pbmc , resolution = 1)
pbmc <- RunUMAP(pbmc , dims = 1:30 , reduction = "pca" , reduction.name = "umap.unintegrated")
plot2<- DimPlot(pbmc , reduction = "umap.unintegrated" , group.by = c("group" , "seurat_clusters"))
plot2

install.packages("harmony")
library(harmony)
scobj <- RunHarmony(pbmc,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:30,reduction.name = "umap_harmony")
plot1<- DimPlot (scobj , reduction = "umap_harmony" , group.by = c("group" , "seurat_clusters"))
plot1
scobj <- IntegrateLayers(object = scobj , method = CCAIntegration , orig.reduction = "harmony" , new.reduction = "integrated.cca" ,
                        verbose = F)
scobj[["RNA"]] <- JoinLayers(scobj[["RNA"]])
scobj <- FindNeighbors(scobj , reduction = "integrated.cca" , dims = 1:30)
scobj <- FindClusters(scobj , resolution = 1)
scobj <- RunUMAP(scobj , dims = 1:30 , reduction = "integrated.cca")
plot1<- DimPlot (scobj , reduction = "umap" , group.by = c("group" , "seurat_clusters"))
plot1


SUMMARY_seurat <- scobj

save(SUMMARY_seurat , file = "sum.Rdata")

