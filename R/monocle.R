library(Seurat)
library(tidyverse)
library(devtools)
#install.packages("qlcMatrix")
#devtools::install_github("cole-trapnell-lab/leidenbase")
#devtools::install_github("cole-trapnell-lab/monocle3")
#library(monocle3)
#BiocManager::install("monocle3")
#install.packages("D:/BaiduNetdiskDownload/monocle_2.24.0.tar.gz"  , repos = NULL , type = "source")
library(monocle)
library(grid)
#options(timeout = 600000000)
#devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
library(CytoTRACE2)
library(gridExtra)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(igraph)
library(ggrastr)
library(tidydr)
#install.packages("tidydr")
library(viridis)
library(ggforce)
# install.packages('ggforce')
#install.packages("C:/Users/31101/Downloads/igraph_2.0.3.tar.gz", repos = NULL, type = "source")




result_sce <- cytotrace2(Rod ,
                         is_seurat = TRUE,
                         slot_type = "counts",
                         species = "human",
                         seed = 1234)
head(result_sce , 2)
annotation <- data.frame(phenotype=result_sce@meta.data$seurat_clusters) %>% set_rownames(., colnames(result_sce))
plots <- plotData(cytotrace2_result=result_sce, annotation=annotation, is_seurat=TRUE)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
(p1+p2+p3+p4) + plot_layout(ncol = 2)






expr_matrix = GetAssayData(Rod, slot = "counts")


p_data <- Rod@meta.data
head(p_data)
pd <- new('AnnotatedDataFrame', data = p_data)


f_data <- data.frame(gene_short_name = row.names(Rod),
                     row.names = row.names(Rod))
head(f_data)
fd <- new('AnnotatedDataFrame', data = f_data)


cds_pre <- newCellDataSet(expr_matrix,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size())


##预处理
# Add Size_Factor文库因子
cds_pre <- estimateSizeFactors(cds_pre)
cds_pre$Size_Factor %>% head()

# 计算基因表达量的离散度
cds_pre <- estimateDispersions(cds_pre)
head(dispersionTable(cds_pre))
cds_pre


### 策略1：marker gene by Seurat
Idents(Muller) = "zhushi"
gene_FAM = FindAllMarkers(Muller)
gene_sle = gene_FAM %>%
  dplyr::filter(p_val<0.01) %>%
  pull(gene) %>% unique()

### 策略2：high dispersion gene by monocle
gene_Disp = dispersionTable(cds_pre)
gene_sle = gene_Disp %>%
  dplyr::filter(mean_expression >= 0.1,
                dispersion_empirical >= dispersion_fit) %>%
  pull(gene_id) %>% unique()

### 策略3：variable(high dispersion) gene by Seurat
gene_sle <- VariableFeatures(sce)


cds <- setOrderingFilter(cds_pre, gene_sle)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds , root_state = 3)

?orderCells

table(cds$State, cds$zhushi)
# 按照state分组
plot(cds$State, cds$Pseudotime)
# state分组的cell
plot_cell_trajectory(cds, color_by = "State")
# 树状图
plot_complex_cell_trajectory(cds, color_by = 'seurat_clusters')
# 分组的细胞表示
plot_cell_trajectory(cds, color_by = "seurat_clusters") +
  facet_wrap("~seurat_clusters", nrow = 1)
# 时间表示
plot_cell_trajectory(cds, color_by = "Pseudotime")
# 基因变化
gene_key = c("MYL6")
plot_genes_jitter(cds[gene_key,],
                  grouping = "seurat_clusters", color_by = "seurat_clusters")
# 小提琴图
plot_genes_violin(cds[gene_key,],
                  grouping = "seurat_clusters", color_by = "seurat_clusters")
# 基因变化
plot_genes_in_pseudotime(cds[gene_key,], color_by = "seurat_clusters")
#theme_bw() +
#scale_y_continuous(limits = c(1, 20, 10)) +
#scale_color_manual(values = c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b"))
?plot_genes_in_pseudotime
# 基因表达差异
diff_pseudo <- differentialGeneTest(cds, cores = 1,
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
BEAM_res <- BEAM(cds, branch_point = 3, cores = 4, progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
BEAM_res_2 <- BEAM_res[1:100, ]

plot_genes_branched_heatmap(cds[rownames(BEAM_res_2), ],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)



head(diff_pseudo)
table(diff_pseudo$qval<0.05)

# 差异基因可视化
library(pheatmap)
library(grid)
peu_gene <- diff_pseudo[which(diff_pseudo$qval<0.01) ,]

index <- grep("^ENSCPOG" , peu_gene$gene_short_name)

rownames(peu_gene) = peu_gene$gene_short_name
peu_gene <- peu_gene[-index ,]
peu_gene %>% arrange(qval)  -> peu_gene#按照qval排个序
peu_gene <- peu_gene[1:100,] # 50

plot_pseudotime_heatmap(cds[peu_gene$gene_short_name,],
                        num_clusters = 3,
                        cores = 2,
                        show_rownames = T,return_heatmap =T,
                        hmcols = viridis(256),
                        use_gene_short_name = T )

#clustercolor  = c("#d2981a", "#a53e1f", "#457277", "#8f657d"))


save(cds, file = 'cds.Rdata')
save(diff_pseudo, file = 'diff_pseudo.Rdata')


#####
# 提取点
data_df <- t(reducedDimS(cds)) %>%
  as.data.frame() %>%
  select('Component 1' = 1, 'Component 2' = 2) %>%
  rownames_to_column("Cells") %>%
  mutate(pData(cds)$State,
         pData(cds)$Pseudotime,
         pData(cds)$orig.ident,
         pData(cds)$stim,
         pData(cds)$seurat_clusters )
colnames(data_df) <- c("cells","Component_1","Component_2","State",

                       "Pseudotime","orig.ident",'stim',"celltype")

# 提取线
reduced_dim_coords <- reducedDimK(cds)
ca_space_df <- Matrix::t(reduced_dim_coords) %>%
  as.data.frame() %>%
  select(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))
dp_mst <- minSpanningTree(cds)
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  select(source = "from", target = "to") %>%
  left_join(ca_space_df %>% select(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
  left_join(ca_space_df %>% select(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")
head(edge_df)


ggplot() +
  geom_point_rast(data = data_df, aes(x = Component_1,
                                      y = Component_2,
                                      color =stim
  )) +
  #scale_color_viridis()+
  scale_color_manual(values = c("#80c97f","#a68dc8",
                                "#ffc000")) +
  theme_bw()+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0.8, "lines"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size=15),
  ) +
  geom_segment(aes_string(x="source_prin_graph_dim_1",
                          y="source_prin_graph_dim_2",
                          xend="target_prin_graph_dim_1",
                          yend="target_prin_graph_dim_2"),
               size=0.75, linetype="solid", na.rm=TRUE, data=edge_df) +
  geom_arc_bar(data=subset(Cellratio,State=='7'),stat = "pie",
               aes(x0=1.8,y0=4,r0=0,r=0.8,amount=Freq,fill=stim
               )) +
  geom_arc_bar(data=subset(Cellratio,State=='4'),stat = "pie",
               aes(x0=-10,y0=0,r0=0,r=0.8,amount=Freq,fill=stim )) +
  scale_fill_manual(values = c("#80c97f","#a68dc8",
                               "#ffc000"))
#coord_fixed()



Cellratio <- prop.table(table(data_df$State, data_df$stim), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"stim", 'Freq')

ggplot() +
  geom_point_rast(data = data_df, aes(x = Component_1,
                                      y = Component_2,
                                      color =Pseudotime
  )) +
  scale_color_viridis() +
  #scale_color_gradient(low = "#73a9ad",high = "#f5f0bb") +
  #scale_color_gradientn(values = seq(0,20,5),
  #                      colours = c('#73a9ad','#90c8ac','#c4dfaa','#f5f0bb')) +
  #scale_color_manual(values = c("#80c97f","#a68dc8",
  # "#ffc000")) +
  theme_bw()+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0.8, "lines"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size=15),
  ) +
  geom_segment(aes_string(x="source_prin_graph_dim_1",
                          y="source_prin_graph_dim_2",
                          xend="target_prin_graph_dim_1",
                          yend="target_prin_graph_dim_2"),
               size=0.75, linetype="solid", na.rm=TRUE, data=edge_df)
#  geom_arc_bar(data=subset(Cellratio,State=='6'),stat = "pie",
#               aes(x0=1.8,y0=4,r0=0,r=0.8,amount=Freq,fill=stim
#               )) +
#  geom_arc_bar(data=subset(Cellratio,State=='1'),stat = "pie",
#               aes(x0=-10,y0=0,r0=0,r=0.8,amount=Freq,fill=stim )) +
#  scale_fill_manual(values = c("#80c97f","#a68dc8",
#                               "#ffc000"))
#coord_fixed()

scale_color_gradientn(values = seq(0,20,5),
                      colours = c('#73a9ad','#90c8ac','#c4dfaa','#f5f0bb'))


# 7.04, 4.23



data_df$celltype2 <- recode(data_df$celltype,
                            "0" = "MC_01" ,
                            "1" = "MC_02" ,
                            "2" = "MC_03" ,
                            "4" = "MC_05" ,
                            "3" = "MC_04" )
Cellratio2$celltype <- recode(Cellratio2$celltype,
                              "0" = "MC_01" ,
                              "1" = "MC_02" ,
                              "2" = "MC_03" ,
                              "4" = "MC_05" ,
                              "3" = "MC_04" )



ggplot() +
  geom_point_rast(data = data_df, aes(x = Component_1,
                                      y = Component_2,
                                      color =celltype
  )) +
  #scale_color_viridis()+
  scale_color_manual(values = c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b")) +
  theme_bw()+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0.8, "lines"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size=15),
  ) +
  geom_segment(aes_string(x="source_prin_graph_dim_1",
                          y="source_prin_graph_dim_2",
                          xend="target_prin_graph_dim_1",
                          yend="target_prin_graph_dim_2"),
               size=0.75, linetype="solid", na.rm=TRUE, data=edge_df) +
  geom_arc_bar(data=subset(Cellratio2,State=='6'),stat = "pie",
               aes(x0=1.8,y0=4,r0=0,r=0.8,amount=Freq,fill=celltype
               )) +
  geom_arc_bar(data=subset(Cellratio2,State=='1'),stat = "pie",
               aes(x0=-10,y0=0,r0=0,r=0.8,amount=Freq,fill=celltype )) +
  geom_arc_bar(data=subset(Cellratio2,State=='7'),stat = "pie",
               aes(x0=4,y0=-3,r0=0,r=0.8,amount=Freq,fill=celltype )) +
  scale_fill_manual(values = c("#a2d2e7","#ffc17f","#cf9f88","#b3e19b","#6fb3a3"))

DimPlot(Rod, reduction = 'umap',group.by = 'seurat_clusters')
FeaturePlot(EYE, features = 'SLC24A1', reduction = 'umap')
