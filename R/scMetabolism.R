#install.packages(c("devtools", "data.table", "wesanderson", "Seurat",  "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
#devtools::install_github("YosefLab/VISION@v2.1.0")
#getOption('timeout')
#options(timeout=10000)
#devtools::install_github("wu-yc/scMetabolism")
library(scMetabolism)
library(remotes)
library(loe)
library(ggplot2)
library(rsvd)
library(pheatmap)
library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)
library(reshape2)

#install_version("loe" , version = "1.1" , repos = "http://cran.us.r-project.org")###loe被从

countexp.Seurat$stim <- recode(countexp.Seurat$orig.ident ,
                   "NC1" = "NC" ,
                   "NC2" = "NC" ,
                   "LIM1" = "LIM" ,
                   "LIM2" = "LIM" ,
                   "EA1" = "EA" ,
                   "EA2" = "EA")
save(EYE, file = './EYE_zuizhong.Rdata')



countexp.Seurat <- sc.metabolism.Seurat(obj = EYE , method = "AUCell" , imputation = F , ncores = 2 , metabolism.type = "KEGG")
#sce2@meta.data$CB <- rownames(sce2@meta.data)

library(Seurat)


EYE@assays$RNA@counts <- exprMat

umap.loc <- countexp.Seurat@reductions$umap@cell.embeddings
colnames(umap.loc) <- c("UMAP_1" , "UMAP_2")
countexp.Seurat@reductions$umap@cell.embeddings <- umap.loc


pathways <- countexp.Seurat@assays$METABOLISM$score

head(rownames(pathways))
DimPlot.metabolism(obj = countexp.Seurat ,
                   pathway = "Glycolysis / Gluconeogenesis" ,
                   dimention.reduction.type = "umap" ,
                   dimention.reduction.run = F , size = 1 )

input.pathway = rownames(countexp.Seurat@assays$METABOLISM$score)

input.pathway = c("Alanine, aspartate and glutamate metabolism" ,
                  "Arginine and proline metabolism" ,
                  "Arginine biosynthesis" ,
                  "Cysteine and methionine metabolism" ,
                 # "Glycosaminoglycan biosynthesis - heparan sulfate / heparin" ,
                  "Oxidative phosphorylation" ,
                  "Propanoate metabolism" ,
                  "Phenylalanine, tyrosine and tryptophan biosynthesis" ,
                  "Purine metabolism" ,
                  "Propanoate metabolism" ,
                  "Pyrimidine metabolism")



input.pathway = c("Citrate cycle (TCA cycle)" ,
                  "Fatty acid elongation" ,
                  "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis" ,
                  "Glycerophospholipid metabolism" ,
                  "Glyoxylate and dicarboxylate metabolism" ,
                  "Pyrimidine metabolism" ,
                  "Terpenoid backbone biosynthesis")


p1 <- DotPlot.metabolism(obj = countexp.Seurat ,
                   pathway = input.pathway ,
                   phenotype = "stim" , #这个参数需按需修改
                   norm = "y") +
  scale_color_gradient(high = "#FF8831",low = "#50AA4B") +
  theme_bw() +
  labs(x="", y = " ",title = " ") +
  theme(
    axis.text.y = element_text(size = 10,
                               face = 'bold',
                               colour = 'black'),
    axis.text.x = element_text(size = 10,
                               face = 'bold',
                               colour = 'black'),
  )

p1


#remove.packages("scMetabolism")

#devtools::install_local("E:\\R\\R.project\\eye汇总\\scMetabolism.zip")






#
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
#将score中barcode的点转为下划线
score_change <- score %>%
  select_all(~str_replace_all(., "\\.", "-"))  #基因ID不规范会报错,下划线替换-
#确定细胞barcode椅子
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
#[1] TRUE
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change) )




###
df <- countexp.Seurat@meta.data
#19列开始是代谢通路的得分，按照celltype计算均值
avg_df = aggregate(df[,10:ncol(df)],
                   list(df$stim),
                   mean)


avg_df <- avg_df %>%
  #select(1:20) %>% #展示前20个
  column_to_rownames("Group.1")
avg_df[1:4,1:4]

exp <- apply(avg_df, 2, scale)
rownames(exp) <- rownames(avg_df)
# 组件


h_state <- Heatmap(t(exp),
                   column_title = "state_gsva",
                   col = colorRampPalette(c('#1A5592','white',"#FF8831"))(100),
                   name= "gsva ",
                   show_row_names = TRUE,
                   show_column_names = TRUE)
h_state
class( h_state )


exp <- t(exp)

exp <- as.data.frame(exp) %>%   mutate(mtxars=row.names(.)) %>% melt()



p1<-ggplot(exp,aes(x=mtxars,y=variable,fill=value)) #热图绘制


p2 <- p1+geom_raster()+scale_fill_gradient2(low="#50AA4B", high="#FF8831", mid="white")
p2
p3 <- p2 +
  theme_bw()+
  theme(axis.text.x =element_text(angle =90,hjust= 1,vjust = 0.5))+
  xlab(NULL) + ylab(NULL) +
  theme(
    axis.text.y = element_text(size = 10,
                               face = 'bold',
                               colour = 'black'),
    axis.text.x = element_text(size = 10,
                               face = 'bold',
                               colour = 'black'),
    legend.key.size = unit(10, "pt")
  )


p3   #1500  600

library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("EA" , "LIM") , c("LIM" , "NC"))
BoxPlot.metabolism(obj = countexp.Seurat,
                   pathway = 'Glycolysis / Gluconeogenesis',
                   phenotype = "stim",
                   ncol = 2) +
  scale_fill_nejm() +
  stat_compare_means(method = "t.test" , hide.ns = F ,
                     comparisons = my_comparisons ,
                     label = "p.signif" ,
                     bracket.size = 0.8 ,
                     label.y = c(0.45, 0.5),
                     tip.length = 0 ,
                     size = 6) +
  ylim(0, 0.6) +
  theme(axis.title.x = element_blank() ,
        axis.text.x = element_text(colour = "black" , face = "bold" , size = 12) ,
        axis.text.y = element_text(colour = "black" , face = "bold") ,
        axis.title.y = element_text(colour = "black" , face = "bold" , size = 15) ,
        panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank() ,
        panel.border = element_rect(color = "black" , size = 1.2 , linetype = "solid") ,
        panel.spacing = unit(0.12 , "cm") ,
        plot.title = element_text(hjust = 0.5 , face = "bold.italic") ,
        legend.position = "none")

save(countexp.Seurat, file = 'Cone_代谢.Rdata')
