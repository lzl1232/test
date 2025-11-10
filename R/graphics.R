library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(scRNAtoolVis)
library(tidyverse)
library(ggVolcano)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(data.table)


load("EYE.Rdata")


#基础绘图
DimPlot(EYE,                           #seurat????(????)
        reduction = 'umap',                  #??ά???ͣ?umap??tsne??pca
        group.by = 'celltype',               #????????
        pt.size = 1,                         #ͼ?е??Ĵ?С
        #split.by = 'seurat_clusters',       #???ֻ?ͼ??????
        label = T)

####更改颜色####
tsnedata <- EYE@reductions$umap@cell.embeddings
tsnedata <- as.data.frame(tsnedata)
tsnedata %>%
  cbind(celltype=EYE@meta.data$celltype) ->tsnedata
tsnedata %>%
  cbind(stim=EYE@meta.data$stim) ->tsnedata
?geom_point
#umap图
celltype_med <- tsnedata %>% group_by(celltype) %>% summarise(umap_1 = median(umap_1) , umap_2 = median(umap_2))
ggplot(tsnedata,aes(umap_1,umap_2,colour=celltype))+
  geom_point(size=0.5)+
  #labs(title = "umap")+
  theme_bw()+
  xlab("")+
  ylab("")+
  #stat_ellipse(type = "norm",linetype=2)+
  scale_color_manual(values =c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b","#704ba3"
                               ))+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 32,
                                  face = "bold",
                                  hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 20,
                                   face = "bold"),
        legend.background = element_blank(),
        legend.key.height = unit(40,"pt"),
        legend.position = "none") +
  geom_label_repel(aes(label = celltype) , fontface = "bold" , data = celltype_med ,
                   point.padding = unit(0.5 , "lines")) +
  geom_segment(aes(x = min(tsnedata$umap_1) , y = min(tsnedata$umap_2) ,
                   xend = min(tsnedata$umap_1) + 2 , yend = min(tsnedata$umap_2) ) ,
               colour = "black" , size = 1 , arrow = arrow(length = unit(0.3 , "cm"))) +
  geom_segment(aes(x = min(tsnedata$umap_1) , y = min(tsnedata$umap_2) ,
                   xend = min(tsnedata$umap_1) , yend = min(tsnedata$umap_2) + 2) ,
               colour = "black" , size = 1 , arrow = arrow(length = unit(0.3 , "cm"))) +
  annotate("text" , x = min(tsnedata$umap_1) + 1 , y = min(tsnedata$umap_2) - 1 , label = "UMAP-1" ,
           color = "black" , size = 3 , fontface = "bold") +
  annotate("text" , x = min(tsnedata$umap_1) - 1 , y = min(tsnedata$umap_2) + 1 , label = "UMAP-2" ,
           color = "black" , size = 3 , fontface = "bold" , angle = 90)
ggsave("umap_stim.jpg",width = 18,height = 18,units = "cm",dpi = 1000)

####细胞注释基因####
gene <- c("VIM" , "CLU" , "PDE6B" , "PDE6A" , "PDE6C" , "GNGT2" , "PCP2" , "TRPM1" , "GAD1" , "TFAP2B" ,'CNTN2', "LAMP5" , "CALB2" , "CD74" , "AIF1" , "RPE65" )
data <- FetchData(EYE , c(gene , "celltype" , "seurat_clusters"))
data$cell <- rownames(data)
data.melt <- reshape2::melt(data , id.vars = c("cell" , "celltype") ,
                            measure.vars = gene ,
                            variable.name = "gene" ,
                            value.name = "Expr") %>%
  group_by(celltype , gene) %>%
  mutate(fillcolor = mean(Expr))

#install.packages("cowplot")
library(cowplot)
ggplot(data.melt , aes(factor(celltype) , Expr , fill = gene)) +
  geom_violin(scale = "width" , adjust = 1 , trim = T) +
  facet_grid(rows = vars(gene) , scales = "free" , switch = "y") +
  scale_y_continuous(expand = c(0 , 0) , position = "right" , labels = function(x)
    c(rep(x = "" , times = length(x) - 2) , x[length(x) - 1] , "")) +
  scale_fill_manual(values = c("#a2d2e7","#a2d2e7","#ffc17f","#ffc17f","#cf9f88","#cf9f88","#6fb3a3","#6fb3a3","#b3e19b","#b3e19b","#ff9d9f","#cdb6da","#cdb6da",
                               "#ff8831","#ff8831","#dba9a8")) +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none" , panel.spacing = unit(0 , "lines") ,
        plot.title = element_text(hjust = 0.5) ,
        panel.background = element_rect(fill = NA , colour = "black") ,
        plot.margin = margin(7 , 7 , 0 , 7 , "pt") ,
        strip.background = element_blank() ,
        strip.text = element_text(face = "bold") ,
        strip.text.y.left = element_text(angle = 0) ,
        axis.title.x = element_blank() ,
        axis.ticks.x = element_blank() ,
        axis.text.x = element_text(angle = 90 , hjust = 1 , vjust = 0.5 , color = "black")) +
  ggtitle("Feature on x-axis with annotation") +
  ylab("Expression Level")
ggsave("vln_stim.jpg",width = 18,height = 18,units = "cm",dpi = 1000)
dev.off()
par(mfrow = c(3 , 3) , xpd = T)
FeaturePlot(EYE , features = "VIM" ,reduction = "umap", cols = c("lightgrey" , "#A2D2E7")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("VIM_stim.jpg",width = 12,height = 12,units = "cm",dpi = 1000)

FeaturePlot(EYE , features = "PDE6A" ,reduction = "umap", cols = c("lightgrey" , "#FFC17F")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("PDE6A_stim.jpg",width = 12,height = 12,units = "cm",dpi = 1000)

FeaturePlot(EYE , features = "PDE6C" ,reduction = "umap", cols = c("lightgrey" , "#CF9F88")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("PDE6C_stim.jpg",width = 12,height = 12,units = "cm",dpi = 1000)


FeaturePlot(EYE , features = "TRPM1" ,reduction = "umap", cols = c("lightgrey" , "#6FB3A3")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("TRPM1_stim.jpg",width = 12,height = 12,units = "cm",dpi = 1000)

FeaturePlot(EYE , features = "GAD1" ,reduction = "umap", cols = c("lightgrey" , "#B3E19B")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("GAD1_stim.jpg",width = 12,height = 12,units = "cm",dpi = 1000)

FeaturePlot(EYE , features = "CNTN2" ,reduction = "umap", cols = c("lightgrey" , "#FF9D9F")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("CNTN2_stim.jpg",width = 12,height = 12,units = "cm",dpi = 1000)

FeaturePlot(EYE , features = "CALB2" ,reduction = "umap", cols = c("lightgrey" , "#CDB6DA")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("CALB2_stim.jpg",width = 12,height = 12,units = "cm",dpi = 1000)

FeaturePlot(EYE , features = "AIF1" ,reduction = "umap", cols = c("lightgrey" , "#FF8831")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("AIF1_stim.jpg",width = 12,height = 12,units = "cm",dpi = 1000)


FeaturePlot(EYE , features = "RPE65" ,reduction = "umap", cols = c("lightgrey" , "#DBA9A8")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("RPE65_stim.jpg",width = 12,height = 12,units = "cm",dpi = 1000)


FeaturePlot(EYE , features = "CNTN2" ,reduction = "umap", cols = c("lightgrey" , "red")) & NoLegend() & NoAxes() & theme(
  panel.border = element_rect(color = "black" , size = 1)
)
ggsave("CNTN2.jpg",width = 12,height = 12,units = "cm",dpi = 1000)












#分组细胞种类占比图
group=str_split(colnames(yizhi_seu@assays$RNA),"_",simplify=T)[,1]
data.frame(group)->group
yizhi_seu@meta.data$group.xi<-group
tongji <- table(yizhi_seu$celltype,yizhi_seu$group.xi) %>% reshape2::melt()
colnames(tongji) <- c("Cluster","Sample","Number")
tongji$Cluster <- factor(tongji$Cluster)                    #ת??Ϊ????
tongji$Percentage <- ave(tongji$Number,tongji$Sample,
                         FUN = function(x) x/sum(x))


prop_df <- table(xuanye$celltype , xuanye$stim) %>% reshape2::melt()
colnames(prop_df) <- c("Cluster" , "Sample" , "Number")
prop_df$Cluster <- factor(prop_df$Cluster)
prop_df$Proportion <- ave(prop_df$Number , prop_df$Sample , FUN = function(x) x/sum(x))


ggplot(prop_df,aes(Sample,Number,fill = Cluster))+
  geom_bar(stat = "identity",width = 0.8, position = "fill")+
  scale_fill_manual(values =c("#c1402b","#2a3663","#d9c5a8","#026868","#935741","#b9c4c7","#68b5be",
                              "#c37c8a","#454d58"))+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(size = 24,color = "black"))+
  theme(axis.text.y = element_text(size = 24,color = "black"))+
  theme(legend.text = element_text(size = 24,color = "black"),
        #legend.key.height = unit(80,"pt"),
        legend.position = "top",legend.box = "horizontal")+
  guides(fill=guide_legend(nrow = 1, byrow=T))
ggsave("percent.jpg",width = 48,height = 48,units = "cm",dpi = 1000)


####每组细胞的高表达基因####
Idents(xuanye) <- "celltype"
Marker <- FindAllMarkers(xuanye , logfc.threshold = 0.25 , only.pos = F)
colors <- c("#c1402b","#2a3663","#d9c5a8","#026868","#935741","#b9c4c7","#68b5be","#c37c8a","#454d58")
jjVolcano(diffData = marker2 , tile.col = colors ,
          aesCol = c("#c1402b" , "#2a3663"))

index <- grep("^ENSCPOG" , Marker$gene)
marker1 <- Marker
marker2 <- marker1[-index ,]
save(marker2 , file = "marker.Rdata")

microglia <- subset(marker2 , subset = marker2$cluster == "Microglia")


xuanye@meta.data$fenzu <- recode(xuanye$celltype ,
                                 "Microglia" = "Micro" ,
                                 "Astrocyte" = "other" ,
                                 "Oligodendrocyte" = "other" ,
                                 "Endothelial" = "other" ,
                                 "Red Blood" = "other" ,
                                 "Monocyte" = "other" ,
                                 "Neuron" = "other" ,
                                 "Ependymal" = "other" ,
                                 "T Cell" = "other")
my_comparisons <- list(c("other" , "Micro"))
VlnPlot(xuanye , features = "IL18" , group.by = "fenzu" , pt.size = 0) +
  theme_bw() +
  theme(axis.title.x = element_blank() ,
        axis.text.x = element_text(colour = "black" , face = "bold" , size = 12) ,
        axis.text.y = element_text(colour = "black" , face = "bold") ,
        axis.title.y = element_text(colour = "black" , face = "bold" , size = 15) ,
        panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank() ,
        panel.border = element_rect(color = "black" , size = 1.2 , linetype = "solid") ,
        panel.spacing = unit(0.12 , "cm") ,
        plot.title = element_text(hjust = 0.5 , face = "bold.italic") ,
        legend.position = "none") +
  stat_compare_means(method = "t.test" , hide.ns = F ,
                     comparisons = my_comparisons ,
                     label = "p.signif" ,
                     bracket.size = 0.8 ,
                     tip.length = 0 ,
                     size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.05 , 0.1))) +
  scale_fill_manual(values = c("#c1402b" , "#2a3663" ))

####亚群绘制####
load("Micro.Rdata")
Micro.tuxing <- Micro@reductions$umap@cell.embeddings
Micro.tuxing <- as.data.frame(Micro.tuxing)
Micro.tuxing %>%
  cbind(celltype=Micro@meta.data$seurat_clusters) -> Micro.tuxing
Micro.tuxing %>%
  cbind(stim=Micro@meta.data$stim) -> Micro.tuxing
Micro.tuxing$fenqun = recode(Micro.tuxing$celltype ,
                             "0" = "MG0" ,
                             "1" = "MG1" ,
                             "2" = "MG2" ,
                             "3" = "MG3" ,
                             "4" = "MG4" ,
                             "5" = "MG5" ,
                             "6" = "MG6" ,
                             "7" = "MG7")


####micro差异基因####
load("mi_LIM_NC.Rdata")
load("mi_EA_LIM.Rdata")

index <- grep("^ENSCPOG" , Mi_LIM_NC$GeneName)
Mi_LIM_NC <- Mi_LIM_NC[-index ,]
index <- grep("^ENSCPOG" , Mi_EA_LIM$GeneName)
Mi_EA_LIM <- Mi_EA_LIM[-index ,]

Mi_LIM_NC $ logp <- -log10(Mi_LIM_NC$p_val_adj)
Mi_EA_LIM $ logp <- -log10(Mi_EA_LIM$p_val_adj)

Mi_LIM_NC$group <- as.factor(ifelse(
  Mi_LIM_NC$p_val_adj <0.01 & abs(Mi_LIM_NC$avg_log2FC)>=0.4,
  ifelse(Mi_LIM_NC$avg_log2FC>=0.4,"up","down"),"NS"
))

Mi_EA_LIM$group <- as.factor(ifelse(
  Mi_EA_LIM$p_val_adj <0.01 & abs(Mi_EA_LIM$avg_log2FC)>=0.4,
  ifelse(Mi_EA_LIM$avg_log2FC>=0.4,"up","down"),"NS"
))

up_data <- filter(Mi_LIM_NC, group == 'up') %>%
  distinct(GeneName , .keep_all = TRUE) %>%
  top_n(10, logp)
down_data <- filter(Mi_LIM_NC, group == 'down') %>%
  distinct(GeneName , .keep_all = TRUE) %>%
  top_n(10, logp)

ggplot(Mi_LIM_NC,aes(avg_log2FC,logp),)+
  geom_point(aes(color=group),size=5)+
  scale_color_manual(values = c("#2a3663","#d9c5a8","#c1402b"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 24,color = "black"))+
  theme(axis.text.y = element_text(size = 24,color = "black"))+
  theme(legend.text = element_text(size = 24,color = "black"),
        #legend.key.height = unit(80,"pt"),
        legend.position = "top",legend.box = "horizontal")+
  guides(fill=guide_legend(nrow = 1, byrow=T))+
  theme(legend.title = element_blank())+

  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=3,col="black",lwd=0.8) +
  geom_point(data = up_data,
             aes(x = avg_log2FC, y = logp),
             color = 'red3', size = 4.5, alpha = 0.2) +
  geom_label_repel(data = up_data,
                   aes(x = avg_log2FC, y = logp, label = GeneName),
                   seed = 233,
                   size = 3.5,
                   color = 'black',
                   min.segment.length = 0,
                   force = 2,
                   force_pull = 2,
                   box.padding = 0.4,
                   max.overlaps = Inf) +
  geom_point(data = down_data,
             aes(x = avg_log2FC, y = logp),
             color = 'blue4', size = 4.5, alpha = 0.2) +
  geom_label_repel(data = down_data,
                   aes(x = avg_log2FC, y = logp, label = GeneName),
                   seed = 233,
                   size = 3.5,
                   color = 'black',
                   min.segment.length = 0,
                   force = 2,
                   force_pull = 2,
                   box.padding = 0.4,
                   max.overlaps = Inf)
ggsave("micro_huoshan1.jpg",width = 48,height = 48,units = "cm",dpi = 1000)

up_data2 <- filter(Mi_EA_LIM, group == 'up') %>%
  distinct(GeneName , .keep_all = TRUE) %>%
  top_n(10, logp)
down_data2 <- filter(Mi_EA_LIM, group == 'down') %>%
  distinct(GeneName , .keep_all = TRUE) %>%
  top_n(10, logp)

ggplot(Mi_EA_LIM,aes(avg_log2FC,logp),)+
  geom_point(aes(color=group),size=5)+
  scale_color_manual(values = c("#2a3663","#d9c5a8","#c1402b"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 24,color = "black"))+
  theme(axis.text.y = element_text(size = 24,color = "black"))+
  theme(legend.text = element_text(size = 24,color = "black"),
        #legend.key.height = unit(80,"pt"),
        legend.position = "top",legend.box = "horizontal")+
  guides(fill=guide_legend(nrow = 1, byrow=T))+
  theme(legend.title = element_blank())+

  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=3,col="black",lwd=0.8) +
  geom_point(data = up_data2,
             aes(x = avg_log2FC, y = logp),
             color = 'red3', size = 4.5, alpha = 0.2) +
  geom_label_repel(data = up_data2,
                   aes(x = avg_log2FC, y = logp, label = GeneName),
                   seed = 233,
                   size = 3.5,
                   color = 'black',
                   min.segment.length = 0,
                   force = 2,
                   force_pull = 2,
                   box.padding = 0.4,
                   max.overlaps = Inf) +
  geom_point(data = down_data2,
             aes(x = avg_log2FC, y = logp),
             color = 'blue4', size = 4.5, alpha = 0.2) +
  geom_label_repel(data = down_data2,
                   aes(x = avg_log2FC, y = logp, label = GeneName),
                   seed = 233,
                   size = 3.5,
                   color = 'black',
                   min.segment.length = 0,
                   force = 2,
                   force_pull = 2,
                   box.padding = 0.4,
                   max.overlaps = Inf)
ggsave("micro_huoshan2.jpg",width = 48,height = 48,units = "cm",dpi = 1000)


####muller亚群####
tsnedata <- Muller@reductions$umap@cell.embeddings
tsnedata <- as.data.frame(tsnedata)
tsnedata %>%
  cbind(celltype=Muller@meta.data$zhushi) ->tsnedata
tsnedata %>%
  cbind(stim=Muller@meta.data$stim) ->tsnedata
?geom_point
#umap图
celltype_med <- tsnedata %>% group_by(celltype) %>% summarise(umap_1 = median(umap_1) , umap_2 = median(umap_2))
ggplot(tsnedata,aes(umap_1,umap_2,colour=celltype))+
  geom_point(size=0.5)+
  #labs(title = "umap")+
  theme_bw()+
  xlab("")+
  ylab("")+
  #stat_ellipse(type = "norm",linetype=2)+
  scale_color_manual(values =c("#a2d2e7","#ffc17f","#cf9f88","#6fb3a3","#b3e19b","#ff9d9f","#cdb6da","#ff8831","#dba9a8","#50aa4b","#704ba3"
  ))+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 32,
                                  face = "bold",
                                  hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 20,
                                   face = "bold"),
        legend.background = element_blank(),
        legend.key.height = unit(40,"pt"),
        legend.position = "none") +
  geom_label_repel(aes(label = celltype) , fontface = "bold" , data = celltype_med ,
                   point.padding = unit(0.5 , "lines")) +
  geom_segment(aes(x = min(tsnedata$umap_1) , y = min(tsnedata$umap_2) ,
                   xend = min(tsnedata$umap_1) + 3 , yend = min(tsnedata$umap_2) ) ,
               colour = "black" , size = 1 , arrow = arrow(length = unit(0.3 , "cm"))) +
  geom_segment(aes(x = min(tsnedata$umap_1) , y = min(tsnedata$umap_2) ,
                   xend = min(tsnedata$umap_1) , yend = min(tsnedata$umap_2) + 3) ,
               colour = "black" , size = 1 , arrow = arrow(length = unit(0.3 , "cm"))) +
  annotate("text" , x = min(tsnedata$umap_1) + 1.5 , y = min(tsnedata$umap_2) - 1 , label = "UMAP-1" ,
           color = "black" , size = 3 , fontface = "bold") +
  annotate("text" , x = min(tsnedata$umap_1) - 1 , y = min(tsnedata$umap_2) + 1.5 , label = "UMAP-2" ,
           color = "black" , size = 3 , fontface = "bold" , angle = 90)
ggsave("umap_stim.jpg",width = 18,height = 18,units = "cm",dpi = 1000)

