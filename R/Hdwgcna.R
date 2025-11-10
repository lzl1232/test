#### 安装包 ####

# options(timeout=10000)
# BiocManager::install(c(
#   'JASPAR2020',
#   'EnsDb.Hsapiens.v86',
#   'BSgenome.Hsapiens.UCSC.hg38'
# ))
# BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
# BiocManager::install(c(
#   'motifmatchr',
#   'TFBSTools',
#   'GenomicRanges'
# ))
# install.packages('xgboost')
#
# devtools::install_github('smorabit/hdWGCNA', ref='dev')
# BiocManager::install("GeneOverlap")
# remotes::install_github('satijalab/seurat-wrappers')
library(Seurat)
library(CytoTRACE2)
library(gridExtra)
library(patchwork)
library(ggplot2)
library(ggpubr)
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
library(motifmatchr)
library(TFBSTools)
library(GenomicRanges)
library(GeneOverlap)
library(JASPAR2020)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(xgboost)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(monocle3)
library(Seurat)
library(SeuratWrappers)




Seurat_obj <- SetupForWGCNA(
  Muller,
  gene_select = 'fraction',
  fraction = 0.05,
  wgcna_name = 'tutorial'
)

Seurat_obj <- MetacellsByGroups(
  seurat_obj = Seurat_obj,
  group.by = c('zhushi', 'stim'),
  k = 25,
  # target_metacells = 250,
  ident.group = 'zhushi',
  # min_cells = 70,
  max_shared = 10
)

Seurat_obj <- NormalizeMetacells(Seurat_obj)
metacell_obj <- GetMetacellObject(Seurat_obj)
Seurat_obj <- ScaleMetacells(Seurat_obj, features = VariableFeatures(Seurat_obj))
Seurat_obj <- RunPCAMetacells(Seurat_obj, features = VariableFeatures(Seurat_obj))
Seurat_obj <- RunHarmonyMetacells(Seurat_obj, group.by.vars = 'stim')
Seurat_obj <- RunUMAPMetacells(Seurat_obj, reduction = 'harmony', dims = 1:30)
DimPlotMetacells(Seurat_obj, group.by = 'zhushi')
DimPlotMetacells(Seurat_obj, group.by = 'stim')



Seurat_obj <- SetDatExpr(
  Seurat_obj,
  group.by = 'zhushi',
  group_name = '1',
  assay = 'RNA',
  slot = 'data'
)

Seurat_obj <- TestSoftPowers(
  Seurat_obj,
  networkType = 'signed')
plot_list <- PlotSoftPowers(Seurat_obj)
wrap_plots(plot_list, ncol = 2)
power_table <- GetPowerTable(Seurat_obj)
head(power_table)


Seurat_obj <- ConstructNetwork(
  Seurat_obj,
  soft_power = 10,
  tom_name = '1',
  overwrite_tom = T,
  SetDatExpr = F
)
## 图一
PlotDendrogram(Seurat_obj, main = 'hdWGCNA Dendrogram')


# 基因模块组成
Seurat_obj@misc$tutorial$wgcna_modules %>% head
table(Seurat_obj@misc$tutorial$wgcna_modules$module)

TOM <- GetTOM(Seurat_obj)
Seurat_obj <- ScaleData(Seurat_obj, features = VariableFeatures(Seurat_obj))
Seurat_obj <- ModuleEigengenes(
  Seurat_obj,
  group.by.vars = 'stim'
)

hMES <- GetMEs(Seurat_obj)
head(hMES)
hMES <- GetMEs(Seurat_obj, harmonized = F)
Seurat_obj <- ModuleConnectivity(
  Seurat_obj,
  group.by = 'zhushi',
  group_name = '1',
)


## 图2
PlotKMEs(Seurat_obj, ncol = 3) + theme_bw()

##
modules <- GetModules(Seurat_obj)
save(modules, file = 'modules.Rdata')
head(modules[, 1:6])
##hub
hub_df <- GetHubGenes(Seurat_obj, n_hubs = 40)  ## kme > 6.0
head(hub_df)

Seurat_obj <- ModuleExprScore(
  Seurat_obj,
  n_genes = 25,
  method = 'Seurat')


## 图3
plot_list2 <- ModuleFeaturePlot(
  Seurat_obj,
  features = 'hMEs',
  order = T
)
wrap_plots(plot_list2)
## 图4
plot_list3 <- ModuleFeaturePlot(
  Seurat_obj,
  features = 'scores',
  order = 'shuffle',
  ucell = T
)
wrap_plots(plot_list3)


hMES <- GetMEs(Seurat_obj, harmonized = T)
mods <- colnames(hMES)
Seurat_obj@meta.data <- cbind(Seurat_obj@meta.data, hMES)

Seurat_obj$zhushi <- recode(Seurat_obj@meta.data$zhushi ,
                        "0" = "MC_01" ,
                        "1" = "MC_02" ,
                        "2" = "MC_03" ,
                        "3" = "MC_04" ,
                        "4" = "MC_05" )
## 图5
DotPlot(Seurat_obj, features = mods, group.by = 'zhushi')+
  theme_bw()+
  theme( axis.text.x=element_text(hjust = 1,vjust=0.5, angle = 90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
## 图6
DotPlot(Seurat_obj, features = mods, group.by = 'stim') +
  theme_bw()+
  theme( axis.text.x=element_text(hjust = 1,vjust=0.5, angle = 90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))



