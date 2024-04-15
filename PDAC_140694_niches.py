library(Seurat)             
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(NICHES)
library(viridis)

############## Niches on PDAC 140694 #######################
data_dir <- '/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/exp1_C1/outs/'
list.files(data_dir)
seurat_object <- Load10X_Spatial(data.dir = data_dir)
pdac <- SCTransform(seurat_object, assay = "Spatial", verbose = FALSE)
pdac <- RunPCA(pdac, assay = "SCT", verbose = FALSE)
pdac <- FindNeighbors(pdac, reduction = "pca", dims = 1:30)
pdac <- FindClusters(pdac, verbose = FALSE)
pdac <- RunUMAP(pdac, reduction = "pca", dims = 1:30)

p1 <- DimPlot(pdac, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(pdac, label = TRUE,  label.size = 3)
# Error in FUN(left, right) : non-numeric argument to binary operator
ggsave("/cluster/home/t116508uhn/myplot.png", plot = p2)

pdac@meta.data$x <- pdac@images$slice1@coordinates$row
pdac@meta.data$y <- pdac@images$slice1@coordinates$col

DefaultAssay(pdac) <- "Spatial"
pdac <- NormalizeData(pdac)

pdac <- SeuratWrappers::RunALRA(pdac)
lr_db <- read.csv("/cluster/home/t116508uhn/64630/lr_cellchat_nichenet.csv")
NICHES_output <- RunNICHES(object = pdac,
                           LR.database = "custom",
                           custom_LR_database = lr_db,
                           species = "human",
                           assay = "alra",
                           position.x = 'x',
                           position.y = 'y',
                           k = 12, 
                           cell_types = "seurat_clusters",
                           min.cells.per.ident = 0,
                           min.cells.per.gene = NULL,
                           meta.data.to.map = c('orig.ident','seurat_clusters'),
                           CellToCell = F,CellToSystem = F,SystemToCell = F,
                           CellToCellSpatial = T,CellToNeighborhood = F,NeighborhoodToCell = F)
                           
                           
niche <- NICHES_output[['CellToCellSpatial']]
Idents(niche) <- niche[['ReceivingType']]

cc.object <- NICHES_output$CellToCellSpatial #Extract the output of interest
cc.object <- ScaleData(cc.object) #Scale
cc.object <- FindVariableFeatures(cc.object,selection.method="disp") #Identify variable features
cc.object <- RunPCA(cc.object,npcs = 100) #RunPCA
cc.object <- RunUMAP(cc.object,dims = 1:100)

Idents(cc.object) <- cc.object[['ReceivingType']]
ec.network <- subset(cc.object, idents ='7')
Idents(ec.network) <- ec.network[['VectorType']]
mark.ec <- FindAllMarkers(ec.network,
                          logfc.threshold = 1,
                          min.pct = 0.5,
                          only.pos = T,
                          test.use = 'roc')
# Pull markers of interest to plot
mark.ec$ratio <- mark.ec$pct.1/mark.ec$pct.2
marker.list.ec <- mark.ec %>% group_by(cluster) %>% top_n(5,avg_log2FC)
p <- DoHeatmap(ec.network,features = marker.list.ec$gene,cells = WhichCells(ec.network,downsample = 100))
ggsave("/cluster/home/t116508uhn/64630/myplot.png", plot = p)




