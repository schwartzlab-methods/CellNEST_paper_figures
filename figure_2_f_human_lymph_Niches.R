library(Seurat)             
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(NICHES)
library(viridis)
############## Niches on Lymph Node #######################
data_dir <- 'data/V1_Human_Lymph_Node_spatial/'
list.files(data_dir)
seurat_object <- Load10X_Spatial(data.dir = data_dir)
lymph <- SCTransform(seurat_object, assay = "Spatial", verbose = FALSE)
lymph <- RunPCA(lymph, assay = "SCT", verbose = FALSE)
lymph <- FindNeighbors(lymph, reduction = "pca", dims = 1:30)
lymph <- FindClusters(lymph, verbose = FALSE)
lymph <- RunUMAP(lymph, reduction = "pca", dims = 1:30)
p1 <- DimPlot(lymph, reduction = "umap",group.by = 'seurat_clusters', label = TRUE)
p2 <- SpatialDimPlot(lymph, label = TRUE,group.by = 'seurat_clusters', label.size = 3)
ggsave("Niches_V1_Human_Lymph_Node_spatial_myplot.png", plot = (p1+p2))

lymph@meta.data$x <- lymph@images$slice1@coordinates$row
lymph@meta.data$y <- lymph@images$slice1@coordinates$col

DefaultAssay(lymph) <- "Spatial"
lymph <- NormalizeData(lymph)

lymph <- SeuratWrappers::RunALRA(lymph)
lr_db <- read.csv("NEST_figures_input_PDAC/lr_cellchat_nichenet.csv")
NICHES_output <- RunNICHES(object = lymph,
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
ec.network <- subset(cc.object,idents ='3')
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
ggsave("Niches_V1_Human_Lymph_Node_spatial_LR_pairs.png", plot = p)

# https://msraredon.github.io/NICHES/articles/03%20Rat%20Alveolus.html
# https://msraredon.github.io/NICHES/articles/01%20NICHES%20Spatial.html

