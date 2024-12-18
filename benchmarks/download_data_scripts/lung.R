library(SeuratDisk)

# write seurat obj of LUAD visium TD1 to anndata object 
dir <- "/mnt/data0/dpaliwal/software/NEST_paper_figures/benchmarks/data/lung"
setwd(dir)
basename = "seurat_ST1"
load(paste0(basename, ".rda"))
saveRDS(seurat, paste0(basename, ".rds"))
seurat_obj <- readRDS(basename, ".rds")
SaveH5Seurat(seurat_obj, filname = paste0(basename, ".h5Seurat"))
Convert(paste0(basename, ".h5Seurat"), dest = "h5ad")
