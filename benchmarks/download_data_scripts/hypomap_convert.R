library(Seurat)
library(SeuratDisk)
library(tools)

source("/cluster/projects/schwartzgroup/deisha/projects/rnaSeqHelpers/helpers.R")
rds_path <- "data/merfish/hypomap.R"

if (file.exists(rds_path)) {
    rds_to_adata(rds_path)
} else {
    print("Conversion failed")
}

