library(yaml)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(scuttle)

# atlas data acquired: https://github.com/rmoffitt/scOh 
# paper ref: https://www.nature.com/articles/s41467-023-40895-6

process_epi <- function(moffitt_data_dir, moffitt_repo, procd_dir) {
    epi_seur <- readRDS(file.path(moffitt_data_dir, "Epi_Subset2021.rds"))
    load(file.path(moffitt_data_dir, "gene_lists.RData"))
    print(gene_lists)
    source(file.path(moffitt_repo, "R", "score_sig.R"), echo = TRUE, local = TRUE)
    source(file.path(moffitt_repo, "R", "pseudobulk_score.R"), echo = TRUE, local = TRUE)
    source(file.path(moffitt_repo, "R", "pseudobulk.R"), echo = TRUE, local = TRUE)
    epi_scored <- score_sig(seurat = epi_seur, gene_lists = gene_lists)
    saveRDS(epi_scored, file.path(procd_dir, "Epi_Subset2021_Scored.rds"))
}

merge_atlas <- function(moffitt_data_dir, procd_dir, atlas_rds) {
    end <- readRDS(file.path(moffitt_data_dir, "Endo_Subset2021.rds"))
    lym <- readRDS(file.path(moffitt_data_dir, "Lymphocytes_Subset2021.rds"))
    myl <- readRDS(file.path(moffitt_data_dir, "Myeloid_Subset2021.rds"))
    str <- readRDS(file.path(moffitt_data_dir, "Stroma_Subset2021.rds"))
    epi <- readRDS(file.path(procd_dir, "Epi_Subset2021_Scored.rds"))
    combined <- merge(end, y = c(epi, lym, myl, str))
    saveRDS(combined, file = file.path(procd_dir, atlas_rds))
}

assign_neoplastic <- function(procd_dir, atlas_rds) {
    obj <- readRDS(file.path(procd_dir, atlas_rds))
    subtype_dict <- list()
    subtype_dict[["moffitt"]] <- c("Moffitt.Basal.25", "Moffitt.Classical.25")
    subtype_dict[["notta"]] <- c("Notta.ClassicalA", "Notta.ClassicalB", "Notta.BasalA", "Notta.BasalB")
    subtype_dict[["raghavan"]] <- c("RaghavanIntermediate", "RaghavanBasal", "RaghavanClassical")
    subtype_dict[["puleo"]] <- c("Puleo.PureClassical", "Puleo.ImmuneClassical", "Puleo.PureBasal")
    neoplastic_indices <- which(obj$CellType2 == "Neoplastic")
    for (scheme in names(subtype_dict)) {
        new_col <- paste0("max_", scheme)
        obj[[new_col]] <- obj$CellType2
        cols <- subtype_dict[[scheme]]
        values <- obj@meta.data[neoplastic_indices, cols]
        max_col_index <- apply(values, 1, which.max)
        max_col_names <- cols[max_col_index]
        obj@meta.data[[new_col]] <- as.character(obj@meta.data[[new_col]])
        obj@meta.data[neoplastic_indices, new_col] <- max_col_names
        print(table(obj[[new_col]]))
    }
    file_name_no_ext <- tools::file_path_sans_ext(atlas_rds)
    SaveH5Seurat(obj, filename = file.path(procd_dir, paste0(file_name_no_ext, ".h5Seurat")))

    # save to h5ad & rds
    Convert(file.path(procd_dir, paste0(file_name_no_ext, ".h5Seurat")), dest = "h5ad")
    saveRDS(obj, file.path(procd_dir, atlas_rds))
}


main <- function() {
    config <- yaml::read_yaml("../config.yml")
    raw_dir <- config$directories$raw
    procd_dir <- config$directories$processed
    moffitt_repo <- config$directories$moffitt
    moffitt_data_dir <- file.path(raw_dir, "moffitt")

    atlas_rds <- "Moffit_Combined.rds"

    process_epi(moffitt_data_dir, moffitt_repo, procd_dir)
    merge_atlas(moffitt_data_dir, procd_dir, atlas_rds)
    assign_neoplastic(procd_dir, atlas_rds)
}

main()
