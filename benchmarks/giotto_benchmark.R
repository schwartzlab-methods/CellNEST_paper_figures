# https://giottosuite.readthedocs.io/en/master/giottoworkflowanalyses.html?highlight=spatCellCellcom 
library(Giotto)
library(yaml)
library(reticulate)
library(Matrix)
library(dplyr)
library(tibble)
library(rlang)
library(jsonlite)

load_obj <- function(data_path) {
    counts <- t(as.matrix(read.csv(file.path(data_path, "counts.csv"), row.names = 1, check.names = F)))
    spatial <- read.csv(file.path(data_path, "coords.csv"), row.names = 1, check.names = F)
    giotto_obj <- createGiottoObject(
        # matrix with raw expression counts
        raw_exprs = counts,
        # data.table or data.frame: rownames: spot names, colnames: ("x", "y") 
        spatial_locs = spatial
    )
    giotto_obj@gene_ID <- as.character(giotto_obj@gene_ID)
    # add clustering results 
    meta <- read.csv(file.path(data_path, "leiden.csv"))
    meta$spot <- as.character(meta$spot)
    meta$cluster <- as.character(meta$cluster)
    giotto_obj <- addCellMetadata(gobject = giotto_obj, new_metadata = meta, by_column = TRUE, column_cell_ID = "spot")

    return(giotto_obj)
}

computeMaxDist <- function(
    dataPath,
    # range for communication in uM 
    range = 250
) {
    spot_size = 55
    scale_factors <- fromJSON(file.path(dataPath, "spatial", "scalefactors_json.json"))
    conversion_factor = spot_size / scale_factors$spot_diameter_fullres # 0.6111 um/pixel
    pix_per_um = (1 / conversion_factor)
    dist = pix_per_um * range
    # return distance in pixels
    return (dist)
}


run_giotto <- function(
    data_path,
    database_path,
    synthetic,
    output_path,
    result_dir
) {
    giotto_obj <- load_obj(data_path)
    print(giotto_obj@gene_ID)

    giotto_obj <- normalizeGiotto(gobject = giotto_obj) # default scale factor = 6000
    giotto_obj <- calculateHVG(gobject = giotto_obj)
    giotto_obj <- runPCA(gobject = giotto_obj)
    giotto_obj <- runUMAP(giotto_obj, dimensions_to_use = 1:5)
    # read in database 
    LR_data <- read.csv(database_path, check.names = F)
    cols_to_remove <- c("Annotation", "Reference")
    LR_data <- LR_data[, !(colnames(LR_data) %in% cols_to_remove)]
    colnames(LR_data) <- tolower(colnames(LR_data))
    LR_data$ligand <- as.character(LR_data$ligand)
    LR_data$receptor <- as.character(LR_data$receptor)
    LR_data$ligand_det <- LR_data$ligand %in% giotto_obj@gene_ID
    LR_data$receptor_det <- LR_data$receptor %in% giotto_obj@gene_ID
    print(LR_data)
    LR_data_det <- LR_data[LR_data$ligand_det & LR_data$receptor_det, ]
    select_ligands = LR_data_det$ligand
    select_receptors = LR_data_det$receptor
    if (synthetic) {
        dataset_name <- basename(data_path)
        if (grepl("mixed", dataset_name)) {
            maximum_distance = 30
        } else {
            maximum_distance = 4
        }
    } else {
        maximum_distance = computeMaxDist(data_path)
    }
    giotto_obj = createSpatialNetwork(gobject = giotto_obj, name = 'Delaunay_network', maximum_distance_delaunay = maximum_distance, minimum_k = 2)
    giotto_obj = createSpatialNetwork(gobject = giotto_obj, name = "knn_network", method = "kNN", maximum_distance_knn = maximum_distance)
    # spatial CCC
    spatial_all_scores = spatCellCellcom(
        giotto_obj,
        spatial_network_name = 'Delaunay_network',
        cluster_column = 'cluster',
        random_iter = 500,
        gene_set_1 = select_ligands,
        gene_set_2 = select_receptors,
        adjust_method = 'fdr',
        do_parallel = T,
        cores = 5,
        verbose = 'a lot'
    )
    print(length(spatial_all_scores))
    # select top lr 
    # lig_nr: total number of cells from lig_cell_type that spatially interact with cells from rec_cell_type
    # rec_nr: total number of cells from rec_cell_type that spatially interact with cells from lig_cell_type
    selected_spat <- spatial_all_scores[p.adj <= 0.5 & abs(log2fc) > 0.1 & lig_nr > 2 & rec_nr > 2]
    write.csv(selected_spat, output_path)
}

main <- function() {
    config <- yaml::read_yaml("config.yml")
    data_dir <- file.path(config$directories$data)
    result_dir <- file.path(config$directories$raw_result_user, "giotto")
    dir.create(result_dir, recursive = TRUE)

    datasets = config$params$datasets
    datasets <- datasets[grepl("mechanistic_random_ccc_wo_relay", datasets)]
    # datasets <- c("lymph_node")
    # datasets <- c("mixed_distribution_mechanistic_low_noise")

    for (dataset in datasets) {
        data_path = file.path(data_dir, dataset)
        output_path = file.path(result_dir, paste0(dataset, ".csv"))
        if (dataset == "lymph_node") {
            synthetic = FALSE
            database_path = config$files$nest_database
        } else {
            synthetic = TRUE
            database_path = list.files(path = data_path, pattern = ".*lrdb\\.csv$", full.names = TRUE)[1]
        }
        run_giotto(data_path, database_path, synthetic, output_path, result_dir)
    }
}

main()
