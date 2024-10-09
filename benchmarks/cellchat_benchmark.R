ptm = Sys.time()

library(CellChat)
library(patchwork)
library(Seurat)
library(hdf5r)
library(jsonlite)
library(dplyr)
library(tibble)

# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.html 

create_obj <- function(
    synthetic, 
    data_path
) {
    data.input <- t(as.matrix(read.csv(file.path(data_path, "counts.csv"), row.names = 1, check.names = F)))
    if (synthetic) {
        spatial.locs <- read.csv(file.path(data_path, "coords.csv"), row.names = 1)
        # ratio = 1 - do not scale pixel distance 
        # tol - tolerance factor to increase the robustness when comparing the center-to-center distance against the interaction.range
        # tol recommended to be half value of cell/spot size in the unit of um
        # tol = 1 approximates center-center distance in 150 x 150 pixel array 
        spatial.factors <- data.frame(ratio = 1, tol = 1)
    } else {
        # locations from full (NOT high/low) resolution images
        visium_obj <- Load10X_Spatial(data.dir = data_path)
        spatial.locs <- GetTissueCoordinates(visium_obj, scale = NULL, cols = c("imagerow", "imagecol"))
        spot.size = 55
        scalefactors = fromJSON(file.path(data_path, 'spatial', 'scalefactors_json.json'))
        conversion.factor = spot.size / scalefactors$spot_diameter_fullres
        spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
    }

    # leiden clusters 
    meta <- read.csv(file.path(data_path, "leiden.csv"), row.names = 1, check.names = F)
    meta$cluster <- as.character(meta$cluster)
    # replace "0" with "zero": error - labels cannot = 0
    meta$cluster <- ifelse(meta$cluster == "0", "zero", meta$cluster) 
    
    # normalize & log-transform
    # print(data.input) 
    data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)

    cellchat <- createCellChat(
        object = data.input, # genes x cells  
        meta = meta, 
        group.by = "cluster",
        datatype = "spatial", 
        coordinates = spatial.locs, 
        spatial.factors = spatial.factors
    )
    
    return (cellchat) 
}

run_cellchat <- function(
    data_path,
    database_path,
    synthetic, 
    output_path
) {
    cellchat <- create_obj(synthetic, data_path)

    if (synthetic) {
        # https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChatDB.html#step-2-formulate-the-input-files-to-be-compatible-with-cellchatdb 
        # https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.html 
        db.user <- read.csv(database_path)
        db.user$ligand <- as.character(db.user$ligand)
        db.user$receptor <- as.character(db.user$receptor)
        db.user$interaction_name <- paste0(db.user$ligand, "_", db.user$receptor)
        db.user$cofactor <- ""
        db.user$complex <- ""
        # assume all synthetic interactions are secreted signaling, not contact-based 
        db.user$annotation <- "Secreted Signaling"
        print(head(db.user))
        db_genes <- union(db.user$ligand, db.user$receptor)
        print(db_genes)
        gene_info <- data.frame(Symbol = db_genes)
        print("printing gene info")
        print(head(gene_info))
        # update database
        db.new <- updateCellChatDB(db = db.user, gene_info = gene_info)
        CellChatDB <- db.new
        cellchat@DB <- CellChatDB
        cellchat <- subsetData(cellchat, features = rownames(cellchat@data))        
    } else {
        CellChatDB <- CellChatDB.human
        # exclude non-protein CCC
        CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), key = "annotation")
        cellchat@DB <- CellChatDB.use
        cellchat <- subsetData(cellchat)
    }
    dplyr::glimpse(CellChatDB$interaction)
    future::plan("multisession", workers = 4) 
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
    if (synthetic) {
        interaction.range = 4
        # scale.distance: we choose this values such that the minimum value of the scaled distances is [1,2]
        # this value can be 1, 0.1, 0.01, 0.001
        scale.distance = 1
        contact.dependent = FALSE
        contact.range = NULL
    } else {
        interaction.range = 250
        scale.distance = 0.01
        contact.dependent = TRUE 
        contact.range = 100
    }
    # can only set number of k nn as an alternative to contact.range, not for interaction.range 
    # address for mixed distribution  
    cellchat <- computeCommunProb(
        cellchat, 
        type = "triMean", 
        raw.use = TRUE, 
        distance.use = TRUE, 
        interaction.range = interaction.range, 
        scale.distance = scale.distance,
        contact.dependent = contact.dependent, 
        contact.range = contact.range
    )
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    # subset inferred ccc of interest using default p-value threshold of 0.05
    df.net <- subsetCommunication(cellchat, thresh = 0.05)
    df.net$source <- as.character(df.net$source)
    df.net$target <- as.character(df.net$target)
    df.net$source[df.net$source == "zero"] <- 0
    df.net$target[df.net$target == "zero"] <- 0
    # filter top prob (strength) 20%
    # prob.thresh <- quantile(df.net$prob, 0.80)
    # df.net.top <- df.net %>% filter(prob >= prob.thresh)
    # write results
    write.csv(df.net, output_path)
}

main <- function() {
    config <- yaml::read_yaml("config.yml")
    data_dir <- file.path(config$directories$data)
    result_dir <- file.path(config$directories$raw_result_user, "cellchat")
    dir.create(result_dir, recursive = TRUE)

    datasets <- config$params$datasets

    for (dataset in datasets) {
        data_path = file.path(data_dir, dataset)
        output_path = file.path(result_dir, paste0(dataset, ".csv"))
        if (dataset == "lymph_node") {
            synthetic = FALSE
            database_path = NA # use default CellChat DB
        } else {
            synthetic = TRUE
            database_path = list.files(path = data_path, pattern = ".*lrdb\\.csv$", full.names = TRUE)[1]
        }
        run_cellchat(data_path, database_path, synthetic, output_path) 
    }
}

main()
