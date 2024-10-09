library(TWCOM)
library(Seurat)
library(yaml)
library(jsonlite)
library(parallel)

prepData <- function(
    dataPath
) {
    # read genes x spots/cells 
    STdata <- t(as.matrix(read.csv(file.path(dataPath, "counts.csv"), row.names = 1, check.names = F)))
    colnames(STdata) <- as.character(colnames(STdata))
    dat <- CreateSeuratObject(
        counts = STdata 
        # meta.data = STmeta
    )
    dat <- NormalizeData(
        dat,
        assay = "RNA", 
        scale.factor = 10000
    )
    stdat <- dat@assays$RNA$data
}

# find max distance in pixels for juxtacrine communication in Visium data 
computeMaxDist <- function(
    dataPath, 
    # contact range for communication in uM 
    contact_range = 250
) {
    spot_size = 55
    scale_factors <- fromJSON(file.path(dataPath, "spatial", "scalefactors_json.json"))
    conversion_factor = spot_size / scale_factors$spot_diameter_fullres # 0.6111 um/pixel 
    pix_per_um = (1 / conversion_factor)
    dist = pix_per_um * contact_range
    # return distance in pixels
    return (dist)
}

runTWCOM <- function(
    dataPath,
    dbPath,
    annotPath,
    outputPath,
    synthetic
) {
    # normalized matrix 
    expMat <- prepData(dataPath)
    # coordinates 
    coords <- read.csv(file.path(dataPath, "coords.csv"), row.names = 1, check.names = F)
    # order rows of dataframe to match expression matrix column order
    print(head(coords))
    STmeta <- coords[match(colnames(expMat), rownames(coords)), ]
    print(head(STmeta))
    annotDf <- read.csv(annotPath, row.names = 1, check.names = F)
    print(annotDf)
    # rownames(annotDf) <- as.character(rownames(annotDf))
    if (synthetic) {
        annotDf$cluster <- as.character(annotDf$cluster)
    }
    # load ligand-receptor db 
    dbDf <- read.csv(dbPath)
    cols_to_remove <- c("Annotation", "Reference")
    dbDf <- dbDf[, !(colnames(dbDf) %in% cols_to_remove)]
    colnames(dbDf) <- tolower(colnames(dbDf))
    dbDf$ligand <- as.character(dbDf$ligand)
    dbDf$receptor <- as.character(dbDf$receptor)
    ligands <- dbDf$ligand
    receptors <- dbDf$receptor

    if (synthetic) {
        # single-cell annotations for single-cell data
        # creates matrix, preserving row ordering from stmeta
        M <- CellTypeMatrix(
            stmeta = annotDf, 
            celltype = "cluster"
        )
        # neighbour distance
        maxDist = 4
    } else {
        print(head(annotDf))
        # cell type proportions for spot-based data
        M <- as.matrix(annotDf)
        print(M[1:10, 1:10])
        # neighbour distance in pixels (250 uM)
        maxDist = computeMaxDist(dataPath, contact_range = 250)
    }

    print("M created")
    rownames(M) <- rownames(annotDf)
    print(M)
    print("cell type matrix")
    print("expression matrix")
    print(expMat[1:10, 1:10])
    # order rows of matrix to match expression matrix column order
    M <- M[match(colnames(expMat), rownames(M)), ]

    print(dim(expMat)) 
    print(dim(M))
    print("cell type matrix after re-ordering")
    # print(M[1:10, 1:10])

    restab <- FRETCOM(
        stdat = expMat,
        M = M, 
        ligands = ligands, 
        receptors = receptors,
        coordx = STmeta$x, 
        coordy = STmeta$y, 
        maxdist = maxDist
    )
    # filter by default q-value = 0.2 
    filtered <- subset(restab, Qvalue < 0.2)
    filtered$lri <- paste(filtered$Ligand, filtered$Receptor, sep = "-") 
    write.csv(filtered, outputPath)
}

process_dataset <- function(dataset, config) {
    data_dir <- file.path(config$directories$data)
    result_dir <- file.path(config$directories$raw_result, "TWCOM")
    dataPath <- file.path(data_dir, dataset)
    outputPath <- file.path(result_dir, paste0(dataset, ".csv"))
    
    if (dataset == "lymph_node") {
        synthetic <- FALSE
        dbPath <- config$files$nest_database
        annotPath <- file.path(data_dir, "cytospace", dataset, "fractional_abundances_by_spot.csv")
    } else {
        synthetic <- TRUE
        dbPath <- list.files(path = dataPath, pattern = ".*lrdb\\.csv$", full.names = TRUE)[1]
        annotPath <- file.path(dataPath, "leiden.csv")
    }
    
    runTWCOM(dataPath, dbPath, annotPath, outputPath, synthetic)
}

main <- function() {
    config <- yaml::read_yaml("config.yml")
    result_dir <- file.path(config$directories$raw_result, "TWCOM")
    dir.create(result_dir, recursive = TRUE)

    datasets <- config$params$datasets
    
    # Use mclapply for parallel processing
    mclapply(datasets, process_dataset, config, mc.cores = detectCores() - 1)  # use all but one core
}

main()









