library("DESeq2", quietly = TRUE)
library(ggplot2)
library(dplyr)
library(argparse)

run_deseq2 <- function(data.dir,
                       counts.file,
                       meta.file,
                       data.type,
                       genes) {
  # convert genes string to column vector
  genes_vec <- unlist(strsplit(genes, ","))
  print(genes_vec)

  mat.path <- file.path(data.dir, counts.file)
  col.path <- file.path(data.dir, meta.file)

  # load counts
  counts <- as.matrix(read.csv(mat.path,
    header = TRUE,
    sep = ",",
    row.names = "genes"
  ))
  # load metadata
  coldata <- read.csv(col.path,
    row.names = 1
  )

  # compare Classical vs. Basal-like 
  coldata[["subtype"]] <- factor(coldata[["subtype"]],
      levels = c("Basal-like", "Classical")
  )

  # check formatting
  print(all(rownames(coldata) %in% colnames(counts)))
  print(all(rownames(coldata) == colnames(counts)))
  counts <- counts[, rownames(coldata)]

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~subtype
  )

  dds <- DESeq(dds)
  res <- results(dds)
  print(res[genes_vec, ])

  # extract the normalized data (not stored after dds)
  normalized.data <- counts(
    dds,
    normalized = T
  )
  normalized.data <- t(normalized.data)

  # write the normalized counts to an output file for the heatmap
  write.csv(
    normalized.data,
    file.path(data.dir, paste0(data.type, "_deseq2_normalized.csv")),
    row.names = TRUE,
    quote = FALSE
  )

  normalized.subset <- as.data.frame(normalized.data[, genes_vec])

  normalized.with.meta <- merge(
    normalized.subset,
    coldata,
    by = 0,
    all = TRUE
  )
  colnames(normalized.with.meta)[1] <- "sample"

  # write the normalized data for the genes of interest for the boxplots
  write.csv(
    normalized.with.meta,
    file.path(data.dir, paste0(data.type, "_goi_deseq2_normalized.csv")),
    row.names = FALSE
  )

  # remove intermediate files
  file.remove(mat.path)
  file.remove(col.path)
}

parser <- ArgumentParser(description = "Run DESeq2 normalization")
parser$add_argument("--data_dir", dest = "data_dir", help = "Data directory")
parser$add_argument("--ccc_genes", dest = "ccc_genes", help = "Genes of interest for the CCC analysis")
args <- parser$parse_args()

run_deseq2(
  data.dir = args$data_dir,
  counts.file = "pre_deseq2_pdo_raw_subset.csv",
  meta.file = "pre_deseq2_pdo_meta.csv",
  data.type = "pdo",
  genes = args$ccc_genes
)
