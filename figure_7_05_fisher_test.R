library(coin)
library(argparse)

exact_fisher <- function(data.dir,
                         file.name,
                         outlier,
                         genes) {
  # convert genes string to column vector
  genes_vec <- unlist(strsplit(genes, ","))
  df <- read.csv(
    file.path(data.dir, file.name),
    header = TRUE
  )
  # exclude the myc-amplified outlier
  df <- df[df$sample != outlier, ]
  df$subtype <- as.factor(df$subtype)
  df$type <- as.factor(df$type)

  for (gene in genes_vec) {
    # analyze MET expression between Classical and Basal-like PDOs
    formula <- as.formula(paste(gene, "~ subtype"))

    res <- oneway_test(
      formula,
      data = df,
      distribution = "exact"
    )
    print(res)
  }

  file.remove(file.path(data.dir, file.name))
}

parser <- ArgumentParser(description = "Perform the Exact Two-Sample Fisher-Pitman Permutation Test to compare DESeq2 normalized bulk RNA-seq data")
parser$add_argument("--data_dir", dest = "data_dir", help = "Data directory")
parser$add_argument("--outlier_id", dest = "outlier_id", help = "Sample ID of outlier")
parser$add_argument("--ccc_genes", dest = "ccc_genes", help = "CCC genes of interest")
args <- parser$parse_args()

# Compute the p-value for the Classical vs. Basal-like PDO analysis
exact_fisher(
  data.dir = args$data_dir,
  file.name = "pdo_goi_deseq2_normalized.csv",
  outlier = args$outlier_id,
  genes = args$ccc_genes
)
