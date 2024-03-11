library(argparse)

transpose_raw <- function(data.dir,
                          input.rds,
                          output.csv) {
  obj <- readRDS(file.path(data.dir, input.rds))
  obj <- t(obj)

  write.csv(
    obj,
    file = file.path(data.dir, output.csv),
    row.names = TRUE,
    quote = FALSE
  )
}

parser <- ArgumentParser(description = "Transpose the raw fibroblast & organoid counts .rds file and write to a .csv file")
parser$add_argument("--data_dir", dest = "data_dir", help = "Data directory containing .rds input file")
parser$add_argument("--rds", dest = "rds", help = "Input .rds file containing raw counts for organoids and fibroblasts")
parser$add_argument("--output_csv", dest = "output_csv", help = "Output csv file to write the transposed counts to")
args <- parser$parse_args()

transpose_raw(
  data.dir = args$data_dir,
  input.rds = args$rds,
  output.csv = args$output_csv
)
