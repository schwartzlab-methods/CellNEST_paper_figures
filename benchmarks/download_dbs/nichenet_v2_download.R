library(dplyr)
library(yaml)

config <- yaml.load_file("../config.yml")
input_dir <- file.path(paste0(".", config$directories$data), "protein_database", "NicheNet_V2_raw", "NicheNet_V2")
signaling_dir <- file.path(input_dir, "networks", "data", "signaling")
gr_dir <- file.path(input_dir, "networks", "data", "gene_regulatory")
source_weights_dir <- file.path(input_dir, "evaluation", "optimization", "results")

output_dir <- file.path(paste0(".", config$directories$data), "protein_database", "nichenet_v2")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# signaling & gene regulatory network databases 
database_files <- list(
    file.path(signaling_dir, "signaling_network_human_21122021.rds"),
    file.path(signaling_dir, "signaling_network_mouse_21122021.rds"),
    file.path(gr_dir, "gr_network_human_21122021.rds"),
    file.path(gr_dir, "gr_network_mouse_21122021.rds")
)

for (file in database_files) {
    db_obj <- readRDS(file)
    db_obj_updated <- db_obj %>% 
    group_by(from, to) %>%
    summarise(
        source = paste(unique(source), collapse = ","),
        database = paste(unique(database), collapse = ",")
    ) %>%
    ungroup()
    csv_filename <- sub("\\.rds$", ".csv", basename(file))
    write.csv(
        db_obj_updated, 
        file.path(output_dir, csv_filename), 
        row.names = F
    )
}

# source weights file 
source_weights_fname <- "all_sourceweights_top25_summarized_final.rds"
source_weights_csv <- sub("\\.rds$", ".csv", basename(source_weights_fname))
weights_obj <- readRDS(file.path(source_weights_dir, source_weights_fname))
write.csv(
    weights_obj, 
    file.path(output_dir, source_weights_csv), 
    row.names = F
)
