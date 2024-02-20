## Tranpose the raw fibroblast & organoid counts .rds file and write to a .csv file

transpose_raw <- function(
	data.dir,
	input.rds, 
	output.csv) {

	organoid <- readRDS(file.path(data.dir, input.rds))
	organoid <- t(organoid)
	write.csv(
		organoid, 
		file = file.path(data.dir, output.csv), 
		row.names = TRUE, 
		quote = FALSE
	)
}

transpose_raw(
	data.dir = "/mnt/data0/dpaliwal/software/nest-analyses/data",
	input.rds = "all_organoids_count_mat_raw_final.93.rds",
	output.csv = "transposed_raw_counts.csv"
)


