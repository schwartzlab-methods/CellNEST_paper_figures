library("DESeq2")
library(ggplot2)
library(dplyr)

## Run DESeq2 normalization
run_deseq2 <- function(data.dir,
		      counts.file,
		      meta.file,
                      data.type,
		      genes) {

	mat.path <- file.path(data.dir, counts.file)
	col.path <- file.path(data.dir, meta.file)

	# load counts 
	counts <- as.matrix(read.csv(mat.path, 
				     header = TRUE, 
				     sep = ",", 
				     row.names = "genes"))
        # load metadata
	coldata <- read.csv(col.path, 
			    row.names = 1)

	# compare Classical to Basal-like for analysis 1 and fibroblasts vs. organoids for analysis 2
        if (data.type == "pdo") {
		coldata[["subtype"]] <- factor(coldata[["subtype"]], 
					       levels = c("Basal-like", "Classical"))
	} else if (data.type == "fibro.pdo") {
		coldata[["type"]] <- factor(coldata[["type"]],
					    levels = c("fibroblast", "organoid"))
	}
	
	print(all(rownames(coldata) %in% colnames(counts)))
	print(all(rownames(coldata) == colnames(counts)))
	counts <- counts[, rownames(coldata)]

	if (data.type == "pdo") {
		dds <- DESeqDataSetFromMatrix(countData = counts,
					      colData = coldata,
					      design = ~ subtype)
	} else if (data.type == "fibro.pdo") {
		dds <- DESeqDataSetFromMatrix(countData = counts,
					      colData = coldata,
					      design = ~ type)
	}

	dds <- DESeq(dds)
	res <- results(dds)

	print(res)
	print(res[genes, ])

        # extract the normalized data 
	normalized.data <- counts(dds, 
				  normalized = T)

        normalized.data <- t(normalized.data)

        # write the normalized counts to an output file for the heatmap 
        all.output.file <- paste0(data.type, "_deseq2_normalized.csv")

        write.csv(normalized.data, 
                  file.path(data.dir, all.output.file), 
		  row.names = TRUE, 
		  quote = FALSE)

        subset.matrix <- normalized.data[, genes]
        normalized.df <- as.data.frame(subset.matrix)

        normalized.with.meta <- merge(normalized.df, 
				      coldata, 
				      by = 0, 
				      all = TRUE)

        colnames(normalized.with.meta)[1] <- "sample"

        goi.output.file <- paste0(data.type, "_goi_deseq2_normalized.csv")

	# write the normalized data for the genes of interest for the boxplots 
        write.csv(normalized.with.meta, 
		  file.path(data.dir, goi.output.file), 
		  row.names = FALSE)
}


run_deseq2(
    data.dir = "/mnt/data0/dpaliwal/software/nest-analyses/data", 
    counts.file = "deseq2_pdo_raw_counts.csv",
    meta.file = "pdo_metadata.csv",
    data.type = "pdo",
    genes = c("MET", "PLXNB2")
)

run_deseq2(
    data.dir = "/mnt/data0/dpaliwal/software/nest-analyses/data",
    counts.file = "deseq2_fibro_pdo_raw_counts.csv",
    meta.file = "fibro_pdo_metadata.csv",
    data.type = "fibro.pdo",
    genes = c("MET", "PLXNB2")
)














