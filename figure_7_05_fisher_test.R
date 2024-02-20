library(coin)

## Perform the Exact Two-Sample Fisher-Pitman Permutation Test to compare DESeq2 normalized bulk RNA-seq data 
exact_fisher <- function(
	data.dir,
	file.name,
        data.type) {

	df <- read.csv(
		file.path(data.dir, file.name), 
		header = TRUE
	)

        # exclude the myc-amplified outlier 
        outlier <- "PDA_107407_Lv_M_OBp9M3_pr3"
        df <- df[df$sample != outlier, ]

        # analyze MET expression between Classical and Basal-like PDOs 
        if (data.type == "pdo") {
                df$subtype <- as.factor(df$subtype)
		met_res <- oneway_test(MET ~ subtype, 
			    data = df, 
			    distribution = "exact")
		print(met_res)

        # Analyze MET and PLXNB2 expression between fibroblasts and PDOs
	} else if (data.type == "fibro.pdo") {
                df$type <- as.factor(df$type)
		met_res <- oneway_test(MET ~ type,
			    data = df,
			    distribution = "exact")
		print(met_res)
		plxnb2_res <- oneway_test(PLXNB2 ~ type,
			    data = df,
			    distribution = "exact")
		print(plxnb2_res)
	}	
}

# Compute the p-value for the Classical vs. Basal-like PDO analysis
exact_fisher(data.dir = "/mnt/data0/dpaliwal/software/nest-analyses/data", 
	     file.name = "pdo_goi_deseq2_normalized.csv", 
	     data.type = "pdo")

# Compute the p-value for the fibroblast vs. PDO analysis
exact_fisher(data.dir = "/mnt/data0/dpaliwal/software/nest-analyses/data",
	     file.name = "fibro.pdo_goi_deseq2_normalized.csv", 
	     data.type = "fibro.pdo")



