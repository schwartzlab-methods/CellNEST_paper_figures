# Set directories 
data_dir="../data"
fig_dir="../figure_7_plots"
mkdir -p "$fig_dir"

# Set counts & metadata files provided 
input_rds="fibroblast_pdo_raw_counts.rds"
input_pdo_meta_file="received_pdo_metadata.csv"
input_fibro_meta_file="received_fibroblast_pdo_metadata.csv"

# Sample IDs of pure classical-like and basal PDOs 
pure_pdo_ids=("PDA_101053_Lv_M_OBp11M3_pr3" "PDA_097484_Lv_M_OBp10M3_pr3" "PDA_093429_Lv_M_OBp11M3_pr3" "PDA_091168_Lv_M_OBp11M3_pr3" "PDA_102454_Lv_M_OBp7_pr3" "PDA_106651_Lv_M_OBp9M3_pr3" "PDA_104760_Lv_M_OBp14M3_pr3" "PDA_100809_Lv_M_OBp8M3_pr3" "PDA_107872_Lv_M_OBp6M3_pr3" "PDA_101781_Lv_M_OBp7M3_pr3" "PDA_091416_Lv_M_OBp11M3_pr3" "PDA_091416_Lv_M_OBp11M3_pr3" "PDA_107407_Lv_M_OBp9M3_pr3" "PDA_097768_Lv_M_OBp11M3_pr3")
# CCC genes of interest for gene expression analyses 
ccc_genes=("MET" "LGALS3")
# Myc amplified outlier 
outlier_id="PDA_107407_Lv_M_OBp9M3_pr3"

# 1 - Tranpose the raw fibroblast & organoid counts, and write to a .csv file in the same directory 
Rscript ./figure_7_01_rds_to_csv.R \
    --data_dir "$data_dir" \
    --rds "$input_rds" \
    --output_csv "transposed_raw_counts.csv"

# 2 - Pre-process counts & metadata files for DESEQ2 
python ./figure_7_02_deseq2_preprocess.py \
    --data_dir "$data_dir" \
    --input_pdo_meta_file "$input_pdo_meta_file" \
    --pure_ids  $(IFS=, ; echo "${pure_pdo_ids[*]}") \
    --input_counts_file "transposed_raw_counts.csv" 

# 3 - Run DESeq2 normalization 
Rscript ./figure_7_03_deseq2.R \
    --data_dir "$data_dir" \
    --ccc_genes $(IFS=, ; echo "${ccc_genes[*]}")

# 4 - Create boxplots 
python ./figure_7_04_box_plot.py \
    --data_dir "$data_dir" \
    --outlier_id "$outlier_id" \
    --ccc_genes $(IFS=, ; echo "${ccc_genes[*]}") \
    --fig_dir "$fig_dir" 

# 5 - Perform Exact Fisher-Pitman Permutation test 
Rscript ./figure_7_05_fisher_test.R \
    --data_dir "$data_dir" \
    --ccc_genes $(IFS=, ; echo "${ccc_genes[*]}") \
    --outlier_id "$outlier_id" 

# 6 - Create a boxplot containing MET, and classical-like and basal associated genes defined by Moffit 2023 
python ./figure_7_06_heatmap.py \
    --data_dir "$data_dir" \
    --fig_dir "$fig_dir" \
    --input_pdo_meta_file "$input_pdo_meta_file" \
    --outlier_id "$outlier_id" \




