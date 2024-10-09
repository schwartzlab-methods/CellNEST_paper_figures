output_dir="../data/protein_database/NicheNet_V2_raw"
wget https://zenodo.org/records/8016880/files/NicheNet_V2.zip?download=1 -O NicheNet_V2.zip
unzip NicheNet_V2.zip -d "$output_dir"

data_source_xlsx="https://github.com/saeyslab/nichenetr/blob/master/vignettes/data_sources.xlsx"
wget -nc -O "$output_dir/data_sources.xlsx" "$data_source_xlsx"
# converted to csv on desktop 



