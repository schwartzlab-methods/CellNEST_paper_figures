# download human stringdb db and info files 
outdir="../data/protein_database/stringdb"
wget -O "$outdir/9606.protein.physical.links.detailed.v12.0.txt.gz" "https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/9606.protein.physical.links.detailed.v12.0.txt.gz"
gunzip "$outdir/9606.protein.physical.links.detailed.v12.0.txt.gz"

wget -O "$outdir/9606.protein.info.v12.0.txt.gz" "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
gunzip "$outdir/9606.protein.info.v12.0.txt.gz" 
