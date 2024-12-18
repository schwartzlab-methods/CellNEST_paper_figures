outdir="../data/protein_database/stringdb"
mkdir -p "$outdir"
cd "$outdir"

# human
wget "https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/9606.protein.physical.links.detailed.v12.0.txt.gz"
gunzip "9606.protein.physical.links.detailed.v12.0.txt.gz"
wget "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
gunzip "9606.protein.info.v12.0.txt.gz"

# mouse
wget "https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/10090.protein.physical.links.detailed.v12.0.txt.gz"
gunzip "10090.protein.physical.links.detailed.v12.0.txt.gz"
wget "https://stringdb-downloads.org/download/protein.info.v12.0/10090.protein.info.v12.0.txt.gz"
gunzip "10090.protein.info.v12.0.txt.gz"
