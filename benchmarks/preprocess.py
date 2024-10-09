from utilities import *
import yaml 
import os
import scanpy as sc
import glob

def main():
    with open("config.yml", "r") as f:
         config = yaml.safe_load(f)
    datasets = config["params"]["datasets"]
    pdac_samples = config["params"]["pdac_samples"]
    cytospace_dir = os.path.join(config["directories"]["data"], "cytospace")
    os.makedirs(cytospace_dir, exist_ok = True)
    for dataset in datasets:
        data_dir = os.path.join(config["directories"]["data"], dataset)
        print(f"Searching in directory: {data_dir}")
        # remove gs from synthetic database files 
        database = glob.glob(os.path.join(data_dir, "*lrdb.csv"))[0]
        if os.path.exists(database):
            filterDb(database)
        if dataset == "lymph_node":
            adata = sc.read_visium(data_dir)
            sc.pp.filter_genes(adata, min_cells = 1)
            cytospaceSpatial(adata, outFolder = os.path.join(cytospace_dir, dataset))
            adata_sc = sc.read(f'./data/adata_sc_lymph.h5ad', backup_url='https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad')
            cytospaceRef(adata_sc, outFolder = os.path.join(cytospace_dir, dataset), ct_col = "Subset")
        else:
            adata = loadSynthetic(inputPath = data_dir, countsFilePattern = "*cellvsgene_not_normalized", coordsFilePattern = "*coordinate")
        adata.write_h5ad(os.path.join(data_dir, "adata.h5ad"))
        writeSpatial(adata = adata, dataDir = data_dir)
        leidenCluster(adata, data_dir)

    # pre-process pdac samples for cytospace
    # for sample in pdac_samples:
        # pdac_spatial_adata = sc.read_visium(os.path.join(config["directories"]["pdac_data"], "spaceranger", sample, "outs"))
        # pdac_spatial_adata.var_names_make_unique()
        # cytospaceSpatial(pdac_spatial_adata, outFolder = os.path.join(cytospace_dir, sample))
        # adata_sc = sc.read(os.path.join(config["directories"]["pdac_data"], "Moffit_Combined.h5ad"))
        # adata_sc.var_names_make_unique()
        # cytospaceRef(adata_sc, outFolder = os.path.join(cytospace_dir, sample), ct_col = "max_moffitt")

if __name__ == "__main__":
    main()
