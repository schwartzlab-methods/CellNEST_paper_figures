from utilities import *
import yaml 
import os
import scanpy as sc
import glob


def main():
    with open("config.yml", "r") as f:
         config = yaml.safe_load(f)
    datasets = config["params"]["datasets"]
    # datasets = [dataset for dataset in datasets if "mechanistic_random_ccc_wo_relay" in dataset]
    datasets.append("lung") # don't add to config
    pdac_samples = config["params"]["pdac_samples"]
    cytospace_dir = os.path.join(config["directories"]["data"], "cytospace")
    os.makedirs(cytospace_dir, exist_ok = True)
    for dataset in datasets:
        data_dir = os.path.join(config["directories"]["data"], dataset)
        print(f"Searching in directory: {data_dir}")
        # remove gs from synthetic database files 
        database_file = glob.glob(os.path.join(data_dir, "*lrdb.csv"))
        if database_file:
            db = database_file[0]
            filterDb(db)
        if dataset == "lymph_node":
            adata = sc.read_visium(data_dir)
            if adata.raw is not None:
                adata.layers["raw_counts"] = adata.raw.X.copy()
            sc.pp.filter_genes(adata, min_cells = 1)
            print("done")
            cytospaceSpatial(adata, outFolder = os.path.join(cytospace_dir, dataset))
            adata_sc = sc.read(f'./data/adata_sc_lymph.h5ad', backup_url='https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad')
            cytospaceRef(adata_sc, outFolder = os.path.join(cytospace_dir, dataset), ct_col = "Subset")
            leiden_res = 0.6
        elif dataset == "lung":
            adata = sc.read_10x_mtx(data_dir)
            positions = pd.read_csv(glob.glob(os.path.join(data_dir, "*tissue_positions_list.csv"))[0], header = None)
            positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
            barcodes = adata.obs_names.tolist() 
            positions_filtered = positions[positions["barcode"].isin(barcodes)]
            positions_filtered = positions_filtered.set_index('barcode').reindex(barcodes).reset_index()
            adata.obsm['spatial'] = positions_filtered[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
       	    if adata.raw is not	None:
       	       	adata.layers["raw_counts"] = adata.raw.X.copy()
            sc.pp.filter_genes(adata, min_cells = 1)
            cytospaceSpatial(adata, outFolder = os.path.join(cytospace_dir, dataset))
            adata_sc = sc.read(os.path.join(data_dir, "single_cell_39114284.h5ad"))
            cytospaceRef(adata_sc, outFolder = os.path.join(cytospace_dir, dataset), ct_col = "Anno1")
            leiden_res = 1
        else:
            adata = loadSynthetic(inputPath = data_dir, countsFilePattern = "*not_normalized", coordsFilePattern = "*coordinate")
            leiden_res = 1
        adata.write_h5ad(os.path.join(data_dir, "adata.h5ad"))
        writeSpatial(adata = adata, dataDir = data_dir)
        print("running leiden")
        leidenCluster(adata, data_dir, res = leiden_res)
        checkDistance(data_dir)

    # pre-process pdac samples for cytospace
    for sample in pdac_samples:
        pdac_spatial_adata = sc.read_visium(os.path.join(config["directories"]["pdac_data"], "spaceranger", sample, "outs"))
        pdac_spatial_adata.var_names_make_unique()
        cytospaceSpatial(pdac_spatial_adata, outFolder = os.path.join(cytospace_dir, sample))
        adata_sc = sc.read(os.path.join(config["directories"]["pdac_data"], "Moffit_Combined.h5ad"))
        adata_sc.var_names_make_unique()
        cytospaceRef(adata_sc, outFolder = os.path.join(cytospace_dir, sample), ct_col = "max_moffitt")

if __name__ == "__main__":
    main()
