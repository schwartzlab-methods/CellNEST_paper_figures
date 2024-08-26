# https://github.com/StatBiomed/SpatialDM/blob/main/tutorial/melanoma.ipynb 
# to do: determine how to select radial basis kernel parameter
import yaml
import os
import spatialdm as sdm
import scanpy as sc
from utilities import *

def runSpatialDM(
        inputPath: str,
        rbkParam: int, # radial basis kernel parameter 
        synthetic: bool, 
        countsFile = None,
        coordsFile = None,
        databasePath = None
    ):
        if synthetic:
            adata = loadSynthetic(inputPath, countsFile, coordsFile)
        else:
            adata = loadBiological(inputPath)

        adata.raw = adata
        # cpm normalization & log-transformation
        sc.pp.normalize_total(adata, target_sum = 1e6)
        sc.pp.log1pm(adata)
        # l: radial basis kernel parameter - decay rate of kernel function 
        # small l: localized signaling, large l: signaling can occur over large distances
        # cutoff: signaling interaction negligible beyond cutoff (200 uM, 1 spot away)
        # alternative to cutoff: n_neighbours = 8
        sdm.weight_matrix(adata, l = rbkParam, cutoff = 0.2, single_cell = False)
        # custom DB for synthetic data 
        if synthetic: 
             sdm.extract_lr(adata, "human", datahost = "nest", databasePath = databasePath)
        # CellChatDB db for biological data 
        # min_cell: filter out sparse ligand/receptor genes 
        else:
            sdm.extract_lr(adata, 'human', min_cell = 3)
        # global Moran selection
        sdm.spatialdm_global(adata, 1000, specified_ind = None, method = 'both', nproc = 1) 
        # select significant pairs
        # threshold: FDR < 0.1 
        sdm.sig_pairs(adata, method = 'permutation', fdr = True, threshold = 0.1) 
        # local spot selection significant pairs to identify where lri occurs at single-spot resolution
        sdm.spatialdm_local(adata, n_perm = 1000, method = 'both', specified_ind = None, nproc = 1)  
        sdm.sig_spots(adata, method = 'permutation', fdr = False, threshold = 0.1)
        ## access global results 
        # adata.uns["global_res"]
        ## access local results 
        # adata.uns["local_perm_p"] 
        ## global and local results are easily accessible through global_res and local_perm_p or local_z_p
        ## save in list of list format 
        ## [ [from_cell, to_cell, lr name, strength], …, [from_cell, to_cell, lr name, strength]] 

def main():
    with open("config.yml", "r") as f:
         config = yaml.safe_load(f)
    databasePath = config["files"]["nest_database"]

    # synthetic 
    runSpatialDM(
         inputPath = os.path.join(config["directories"]["synthetic"]["uniform_distribution"], "no_noise")
         countsFile = "uniform_distribution_cellvsgene_not_normalized",
         coordsFile = "uniform_distribution_coordinate",
         # random - need to find how to specify 
         rbkParam = 75, 
         synthetic = True, 
         databasePath = databasePath
    )
    # lymph node 
    # spot-to-spot distance: 100 um; need to account for scaling? 
    runSpatialDM(
         inputPath = config["directories"]["lymph_node_visium"], 
         # random - need to find how to specify 
         rbkParam = 75, 
         synthetic = False
    )

if __name__ == "__main__":
    main()