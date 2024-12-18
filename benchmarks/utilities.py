import os
import numpy as np
import pandas as pd
import gzip
import pickle
import anndata as ad 
import scanpy as sc
import pandas as pd
import scipy
import numpy as np
import glob 

def loadSynthetic(
    inputPath: str,
    countsFilePattern: str,
    coordsFilePattern: str
) -> ad.AnnData:
    countsFile = glob.glob(os.path.join(inputPath, countsFilePattern))[0]
    coordsFile = glob.glob(os.path.join(inputPath, coordsFilePattern))[0]
    with gzip.open(os.path.join(countsFile), "rb") as fp:
        countsMat = pickle.load(fp)
    with gzip.open(os.path.join(coordsFile), "rb") as fp:
        tempX, tempY, _ = pickle.load(fp)
    nCells = tempX.shape[0]
    coords = np.zeros((nCells, 2))
    for i in range(0, nCells):
        coords[i][0] = tempX[i]
        coords[i][1] = tempY[i] 
    print(coords.shape) # (nCells, 2)
    adata = ad.AnnData(X = countsMat)
    print(inputPath)
    min_value = adata.X.min()
    print("size of adata:", adata.X.size)
    num_negatives = (adata.X < 0).sum()
    num_positives = (adata.X > 0).sum()
    print("number of negatives:", num_negatives)
    print("number of positives:", num_positives)
    if min_value < 0:
        adata.X +=(abs(min_value))
    num_negatives = (adata.X < 0).sum()
    num_positives = (adata.X > 0).sum()
    print("number of negatives:", num_negatives)
    print("number of positives:", num_positives)
    print(coords)
    adata.obsm["spatial"] = coords
    sc.pp.filter_cells(adata, min_counts = 1)
    print(f"number of cells after filtering: {adata.n_obs}")
    adata.var_names_make_unique()

    return adata

def loadBiological(
    inputPath: str
) -> ad.AnnData:
    adata = sc.read_visium(inputPath)
    adata.var_names_make_unique()

    return adata

def filterDb(
    dbPath: str
) -> None:
    df = pd.read_csv(dbPath)
    df = df.applymap(lambda x: x.replace('g', '') if isinstance(x, str) and 'g' in x else x)
    df.to_csv(dbPath, index = False)

# check maximum distance between 
# ligand-receptor pairs in the ground truth dictionary 
def checkDistance(dataDir):
    coords = pd.read_csv(os.path.join(dataDir, "coords.csv"), index_col = 0)
    groundTruth = glob.glob(os.path.join(dataDir, "*ground_truth_ccc"))[0]
    with gzip.open(groundTruth, 'rb') as fp:  
        _, TP, random_activation = pickle.load(fp)
    max_dist = 0
    for i, j in TP.items():
        for k in j:
            if k in coords.index:
                i_coords = coords.loc[i].values
                k_coords = coords.loc[k].values 
                distance = ((i_coords[0] - k_coords[0]) ** 2 + (i_coords[1] - k_coords[1]) ** 2) ** 0.5
                max_dist = max(max_dist, distance)
    return max_dist 

def writeSpatial(adata: ad.AnnData, dataDir: str) -> None:
    adata.var_names_make_unique()
    counts = adata.X
    if isinstance(counts, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)):
        counts = counts.toarray()
    countsDf = pd.DataFrame(data = counts, index = adata.obs_names, columns = adata.var_names)
    countsDf.to_csv(os.path.join(dataDir, "counts.csv"))
    coords = adata.obsm["spatial"]
    coordsDf = pd.DataFrame(data = coords, index = adata.obs_names, columns = ["x", "y"])
    coordsDf.to_csv(os.path.join(dataDir, "coords.csv"))

def leidenCluster(adata: ad.AnnData, dataDir: str, res: float) -> None:
    sc.pp.normalize_total(adata, target_sum = 1e6)
    sc.pp.log1p(adata)
    # unit variance & zero mean across all cells 
    sc.pp.scale(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors = 15, n_pcs = 50)
    sc.tl.leiden(adata, resolution = res)
    clusterDf = pd.DataFrame(adata.obs["leiden"])
    clusterDf.reset_index(inplace = True)
    clusterDf.rename(columns = {"index":"spot", "leiden":"cluster"}, inplace = True)
    clusterDf.to_csv(os.path.join(dataDir, "leiden.csv"), index = False)

def getRawCounts(adata: ad.AnnData):
    if "raw_counts" in adata.layers:
        mat = adata.layers["raw_counts"]
    else:
        mat = adata.X
    counts = mat.todense() if hasattr(mat, 'todense') else mat

    return counts

# spatial coordinates & expression
def cytospaceSpatial(adata: ad.AnnData, outFolder: str) -> None:
    adata.var_names_make_unique()
    counts = getRawCounts(adata)
    print(counts.shape)
    counts_df = pd.DataFrame(
        counts.T,
        index = adata.var_names,
        columns = adata.obs_names
    )
    counts_df.index.name = 'GENES'
    os.makedirs(outFolder, exist_ok = True)
    counts_df.to_csv(os.path.join(outFolder, f"spatial_counts.csv"))
    # coords to csv 
    coordinates = adata.obsm['spatial']
    coordinates_df = pd.DataFrame(coordinates, index = adata.obs_names, columns=['row', 'col'])
    coordinates_df.reset_index(inplace = True) 
    coordinates_df.rename(columns={'index': 'SpotID'}, inplace = True)
    coordinates_df.to_csv(os.path.join(outFolder, f"spatial_coordinates.csv"), index = False)

def cytospaceRef(adata: ad.AnnData, outFolder: str, ct_col: str) -> None:
    adata.var_names_make_unique()
    # subset 10 000 cells with an equal distribution across all cell types
    total_cells = 10000
    cell_types = adata.obs[ct_col].value_counts()
    num_cell_types = len(cell_types)
    cells_per_type = total_cells // num_cell_types
    sampled_indices = []
    for cell_type, count in cell_types.items():
        cell_type_indices = adata.obs[adata.obs[ct_col] == cell_type].index
        if count <= cells_per_type:
            sampled_indices.extend(cell_type_indices)
        else:
            sampled_indices.extend(np.random.choice(cell_type_indices, cells_per_type, replace = False))
    adata_subset = adata[sampled_indices, :].copy()
    # raw counts to csv 
    counts = getRawCounts(adata_subset)
    counts_df = pd.DataFrame(counts.T, index = adata_subset.var_names, columns = adata_subset.obs_names)
    adata_subset.var_names_make_unique()
    counts_df.reset_index(inplace = True)
    counts_df.rename(columns={'index': 'GENES'}, inplace = True)
    os.makedirs(outFolder, exist_ok = True)
    counts_df.to_csv(os.path.join(outFolder, f"sc_counts.csv"), index = False)
    # write cell type labels to csv 
    cell_type_df = pd.DataFrame({
    'Cell IDs': adata_subset.obs_names,
    'CellType': adata_subset.obs[ct_col]
    })
    # remove special chars
    cell_type_df['CellType'] = cell_type_df['CellType'].str.replace(r'[^a-zA-Z0-9 ]', ' ', regex = True)
    cell_type_df.to_csv(os.path.join(outFolder, f"sc_cell_types.csv"), index = False)
