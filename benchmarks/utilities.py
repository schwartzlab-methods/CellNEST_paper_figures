import os
import numpy as np
import pandas as pd
import gzip
import pickle
import anndata as ad 
import scanpy as sc

def loadSynthetic(
    inputPath: str,
    countsFile: str,
    coordsFile: str
) -> ad.AnnData:
    with gzip.open(os.path.join(inputPath, countsFile), "rb") as fp:
        countsMat = pickle.load(fp)
    with gzip.open(os.path.join(inputPath, coordsFile), "rb") as fp:
        tempX, tempY, _ = pickle.load(fp)
    nCells = tempX.shape[0]
    coords = np.zeros((nCells, 2))
    for i in range(0, nCells):
        coords[i][0] = tempX[i]
        coords[i][1] = tempY[i] 
    print(coords.shape()) # (nCells, 2)
    adata = ad.AnnData(X = countsMat)
    adata.obsm["spatial"] = coords

    return(adata)

def loadBiological(
    inputPath: str
) -> ad.AnnData:
    adata = sc.read_visium(inputPath)
    
    return (adata)
