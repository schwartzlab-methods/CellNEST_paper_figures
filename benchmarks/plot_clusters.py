from scipy import sparse
import pickle
import scipy.linalg
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_hex, rgb2hex
from typing import List
import altair as alt
from vega_datasets import data
import pandas as pd
import numpy as np
import csv
import pickle
import scipy.io as sio
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import os
import yaml
import altairThemes
alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")

# specific to leiden = 0.6 
colours = [
    "brown",    # Cluster 0
    "red",      # Cluster 1
    "blue",     # Cluster 2
    "green",    # Cluster 3
    "purple",   # Cluster 4 (dark purple)
    "orange",   # Cluster 5
    "pink",     # Cluster 6
    "violet"    # Cluster 7 (light purple)
]


def plotClusters(
        spacerangerDir: str,
        # generatedDataDir: str,
        # minCells: int,
        # minCounts: int,
        # nCompsPCA = 200,
        labelFile: str,
        savePath: str,
    ):
    # os.makedirs(generatedDataDir, exist_ok = True)
    # adata_h5 = st.Read10X(spacerangerDir, count_file = 'filtered_feature_bc_matrix.h5') 
    adata_h5 = sc.read_visium(spacerangerDir, count_file = 'filtered_feature_bc_matrix.h5') 
    # sc.pp.filter_cells(adata_h5, min_counts = minCounts)
    # gene_ids = adata_h5.var['gene_ids']
    coordinates = adata_h5.obsm['spatial']
    cell_barcode = np.array(adata_h5.obs.index)
    # sc.pp.filter_genes(adata_h5, min_cells = minCells)
    # adata_X = sc.pp.normalize_total(adata_h5, target_sum = 1, exclude_highly_expressed = True, inplace = False)['X']
    # adata_X = sc.pp.scale(adata_X)
    # if nCompsPCA > 0:
        # adata_X = sc.pp.pca(adata_X, n_comps = nCompsPCA)
    # features = adata_X
    # with open(generatedDataDir + 'features', 'wb') as fp:
        # pickle.dump(features, fp)
    # np.save(generatedDataDir + 'features.npy', features)
    coordinates = np.array(coordinates)
    cell_barcode = np.array(cell_barcode)
    # np.save(generatedDataDir + 'coordinates.npy', coordinates)   
    # np.save(generatedDataDir + 'barcodes.npy', cell_barcode)

    barcode_info=[]
    for i in range (0, coordinates.shape[0]):
        barcode_info.append([cell_barcode[i], coordinates[i,0],coordinates[i,1],0])

    for i in range(len(barcode_info)):
        barcode_info[i][1], barcode_info[i][2] = (
            barcode_info[i][1],
            - barcode_info[i][2],
        )

    number = 20
    cmap = plt.get_cmap('tab20')
    colors = [cmap(i) for i in np.linspace(0, 1, number)]

    number = 20
    cmap = plt.get_cmap('tab20b')
    colors_2 = [cmap(i) for i in np.linspace(0, 1, number)]

    colors=colors+colors_2

    number = 20
    cmap = plt.get_cmap('tab20c')
    colors_2 = [cmap(i) for i in np.linspace(0, 1, number)]

    colors=colors+colors_2

    number = 8
    cmap = plt.get_cmap('Set2')
    colors_2 = [cmap(i) for i in np.linspace(0, 1, number)]

    colors=colors+colors_2

    number = 12
    cmap = plt.get_cmap('Set3')
    colors_2 = [cmap(i) for i in np.linspace(0, 1, number)]

    colors=colors+colors_2

    for i in range (0, len(colors)): 
        colors[i] = matplotlib.colors.to_hex([colors[i][0], colors[i][1], colors[i][2], colors[i][3]])
    cell_label=[]
    with open(labelFile) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            cell_label.append(line) 
    barcode_label=dict()
    for i in range (1, len(cell_label)):
        if len(cell_label[i])>0 :
            barcode_label[cell_label[i][0]] = cell_label[i][1]
    for i in range (0, len(barcode_info)):
        if barcode_info[i][0] in barcode_label:
            barcode_info[i][3] = barcode_label[barcode_info[i][0]]

    data_list=dict()
    data_list['cluster_label']=[]
    data_list['X']=[]
    data_list['Y']=[]

    for i in range (0, len(barcode_info)):
        data_list['cluster_label'].append(barcode_info[i][3])
        data_list['X'].append(barcode_info[i][1])
        data_list['Y'].append(barcode_info[i][2])

    data_list_pd = pd.DataFrame(data_list)

    chart = alt.Chart(data_list_pd).mark_point(filled=True, size = 20).encode(
        alt.X('X', scale=alt.Scale(zero=False)),
        alt.Y('Y', scale=alt.Scale(zero=False)),
        tooltip=["cluster_label"],
        #alt.Size('pop:Q'),
        color=alt.Color('cluster_label:N', title = "cluster", scale = alt.Scale(range = colours))
    ).configure_legend()

    print('output is saved here: '+savePath)
    chart.save(savePath)

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    dataDir = os.path.join(config["directories"]["data"], "lymph_node")
    plotClusters(
        spacerangerDir = dataDir,
        labelFile = os.path.join(dataDir, "leiden.csv"), 
        savePath = os.path.join(dataDir, "leiden_clusters.html")
    )

if __name__ == "__main__":
    main()


