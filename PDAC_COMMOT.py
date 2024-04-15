# This script will take very high memory to run
import os
#import glob 
import pandas as pd
#import shutil
import copy
import csv
import numpy as np
import sys
import altair as alt
from collections import defaultdict
import stlearn as st
import scanpy as sc
import qnorm
import scipy
import pickle
import gzip
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import euclidean_distances
import scanpy as sc
import commot as ct
import gc
import ot
import anndata
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/' , help='The path to dataset') 
parser.add_argument( '--data_name', type=str, default='PDAC_64630', help='The name of dataset')
args = parser.parse_args()

parser = argparse.ArgumentParser()
parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/exp1_C1/outs/' , help='The path to dataset') 
parser.add_argument( '--data_name', type=str, default='PDAC_140694', help='The name of dataset')
args = parser.parse_args()




adata = st.Read10X(path=args.data_path, count_file='filtered_feature_bc_matrix.h5') 
print(adata)

adata.var_names_make_unique()
adata.raw = adata
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type='Secreted Signaling', database='CellChat')
df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata, min_cell_pct=0.05)
ct.tl.spatial_communication(adata, database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)
print('data write')
adata.write(args.data_name+"_commot_adata.h5ad")
print('data read')
adata = sc.read_h5ad(args.data_name+"_commot_adata.h5ad")


ct.tl.communication_direction(adata, database_name='cellchat',  lr_pair=['PLXNB2','MET'],  k=5) #pathway_name='CCL',
print('Plot the CCC')
ct.pl.plot_cell_communication(adata, database_name='cellchat', lr_pair=['CCL19','CCR7'], plot_method='grid', background_legend=True,
    scale=0.00003, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='leiden', cmap='Alphabet',
    normalize_v = True, normalize_v_quantile=0.995, filename='PSAP_cluster.pdf') #pathway_name='CCL', 


