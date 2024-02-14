print('package loading')
import numpy as np
import pickle
from scipy import sparse
import numpy as np
import qnorm
from scipy.sparse import csr_matrix
from collections import defaultdict
import pandas as pd
import gzip
import argparse
import os
import scanpy as sc
import altair as alt
import altairThemes # assuming you have altairThemes.py at your current directoy or your system knows the path of this altairThemes.py.
import gc
import copy
alt.themes.register("publishTheme", altairThemes.publishTheme)
# enable the newly registered theme
alt.themes.enable("publishTheme")




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ################## Mandatory ####################################################################
    parser.add_argument( '--data_name', type=str, default='LUAD_TD1', help='The name of dataset') # 
    parser.add_argument( '--data_from', type=str, default='data/LUAD_GSM5702473_TD1/', help='Path to the dataset to read from. Space Ranger outs/ folder is preferred. Otherwise, provide the *.mtx file of the gene expression matrix.')
    parser.add_argument( '--filter_min_cell', type=int, default=1 , help='Minimum number of cells for gene filtering') 
    parser.add_argument( '--tissue_position_file', type=str, default='data/LUAD_GSM5702473_TD1/GSM5702473_TD1_tissue_positions_list.csv', help='If your --data_from argument points to a *.mtx file instead of Space Ranger, then please provide the path to tissue position file.')
    parser.add_argument( '--output_name', type=str, default='NEST_figures_output/', help='Output file name prefix according to user\'s choice')
    args = parser.parse_args()

    if args.tissue_position_file == 'None': # Data is available in Space Ranger output format
        adata_h5 = sc.read_visium(path=args.data_from, count_file='filtered_feature_bc_matrix.h5')
        print('input data read done')
        gene_count_before = len(list(adata_h5.var_names) )    
        sc.pp.filter_genes(adata_h5, min_cells=args.filter_min_cell)
        gene_count_after = len(list(adata_h5.var_names) )  
        print('Gene filtering done. Number of genes reduced from %d to %d'%(gene_count_before, gene_count_after))
        gene_ids = list(adata_h5.var_names)
        coordinates = adata_h5.obsm['spatial']
        cell_barcode = np.array(adata_h5.obs.index)
        cell_vs_gene = sparse.csr_matrix.toarray(adata_h5.X)   

    else:
        
        temp = sc.read_10x_mtx(args.data_from)
        print('*.mtx file read done')
        gene_count_before = len(list(temp.var_names) )
        sc.pp.filter_genes(temp, min_cells=args.filter_min_cell)
        gene_count_after = len(list(temp.var_names) )
        print('Gene filtering done. Number of genes reduced from %d to %d'%(gene_count_before, gene_count_after))
        gene_ids = list(temp.var_names) 
        cell_barcode = np.array(temp.obs.index)
        cell_vs_gene = sparse.csr_matrix.toarray(temp.X)
    
        
        # now read the tissue position file. It has the format:     
        df = pd.read_csv(args.tissue_position_file, sep=",", header=None)   
        tissue_position = df.values
        barcode_vs_xy = dict() # record the x and y coordinates for each spot/cell
        for i in range (0, tissue_position.shape[0]):
            barcode_vs_xy[tissue_position[i][0]] = [tissue_position[i][5], tissue_position[i][4]] #for some weird reason, in the .h5 format for LUAD sample, the x and y are swapped
        
        coordinates = np.zeros((cell_barcode.shape[0], 2)) # insert the coordinates in the order of cell_barcodes
        for i in range (0, cell_barcode.shape[0]):
            coordinates[i,0] = barcode_vs_xy[cell_barcode[i]][0]
            coordinates[i,1] = barcode_vs_xy[cell_barcode[i]][1]
    

    ##################### make metadata: barcode_info ###################################
    i=0
    barcode_info=[]
    for cell_code in cell_barcode:
        barcode_info.append([cell_code, coordinates[i,0],coordinates[i,1], 0]) # last entry will hold the component number later
        i=i+1




    ################### ####################
    target_gene_list = ['APOE', 'SDC1', 'FN1', 'RPSA', 'TGFB1', 'ACVRL1', 'TGFBR2']
    target_gene = ''
    
    for target_gene in target_gene_list:
        target_gene_id = -1
        for j in range (0, len(gene_ids)):
            if gene_ids[j]==target_gene: #'PLXNB2':
                target_gene_id = j
                break
        
        i=0
        barcode_info=[]
        for cell_code in cell_barcode:
            barcode_info.append([cell_code, coordinates[i,0], coordinates[i,1], cell_vs_gene[i,target_gene_id]])
            i=i+1
            
        data_list=dict()
        data_list['X']=[]
        data_list['Y']=[]   
        data_list['gene_expression']=[] 
        
        for i in range (0, len(barcode_info)):
            data_list['X'].append(barcode_info[i][1])
            data_list['Y'].append(-barcode_info[i][2])
            data_list['gene_expression'].append(barcode_info[i][3])
        
        
        source= pd.DataFrame(data_list)
        
        chart = alt.Chart(source).mark_point(filled=True).encode(
            alt.X('X', scale=alt.Scale(zero=False)),
            alt.Y('Y', scale=alt.Scale(zero=False)),
            color=alt.Color('gene_expression:Q', scale=alt.Scale(scheme='magma'))
        )
        chart.save(args.output_name + args.data_name + '_heatmap_' + target_gene + '.html')
