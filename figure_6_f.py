import numpy as np
import csv
import pickle
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.colors import  rgb2hex # LinearSegmentedColormap, to_hex,
from scipy.sparse import csr_matrix
from collections import defaultdict
import pandas as pd
import gzip
import argparse
import os
import scipy.stats
from scipy.sparse.csgraph import connected_components
from pyvis.network import Network
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
import altair as alt
import altairThemes # assuming you have altairThemes.py at your current directoy or your system knows the path of this altairThemes.py.
import gc
import copy
alt.themes.register("publishTheme", altairThemes.publishTheme)
# enable the newly registered theme
alt.themes.enable("publishTheme")


#current_directory = ??

##########################################################
# preprocessDf, plot: these two functions are taken from GW's repository                                                                                                                                                                     /mnt/data0/gw/research/notta_pancreatic_cancer_visium/plots/fatema_signaling/hist.py                                                                                                                                                                                         

def preprocessDf(df):
  """Transform ligand and receptor columns."""
  df["ligand-receptor"] = df["ligand"] + '-' + df["receptor"]
  df["component"] = df["component"] #.astype(str).str.zfill(2)

  return df


def plot(df):
  set1 = altairThemes.get_colour_scheme("Set1", len(df["component"].unique()))
  set1[0] = '#000000'
  base = alt.Chart(df).mark_bar().encode(
            x=alt.X("ligand-receptor:N", axis=alt.Axis(labelAngle=45), sort='-y'),
            y=alt.Y("count()"),
            color=alt.Color("component:N", scale = alt.Scale(range=set1)),
            order=alt.Order("component:N", sort="ascending"),
            tooltip=["component"]
        )
  p = base

  return p

####################### Set the name of the sample you want to visualize ###################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_name', type=str, default='PDAC_64630', help='The name of dataset') # 
#    parser.add_argument( '--model_name', type=str, help='Name of the trained model', required=True)
    parser.add_argument( '--top_edge_count', type=int, default=1300 ,help='Number of the top communications to plot. To plot all insert -1') # 
#    parser.add_argument( '--top_percent', type=int, default=20, help='Top N percentage communications to pick')    
#    parser.add_argument( '--metadata_from', type=str, default='metadata/', help='Path to grab the metadata') 
#    parser.add_argument( '--output_path', type=str, default='output/', help='Path to save the visualization results, e.g., histograms, graph etc.')
    parser.add_argument( '--barcode_info_file', type=str, default='NEST_figures_input/PDAC_64630_barcode_info', help='Path to load the barcode information file produced during data preprocessing step')
    parser.add_argument( '--annotation_file_path', type=str, default='NEST_figures_input/PDAC_64630_annotation.csv', help='Path to load the annotation file in csv format (if available) ')
    parser.add_argument( '--selfloop_info_file', type=str, default='NEST_figures_input/PDAC_64630_self_loop_record', help='Path to load the selfloop information file produced during data preprocessing step')
    parser.add_argument( '--top_ccc_file', type=str, default='NEST_figures_input/PDAC_64630_top20percent.csv', help='Path to load the selected top CCC file produced during data postprocessing step')
    parser.add_argument( '--output_name', type=str, default='NEST_figures_output/', help='Output file name prefix according to user\'s choice')
    args = parser.parse_args()


    output_name = args.output_name
    
    ##################### make cell metadata: barcode_info ###################################
    with gzip.open(args.barcode_info_file, 'rb') as fp:  #b, a:[0:5]        
        barcode_info = pickle.load(fp)    

    ###############################  read which spots have self loops ###############################################################
    with gzip.open(args.selfloop_info_file, 'rb') as fp:  #b, a:[0:5]   _filtered
        self_loop_found = pickle.load(fp)

    ####### load annotations ##############################################
    annotation_data = pd.read_csv(args.annotation_file_path, sep=",")
    pathologist_label=[]
    for i in range (0, len(annotation_data)):
        pathologist_label.append([annotation_data['Barcode'][i], annotation_data['IX_annotation'][i]])

    barcode_type=dict() # record the type (annotation) of each spot (barcode)
    for i in range (0, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1]

    ######################### read the NEST output in csv format ####################################################

    inFile = args.top_ccc_file
    df = pd.read_csv(inFile, sep=",")

    #################################################################################################################
    csv_record = df.values.tolist() # barcode_info[i][0], barcode_info[j][0], ligand, receptor, edge_rank, label, i, j, score

    ## sort the edges based on their rank (column 4), low to high, low being higher attention score
    csv_record = sorted(csv_record, key = lambda x: x[4])
    ## add the column names and take first top_edge_count edges
    # columns are: from_cell, to_cell, ligand_gene, receptor_gene, rank, attention_score, component, from_id, to_id
    df_column_names = list(df.columns)
#    print(df_column_names)

    print(len(csv_record))

    if args.top_edge_count != -1:
        csv_record_final = [df_column_names] + csv_record [0:min(args.top_edge_count, len(csv_record))]

    ## add a dummy row at the end for the convenience of histogram preparation (to keep the color same as altair plot)
    i=0
    j=0
    csv_record_final.append([barcode_info[i][0], barcode_info[j][0], 'no-ligand', 'no-receptor', 0, 0, i, j, 0]) # dummy for histogram

    csv_record = 0
    gc.collect()

    ######################## connected component finding #################################
    print('Finding connected component')
    connecting_edges = np.zeros((len(barcode_info),len(barcode_info)))  
    for k in range (1, len(csv_record_final)-1): # last record is a dummy for histogram preparation
        i = csv_record_final[k][6]
        j = csv_record_final[k][7]
        connecting_edges[i][j]=1
            
    graph = csr_matrix(connecting_edges)
    n_components, labels = connected_components(csgraph=graph,directed=True, connection = 'weak',  return_labels=True) # It assigns each SPOT to a component based on what pair it belongs to
    print('Number of connected components %d'%n_components) 

    count_points_component = np.zeros((n_components))
    for i in range (0, len(labels)):
        count_points_component[labels[i]] = count_points_component[labels[i]] + 1

    id_label = 2 # initially all are zero. =1 those who have self edge but above threshold. >= 2 who belong to some component
    index_dict = dict()
    for i in range (0, count_points_component.shape[0]):
        if count_points_component[i]>1:
            index_dict[i] = id_label
            id_label = id_label+1

    print('Unique component count %d'%id_label)

    for i in range (0, len(barcode_info)):
        if count_points_component[labels[i]] > 1:
            barcode_info[i][3] = index_dict[labels[i]] #2
        elif connecting_edges[i][i] == 1 and (i in self_loop_found and i in self_loop_found[i]): # that is: self_loop_found[i][i] do exist 
            barcode_info[i][3] = 1
        else: 
            barcode_info[i][3] = 0

    # update the label based on found component numbers
    #max opacity
    for record in range (1, len(csv_record_final)-1):
        i = csv_record_final[record][6]
        label = barcode_info[i][3]
        csv_record_final[record][5] = label


    ########################## filtering ###########
    
    ## change the csv_record_final here if you want histogram for specific components/regions only. e.g., if you want to plot only stroma region, or tumor-stroma regions etc.    ##
    #region_of_interest = [...] 
    csv_record_final_temp = []
    csv_record_final_temp.append(csv_record_final[0])
    component_dictionary_dummy = dict()
    for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
        i = csv_record_final[record_idx][6]
        j = csv_record_final[record_idx][7]
        if barcode_type[barcode_info[i][0]] == 'tumor' or barcode_type[barcode_info[j][0]] == 'tumor': 
            csv_record_final_temp.append(csv_record_final[record_idx])
        if csv_record_final[record_idx][5] not in component_dictionary_dummy:
            component_dictionary_dummy[csv_record_final[record_idx][5]] = csv_record_final[record_idx]
            
    # insert just one record from each other components so that the color scheme does not change in the altair scatter plot and histogram :-(
    for component_id in component_dictionary_dummy:
        csv_record_final_temp.append(component_dictionary_dummy[component_id])
    
    csv_record_final_temp.append(csv_record_final[len(csv_record_final)-1])
    csv_record_final = copy.deepcopy(csv_record_final_temp)
    
  
    ###################################  Histogram plotting #################################################################################

    df = pd.DataFrame(csv_record_final)
    df.to_csv('temp_csv.csv', index=False, header=False)
    df = pd.read_csv('temp_csv.csv', sep=",")
    os.remove('temp_csv.csv') # delete the intermediate file

    print('len of loaded csv for histogram generation is %d'%len(df))
    df = preprocessDf(df)
    p = plot(df)
    outPath = output_name + args.data_name + '_histogram_test_tumor.html'
    p.save(outPath)	
    print('Histogram plot generation done')


