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
    parser.add_argument( '--barcode_info_file', type=str, default='NEST_figures_input_PDAC/PDAC_64630_barcode_info', help='Path to load the barcode information file produced during data preprocessing step')
    parser.add_argument( '--annotation_file_path', type=str, default='NEST_figures_input_PDAC/PDAC_64630_annotation_ayah_morphology.csv', help='Path to load the annotation file in csv format (if available) ') #_ayah_histology
    parser.add_argument( '--output_name', type=str, default='NEST_figures_output/', help='Output file name prefix according to user\'s choice')
    args = parser.parse_args()


    output_name = args.output_name
    ##################### make cell metadata: barcode_info ###################################
    with gzip.open(args.barcode_info_file, 'rb') as fp:  #b, a:[0:5]        
        barcode_info = pickle.load(fp)    

    ####### load annotations ##############################################
    annotation_data = pd.read_csv(args.annotation_file_path, sep=",")
    pathologist_label=[]
    for i in range (0, len(annotation_data)):
        pathologist_label.append([annotation_data['Barcode'][i], annotation_data['tumour morphology'][i]]) 

    barcode_type=dict() # record the type (annotation) of each spot (barcode)
    for i in range (0, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1]


    #### draw altair plot ###########
    data_list=dict()
    data_list['pathology_label']=[]
    data_list['X']=[]
    data_list['Y']=[]     

    for i in range (0, len(barcode_info)):        
        data_list['pathology_label'].append(barcode_type[barcode_info[i][0]])
        data_list['X'].append(barcode_info[i][1])
        data_list['Y'].append(-barcode_info[i][2])

   
    data_list_pd = pd.DataFrame(data_list)
    category_count = len(list(set(data_list['pathology_label']))) 
    set1 = altairThemes.get_colour_scheme("Set1",category_count)
    set1[0] = '#000000'
    chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
        alt.X('X', scale=alt.Scale(zero=False)),
        alt.Y('Y', scale=alt.Scale(zero=False)),
        #shape = alt.Shape('pathology_label:N'), #shape = "pathology_label",
        color=alt.Color('pathology_label:N', scale=alt.Scale(range=set1)), 
        tooltip=['pathology_label'] #,'opacity'
    )

    chart.save(output_name + args.data_name +'_morphology_altair_plot.html')
    print('Altair plot generation done')
    


