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
    parser.add_argument( '--data_name', type=str, default='V1_Human_Lymph_Node_spatial', help='The name of dataset') # 
    parser.add_argument( '--top_edge_count', type=int, default= -1,help='Number of the top communications to plot. To plot all insert -1') # 
    parser.add_argument( '--barcode_info_file', type=str, default='NEST_figures_input_human_lymph/V1_Human_Lymph_Node_spatial_barcode_info', help='Path to load the barcode information file produced during data preprocessing step')
    parser.add_argument( '--annotation_file_path', type=str, default='NEST_figures_input_human_lymph/V1_Human_Lymph_Node_spatial_annotation.csv', help='Path to load the annotation file in csv format (if available) ')
    parser.add_argument( '--selfloop_info_file', type=str, default='NEST_figures_input_human_lymph/V1_Human_Lymph_Node_spatial_self_loop_record', help='Path to load the selfloop information file produced during data preprocessing step')
    parser.add_argument( '--top_ccc_file', type=str, default='NEST_figures_input_human_lymph/V1_Human_Lymph_Node_spatial_top20percent.csv', help='Path to load the selected top CCC file produced during data postprocessing step')
    parser.add_argument( '--input_edge_file', type=str, default='NEST_figures_input_human_lymph/V1_Human_Lymph_Node_spatial_input_graph.csv', help='Input edge list file name')
    parser.add_argument( '--output_name', type=str, default='NEST_figures_output/', help='Output file name prefix according to user\'s choice')
    parser.add_argument( '--histogram_attention_score', type=int, default=-1, help='Set 1 to plot the histograms based on total attention scores of the ligand-receptor pairs')
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
        pathologist_label.append([annotation_data['Barcode'][i], annotation_data['Type'][i]])

    barcode_type=dict() # record the type (annotation) of each spot (barcode)
    for i in range (0, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1]

    ######################### read the NEST output in csv format ####################################################

    inFile = args.top_ccc_file
    df = pd.read_csv(inFile, sep=",")

    #################################################################################################################
    csv_record = df.values.tolist() # barcode_info[i][0], barcode_info[j][0], ligand, receptor, edge_rank, label, i, j, score
    
    ##########################################################################
    ## sort the edges based on their rank (column 4), low to high, low being higher attention score
    csv_record = sorted(csv_record, key = lambda x: x[4])
    ## add the column names and take first top_edge_count edges
    # columns are: from_cell, to_cell, ligand_gene, receptor_gene, rank, component, from_id, to_id, attention_score
    df_column_names = list(df.columns)
#    print(df_column_names)

    print(len(csv_record))

    csv_record_final = [df_column_names] + csv_record #[0:(len(csv_record)*20)//100] 
      
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
        if barcode_type[barcode_info[i][0]] == 'T-cell' and barcode_type[barcode_info[j][0]] == 'T-cell': 
            csv_record_final_temp.append(csv_record_final[record_idx])
        if csv_record_final[record_idx][5] not in component_dictionary_dummy:
            component_dictionary_dummy[csv_record_final[record_idx][5]] = csv_record_final[record_idx]
            
    # insert just one record from each other components so that the color scheme does not change in the altair scatter plot and histogram :-(
    for component_id in component_dictionary_dummy:
        csv_record_final_temp.append(component_dictionary_dummy[component_id])
    
    csv_record_final_temp.append(csv_record_final[len(csv_record_final)-1])
    csv_record_final = copy.deepcopy(csv_record_final_temp)
    ###################################################  
    component_list = dict()
    for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
        record = csv_record_final[record_idx]
        i = record[6]
        j = record[7]
        component_label = record[5]
        barcode_info[i][3] = component_label #?
        barcode_info[j][3] = component_label #?
        component_list[component_label] = ''

    component_list[0] = ''
    unique_component_count = max(len(component_list.keys()), id_label)


    ##################################### Altair Plot ##################################################################
    ## dictionary of those spots who are participating in CCC ##
    active_spot = defaultdict(list)
    for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
        record = csv_record_final[record_idx]
        i = record[6]
        pathology_label = barcode_type[barcode_info[i][0]]
        component_label = record[5]
        X = barcode_info[i][1]
        Y = -barcode_info[i][2]
        opacity = np.float(record[8])
        active_spot[i].append([pathology_label, component_label, X, Y, opacity])
        
        j = record[7]
        pathology_label = barcode_type[barcode_info[j][0]]
        component_label = record[5]
        X = barcode_info[j][1]
        Y = -barcode_info[j][2]
        opacity = np.float(record[8])   
        active_spot[j].append([pathology_label, component_label, X, Y, opacity])
        ''''''
        
    ######### color the spots in the plot with opacity = attention score #################
    opacity_list = []
    for i in active_spot:
        sum_opacity = []
        for edges in active_spot[i]:
            sum_opacity.append(edges[4])
            
        avg_opacity = np.max(sum_opacity) #np.mean(sum_opacity)
        opacity_list.append(avg_opacity) 
        active_spot[i]=[active_spot[i][0][0], active_spot[i][0][1], active_spot[i][0][2], active_spot[i][0][3], avg_opacity]

    min_opacity = np.min(opacity_list)
    max_opacity = np.max(opacity_list)

    #### making dictionary for converting to pandas dataframe to draw altair plot ###########
    data_list=dict()
    data_list['pathology_label']=[]
    data_list['component_label']=[]
    data_list['X']=[]
    data_list['Y']=[]   
    data_list['opacity']=[]  

    for i in range (0, len(barcode_info)):        
        if i in active_spot:
            data_list['pathology_label'].append(active_spot[i][0])
            data_list['component_label'].append(active_spot[i][1])
            data_list['X'].append(active_spot[i][2])
            data_list['Y'].append(active_spot[i][3])
            data_list['opacity'].append((active_spot[i][4]-min_opacity)/(max_opacity-min_opacity))
            
        else:
            data_list['pathology_label'].append(barcode_type[barcode_info[i][0]])
            data_list['component_label'].append(0) # make it zero so it is black
            data_list['X'].append(barcode_info[i][1])
            data_list['Y'].append(-barcode_info[i][2])
            data_list['opacity'].append(0.1)
            # barcode_info[i][3] = 0



    # converting to pandas dataframe

    data_list_pd = pd.DataFrame(data_list)
    id_label = len(list(set(data_list['component_label']))) # unique_component_count
    set1 = altairThemes.get_colour_scheme("Set1", id_label)
    set1[0] = '#000000'
    chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
        alt.X('X', scale=alt.Scale(zero=False)),
        alt.Y('Y', scale=alt.Scale(zero=False)),
        shape = alt.Shape('pathology_label:N'), #shape = "pathology_label",
        color=alt.Color('component_label:N', scale=alt.Scale(range=set1)),
        #opacity=alt.Opacity('opacity:N'), #"opacity", 
        tooltip=['component_label'] #,'opacity'
    )

    chart.save(output_name + args.data_name +'_Tcell_altair_plot.html')
    print('Altair plot generation done')

    ##################### save the top_edge_count in csv #############
    df = pd.DataFrame(csv_record_final[0:len(csv_record_final)-1])
    df.to_csv(output_name + args.data_name +'_top20p_Tcell.html', index=False, header=False)
    ###################################  Histogram plotting #################################################################################
    
    df = pd.DataFrame(csv_record_final)
    df.to_csv('temp_csv.csv', index=False, header=False)
    df = pd.read_csv('temp_csv.csv', sep=",")
    os.remove('temp_csv.csv') # delete the intermediate file

    print('len of loaded csv for histogram generation is %d'%len(df))
    df = preprocessDf(df)
    p = plot(df)
    outPath = output_name + args.data_name + '_Tcell_histogram_test.html'
    p.save(outPath)	
    print('Histogram plot generation done')

    if args.histogram_attention_score==1:
        lr_score = defaultdict(list)
        for i in range (1, len(csv_record_final)-1):    
            lr_score[csv_record_final[i][2]+'-'+csv_record_final[i][3]].append(csv_record_final[i][8])
        for key in lr_score.keys():
            lr_score[key]=np.sum(lr_score[key])

        # now plot the histograms where X axis will show the name or LR pair and Y axis will show the score.
        data_list=dict()
        data_list['X']=[]
        data_list['Y']=[] 
        for key in lr_score.keys(): #len(two_hop_pattern_distribution)):
            data_list['X'].append(key)
            data_list['Y'].append(lr_score[key])
            
        data_list_pd = pd.DataFrame({
            'Ligand-Receptor Pairs': data_list['X'],
            'Total Attention Score': data_list['Y']
        })
    
        chart = alt.Chart(data_list_pd).mark_bar().encode(
            x=alt.X("Ligand-Receptor Pairs:N", axis=alt.Axis(labelAngle=45), sort='-y'),
            y='Total Attention Score'
        )
    
        chart.save(output_name + args.data_name +'_LRpair_score.html')
        print('Saved at '+output_name + args.data_name +'_LRpair_score.html')
 '''
In [32]: lr_score
Out[32]: 
defaultdict(list,
            {'CCL21-CXCR4': 4225.405765723297,
             'CCL21-CCR7': 2949.097245159767,
             'PCDH7-CXCR4': 0.3704821288748415,
             'TGFB1-ENG': 4.080723634460016,
             'CXCL14-CXCR4': 0.5531352038615932,
             'TGFB1-TGFBR2': 26.625946890996737,
             'HLA-DRA': 0.37394397965082715,
             'CELSR1-HLA': 0.3759620404342901,
             'NPB-CCR7': 0.0929142342989421,
             'PTPRF-RACK1': 0.3762823146663877,
             'NPW-CCR7': 0.0929572614024158,
             'NPB-CXCR4': 0.371176383322931,
             'RET-RACK1': 0.37628637558792494,
             'NPW-CXCR4': 0.3713071203586996,
             'HMGB1-AR': 0.2741797884119278,
             'HMGB1-NTRK1': 0.27418142957807,
             'HMGB1-TLR9': 0.2741830609056528,
             'HMGB1-LY96': 0.2741846204268385,
      
 '''
    ################################ Density Curve #############################
    combined_score_distribution_ccl19_ccr7 = []
    combined_score_distribution = []
    # 63470 is the length of csv_record_final (number of records in Tcell zone)
    for k in range (1, len(csv_record_final)-1): 
        i = csv_record_final[k][6]
        j = csv_record_final[k][7]
        ligand = csv_record_final[k][2]
        receptor = csv_record_final[k][3]
        if ligand =='CCL19' and receptor == 'CCR7':
            combined_score_distribution_ccl19_ccr7.append(csv_record_final[k][8])
        else:
            combined_score_distribution.append(csv_record_final[k][8])
            
    some_dict = dict(A=combined_score_distribution, B=combined_score_distribution_ccl19_ccr7)
    
    df = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in some_dict.items()]))
    
    df = df.rename(columns={'A': 'all_pairs', 'B': 'CCL19_CCR7'})
    
    source = df
    #########################################################################
    chart = alt.Chart(source).transform_fold(
        ['all_pairs',
         'CCL19_CCR7'],
        as_ = ['distribution_type', 'value']
    ).transform_density(
        density = 'value',
        groupby=['distribution_type'],        
        steps=200
    ).mark_area(opacity=0.9).encode(
        alt.X('value:Q'),
        alt.Y('density:Q', stack=None ),
        alt.Color('distribution_type:N')
    )
    ####################### or ###################################### 
    '''
    chart = alt.Chart(source).transform_fold(
        ['all_pairs', 'CCL19_CCR7'],
        as_=['Distribution Type', 'Attention Score']
    ).mark_bar(
        opacity=0.5,
        binSpacing=0
    ).encode(
        alt.X('Attention Score:Q', bin=alt.Bin(maxbins=100)),
        alt.Y('count()', stack='zero'),
        alt.Color('Distribution Type:N')
    )
    '''
    chart.save(output_name + args.data_name +'_Tcell_attention_distribution.html')  
    ##################################################################################################################
    df = pd.read_csv(args.input_edge_file, sep=",", header=None)

    input_score_distribution_ccl19_ccr7 = []
    input_score_distribution = []
    # 63470 is the length of csv_record_final (number of records in Tcell zone)
    for k in range (0, len(df)): 
        i = df[0][k]
        j = df[1][k]
        ligand = df[3][k]
        receptor = df[4][k]
        score = df[2][k]
        if barcode_type[barcode_info[i][0]]=='T-cell' and barcode_type[barcode_info[j][0]]=='T-cell': 
            if ligand =='CCL19' and receptor == 'CCR7':
                input_score_distribution_ccl19_ccr7.append(score)
            else:
                input_score_distribution.append(score)
            
    some_dict = dict(A=input_score_distribution, B=input_score_distribution_ccl19_ccr7)
    
    df = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in some_dict.items()]))
    
    df = df.rename(columns={'A': 'all_pairs', 'B': 'CCL19_CCR7'})
    
    source = df
    #########################################################################
    chart = alt.Chart(source).transform_fold(
        ['all_pairs',
         'CCL19_CCR7'],
        as_ = ['distribution_type', 'value']
    ).transform_density(
        density = 'value',
        groupby=['distribution_type'],        
        steps=0.5
    ).mark_area(opacity=0.9).encode(
        alt.X('value:Q'),
        alt.Y('density:Q', stack=None ),
        alt.Color('distribution_type:N')
    )
    ####################### or ###################################### 
    '''
    chart = alt.Chart(source).transform_fold(
        ['all_pairs', 'CCL19_CCR7'],
        as_=['Distribution Type', 'Attention Score']
    ).mark_bar(
        opacity=0.5,
        binSpacing=0
    ).encode(
        alt.X('Attention Score:Q', bin=alt.Bin(maxbins=100)),
        alt.Y('count()', stack='zero'),
        alt.Color('Distribution Type:N')
    )
    '''
    ###########################################################################
    chart.save(output_name + args.data_name +'_Tcell_input_coexpression_distribution.html')  
    

