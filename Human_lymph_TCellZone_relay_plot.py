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
  df["component"] = df["component"]

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
    parser.add_argument( '--top_edge_count', type=int, default=400000, help='Number of the top communications to plot. To plot all insert -1') # 
    parser.add_argument( '--barcode_info_file', type=str, default='NEST_figures_input/V1_Human_Lymph_Node_spatial_barcode_info', help='Path to load the barcode information file produced during data preprocessing step')
    parser.add_argument( '--annotation_file_path', type=str, default='NEST_figures_input/V1_Human_Lymph_Node_spatial_annotation.csv', help='Path to load the annotation file in csv format (if available) ') #_ayah_histology
    parser.add_argument( '--selfloop_info_file', type=str, default='NEST_figures_input/V1_Human_Lymph_Node_spatial_self_loop_record', help='Path to load the selfloop information file produced during data preprocessing step')
    parser.add_argument( '--top_ccc_file', type=str, default='NEST_figures_input/V1_Human_Lymph_Node_spatial_top20percent.csv', help='Path to load the selected top CCC file produced during data postprocessing step')
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
    if args.annotation_file_path != '':
        annotation_data = pd.read_csv(args.annotation_file_path, sep=",")
        pathologist_label=[]
        for i in range (0, len(annotation_data)):
            pathologist_label.append([annotation_data['Barcode'][i], annotation_data['Type'][i]])
    
        barcode_type=dict() # record the type (annotation) of each spot (barcode)
        for i in range (0, len(pathologist_label)):
            barcode_type[pathologist_label[i][0]] = pathologist_label[i][1]

    else:
        barcode_type=dict() # record the type (annotation) of each spot (barcode)
        for i in range (0, len(barcode_info)):
            barcode_type[barcode_info[i][0]] = ''
      
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

    
	csv_record_final = [df_column_names] + csv_record #[0:min(args.top_edge_count, len(csv_record))]

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
#    print('Number of connected components %d'%n_components) 

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
    
    #####################################
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

    ####################### keep top 5000 ###########################################################################
'''
    ## sort the edges based on their rank (column 4), low to high, low being higher attention score
    csv_record_temp = csv_record_final[1:len(csv_record_final)-1]
    csv_record_temp = sorted(csv_record_temp, key = lambda x: x[4])

    csv_record_temp = csv_record_temp[0:5000] + [csv_record_final[len(csv_record_final)-1]]
    csv_record_temp = [csv_record_final[0]] + csv_record_temp  
    csv_record_final = csv_record_temp 
    '''
####################### pattern finding ##########################################################################
    # make a dictionary to keep record of all the outgoing edges [to_node, ligand, receptor] for each node
 
    each_node_outgoing = defaultdict(list)
    for k in range (1, len(csv_record_final)-1): # last record is a dummy for histogram preparation
        i = csv_record_final[k][6]
        j = csv_record_final[k][7]
        if i == j:
            continue        
        ligand = csv_record_final[k][2]
        receptor = csv_record_final[k][3]
        each_node_outgoing[i].append([j, ligand, receptor, k]) 
    
    # all possible 2-hop pattern count
    pattern_distribution = defaultdict(list)
    pattern_distribution_cell_info = defaultdict(list)
    # pattern_distribution['ligand-receptor to ligand-receptor']=[1,1,1,1, ...]
    edge_list_2hop = []
    target_relay = 'CCL19-CCR7 to CCL21-CXCR4'
    for i in each_node_outgoing:
        for tupple in each_node_outgoing[i]: # first hop
            j = tupple[0]
            lig_rec_1 = tupple[1]+'-'+tupple[2]
            record_id_1 = tupple[3]
            if j in each_node_outgoing:
                for tupple_next in each_node_outgoing[j]: # second hop
                    k = tupple_next[0]
                    if k == i or k == j:
                        continue
                    lig_rec_2 = tupple_next[1]+'-'+tupple_next[2]
                    record_id_2 = tupple_next[3]
                    pattern_distribution[lig_rec_1 + ' to ' + lig_rec_2].append(1)
                    pattern_distribution_cell_info[lig_rec_1 + ' to ' + lig_rec_2].append([barcode_info[i], barcode_info[j], barcode_info[k]]) # a_id, b_id, c_id 
                    relay = lig_rec_1 + ' to ' + lig_rec_2
                    if relay == target_relay:
                        edge_list_2hop.append([record_id_1,record_id_2])

	
    two_hop_pattern_distribution = []
    same_count = 0
    for key in pattern_distribution:
        count = len(pattern_distribution[key])
        two_hop_pattern_distribution.append([key, count]) 
        #if lig_rec_1 == lig_rec_2:
        #    same_count = same_count + 1
    
    two_hop_pattern_distribution = sorted(two_hop_pattern_distribution, key = lambda x: x[1], reverse=True) # high to low
    
    data_list=dict()
    data_list['X']=[]
    data_list['Y']=[] 
    for i in range (0, 100): #len(two_hop_pattern_distribution)):
        data_list['X'].append(two_hop_pattern_distribution[i][0])
        data_list['Y'].append(two_hop_pattern_distribution[i][1])
        
    data_list_pd = pd.DataFrame({
        'Relay Patterns': data_list['X'],
        'Pattern Abundance (#)': data_list['Y']
    })
  
    chart = alt.Chart(data_list_pd).mark_bar().encode(
        x=alt.X("Relay Patterns:N", axis=alt.Axis(labelAngle=45), sort='-y'),
        y='Pattern Abundance (#)'
    )

    chart.save(output_name + args.data_name + '_TcellZone' +'_pattern_distribution.html')
    ######################### save as table format and numpy format ########
    data_list_pd.to_csv(output_name + args.data_name +'_top20p_Tcell_topEdge'+str(args.top_edge_count)+'_relay_count.csv', index=False)
  
    with gzip.open(output_name + args.data_name + 'top20p_Tcell_topEdge'+str(args.top_edge_count)+'_pattern_distribution_cell_info', 'wb') as fp: 
	    pickle.dump(pattern_distribution_cell_info, fp)

    with gzip.open(output_name + args.data_name + 'top20p_Tcell_topEdge'+str(args.top_edge_count)+'_pattern_distribution_cell_info', 'rb') as fp:  
	    pattern_distribution_cell_info = pickle.load(fp)

    pattern_list = list(pattern_distribution_cell_info.keys()) 
    # see which cells are forming the first pattern
    group_list = pattern_distribution_cell_info[pattern_list[0]]
    for group_id in range(0, len(group_list)):
        print('group: %d'%group_id)
        print(group_list[group_id])

    ######################### Plotting the two hops #####################
    set1 = altairThemes.get_colour_scheme("Set1", unique_component_count)
    colors = set1
    colors[0] = '#000000' # black means no CCC
    ids = []
    x_index=[]
    y_index=[]
    colors_point = []
    for i in range (0, len(barcode_info)):    
        ids.append(i)
        x_index.append(barcode_info[i][1])
        y_index.append(barcode_info[i][2])    
        colors_point.append(colors[barcode_info[i][3]]) 
    
    max_x = np.max(x_index)
    max_y = np.max(y_index)



    g = nx.MultiDiGraph(directed=True) 
    for i in range (0, len(barcode_info)):        
        marker_size = 'circle'
        label_str =  str(i)+'_c:'+str(barcode_info[i][3]) #  label of the node or spot is consists of: spot id, component number
        if args.annotation_file_path != '': 
            if 'T-cell' in str(barcode_type[barcode_info[i][0]]): # Tumour
                marker_size = 'box'
                label_str = label_str +'_'+ 'T-cell' # also add the type of the spot to the label if annotation is available 
        
        g.add_node(int(ids[i]), x=int(x_index[i]), y=int(y_index[i]), label = label_str, pos = str(x_index[i])+","+str(-y_index[i])+" !", physics=False, shape = marker_size, color=matplotlib.colors.rgb2hex(colors_point[i]))    


    # scale the edge scores [0 to 1] to make width in the plot work
    score_list = []
    for k in range (1, len(csv_record_final)-1):
        score_list.append(csv_record_final[k][8])

    min_score = np.min(score_list)
    max_score = np.max(score_list)

    count_edges = 0
    for relay in edge_list_2hop:
        print(relay)
        for hop in range (0, len(relay)):
            k = relay[hop] # record id for the hops
            i = csv_record_final[k][6]
            j = csv_record_final[k][7] 
            print('%d %d'%(i, j))
            ligand = csv_record_final[k][2]
            receptor = csv_record_final[k][3]            
            edge_score = csv_record_final[k][8]
            edge_score = (edge_score-min_score)/(max_score-min_score)   
            title_str =  "L:" + ligand + ", R:" + receptor+ ", "+ str(edge_score) #+
            g.add_edge(int(i), int(j), label = title_str, color='#FF0000', value=np.float64(edge_score)) #colors_point[i]
            count_edges = count_edges + 1

    print("total edges plotted: %d"%count_edges)

    nt = Network( directed=True, height='1000px', width='100%') #"500px", "500px",, filter_menu=True     
    nt.from_nx(g)
    nt.save_graph(output_name + args.data_name + '_TcellZone' +'_relay_mygraph.html')
    print('Edge graph plot generation done')
    ########################################################################
    # convert it to dot file to be able to convert it to pdf or svg format for inserting into the paper
    write_dot(g, output_name + args.data_name +  '_TcellZone' + "_relay_test_interactive.dot")
    print('dot file generation done')
