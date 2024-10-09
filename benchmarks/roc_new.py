import os
import pandas as pd
import copy
import csv
import numpy as np
import sys
from collections import defaultdict
import pickle
import gzip
import yaml
import matplotlib.pyplot as plt
import altairThemes
import altair as alt
import glob
alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")

def tool_params(ccc_tool):
    tool_dict = {
        "cellchat": {
            "sender": "source",
            "receiver": "target",
            "ligand": "ligand",
            "receptor": "receptor",
            "lri": "interaction_name",
            "strength": "prob"
        },
        "giotto": {
            "sender": "lig_cell_type",
            "receiver": "rec_cell_type",
            "ligand": "ligand",
            "receptor": "receptor",
            "lri": "LR_comb",
            "strength": "PI"  # PI: significance score: log2fc * -log10(p.adj)
        },
        "TWCOM": {
            "sender": "CellType_Sender",
            "receiver": "CellType_Receiver",
            "lri": "lri",
            "strength": "Estimate"
        }
    }

    return tool_dict[ccc_tool]

# convert ground truth dictionary from cell to cluster resolution 
def convert_ground_truth(
    cell_ground_truth_dict,
    cluster_df: pd.DataFrame
):
    cell_ground_truth_dict = dict(cell_ground_truth_dict)
    mapping_dict = cluster_df.set_index("spot")["cluster"].to_dict() 
    cluster_ground_truth_dict = {}
    for sender_cell, ccc_dict in cell_ground_truth_dict.items():
        sender_cluster = mapping_dict[sender_cell]
        if sender_cluster not in cluster_ground_truth_dict:
            cluster_ground_truth_dict[sender_cluster] = {}
        for receiver_cell, lr_list in ccc_dict.items():
            receiver_cluster = mapping_dict[receiver_cell]
            if receiver_cluster not in cluster_ground_truth_dict[sender_cluster]:
                cluster_ground_truth_dict[sender_cluster][receiver_cluster] = []
            for k in lr_list:
                if k not in cluster_ground_truth_dict[sender_cluster][receiver_cluster]:
                    cluster_ground_truth_dict[sender_cluster][receiver_cluster].append(k)

    return cluster_ground_truth_dict

# filter ligand-receptor interactions by strength (default = + strength)
def filter_ccc_df(
        df, 
        strength_col, 
        strength_thresh = 0
):
    df_filtered = df[df[strength_col] > strength_thresh]

    return df_filtered 

def load_tool_results(
        datapoint_size,
        tool_result, 
        tool,
        lr_db
):

    col_dict = tool_params(tool)
    sender_col = col_dict['sender']
    receiver_col = col_dict['receiver']
    ligand_col = col_dict['ligand']
    receptor_col = col_dict['receptor']
    strength_col = col_dict['strength']

    tool_result = filter_ccc_df(tool_result, strength_col)    

    attention_scores = []
    lig_rec_dict = []
    distribution = []
    for i in range (0, datapoint_size):
        attention_scores.append([])   
        lig_rec_dict.append([])   
        for j in range (0, datapoint_size):	
            attention_scores[i].append([])   
            attention_scores[i][j] = []
            lig_rec_dict[i].append([])   
            lig_rec_dict[i][j] = []
    LR_pairs = lr_db.apply(lambda row: f"{row['ligand']}-{row['receptor']}", axis=1).tolist()
    lri_map = {pair: idx for idx, pair in enumerate(LR_pairs)}
    for i in range(datapoint_size):
        for j in range(datapoint_size):
            result = tool_result[(tool_result[sender_col] == i) & (tool_result[receiver_col] == j)]
            if not result.empty:
                for _, row in result.iterrows():
                    lri = f"{row[ligand_col]}-{row[receptor_col]}"
                    if lri in lri_map:
                        lig_rec_dict[i][j].append(lri_map[lri])
                        attention_scores[i][j].append(row[strength_col])
                        distribution.append(row[strength_col])

    return lig_rec_dict, attention_scores, distribution

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    datasets = config["params"]["datasets"]
    tools = config["params"]["tools"]
    tools.remove("TWCOM") ##### until created ######
    datasets.remove("lymph_node")
    data_dir = config["directories"]["data"]
    raw_result_dir = config["directories"]["raw_result"]
    filtered_result_dir = os.path.join(config["directories"]["filtered_result"])

    for dataset in datasets:
        coordinate_file = glob.glob(os.path.join(data_dir, dataset, "*coordinate"))[0]
        ground_truth_file = glob.glob(os.path.join(data_dir, dataset, "*ground_truth_ccc"))[0]
        with gzip.open(coordinate_file, 'rb') as fp:
            x_index, y_index, no_need = pickle.load(fp)
        with gzip.open(ground_truth_file, 'rb') as fp:  
            _, lig_rec_dict_TP_cell, random_activation = pickle.load(fp)
        lr_database = pd.read_csv(glob.glob(os.path.join(data_dir, dataset, "*lrdb.csv"))[0])
        cluster_df = pd.read_csv(os.path.join(data_dir, dataset, "leiden.csv"))
        # convert ground truth dictionary from cell to cluster resolution
        lig_rec_dict_TP = convert_ground_truth(lig_rec_dict_TP_cell, cluster_df)
        datapoint_size = len(cluster_df["cluster"].unique())

        tp = 0 # true positives 
        for i in lig_rec_dict_TP.keys():
            for j in lig_rec_dict_TP[i].keys():
                tp = tp + len(lig_rec_dict_TP[i][j])
        positive_class = tp

        results_list = []
        for tool in tools: 
            # load tool's results (filtered for p-value already)
            tool_result = pd.read_csv(os.path.join(raw_result_dir, tool, f"{dataset}.csv"))
            lig_rec_dict, attention_scores, distribution = load_tool_results(datapoint_size, tool_result, tool, lr_database)
            distribution = sorted(distribution, reverse = True)
            min_limit =  distribution[len(distribution) - 1] # min score to be considered

            detected_TP = 0
            for i in range (0, datapoint_size):
                for j in range (0, datapoint_size):
                    lr_pair_list = lig_rec_dict[i][j] # list of values 
                    if len(lr_pair_list) > 0:
                        for k in lr_pair_list:  
                            # if attention_scores[i][j][k] < min_limit:
                                # continue # ignore                
                            if i in lig_rec_dict_TP and j in lig_rec_dict_TP[i] and k in lig_rec_dict_TP[i][j]:
                                #print("i=%d j=%d k=%d"%(i, j, k))
                                detected_TP = detected_TP + 1 # detected true positive
                
            negative_class = len(distribution) - detected_TP # WE NEED THIS TO CALCULATE 'FALSE POSITIVE RATE'
            print(negative_class)
            print(detected_TP)

            plot_dict = defaultdict(list)
            for percentile_value in [90, 80, 70, 60, 50, 40, 30, 20, 10, 0]:
                threshold_percentile =  np.percentile(distribution, percentile_value)
                existing_lig_rec_dict = [] # record COMMOT detected edges that are above the threshold percentile attention score
                for i in range (0, datapoint_size):
                    existing_lig_rec_dict.append([])   
                    for j in range (0, datapoint_size):	
                        existing_lig_rec_dict[i].append([])   
                        existing_lig_rec_dict[i][j] = []

                # connecting_edges = np.zeros((datapoint_size, datapoint_size))
                total_edges_count = 0
                for i in range (0, datapoint_size):
                    for j in range (0, datapoint_size):
                        atn_score_list = attention_scores[i][j]
                        for k in range (0, len(atn_score_list)):
                            if attention_scores[i][j][k] >= threshold_percentile: # and attention_scores[i][j][k] <= max_limit: 
                                # connecting_edges[i][j] = 1
                                existing_lig_rec_dict[i][j].append(lig_rec_dict[i][j][k])
                                total_edges_count = total_edges_count + 1
                                
                confusion_matrix = np.zeros((2,2))
                for i in range (0, datapoint_size):
                    for j in range (0, datapoint_size):
                        if len(existing_lig_rec_dict[i][j])>0:
                            for k in existing_lig_rec_dict[i][j]:   
                                if i in lig_rec_dict_TP and j in lig_rec_dict_TP[i] and k in lig_rec_dict_TP[i][j]:
                                    #print("i=%d j=%d k=%d"%(i, j, k))
                                    confusion_matrix[0][0] = confusion_matrix[0][0] + 1
                                else:
                                    confusion_matrix[1][0] = confusion_matrix[1][0] + 1                 
                        
                print('%d, %g, %g'%(percentile_value,  (confusion_matrix[1][0]/negative_class)*100, (confusion_matrix[0][0]/positive_class)*100))    
                FPR_value = (confusion_matrix[1][0]/negative_class)#*100
                TPR_value = (confusion_matrix[0][0]/positive_class)#*100
                plot_dict['FPR'].append(FPR_value)
                plot_dict['TPR'].append(TPR_value)
                plot_dict['Type'].append(tool) 
                results_list.append({'FPR': FPR_value, 'TPR': TPR_value, 'Tool': tool})

            # data_list_pd = pd.DataFrame(plot_dict) 
            results_df = pd.DataFrame(results_list)

            chart = alt.Chart(results_df).mark_line().encode(
                x='FPR:Q',
                y='TPR:Q',
                color='Tool:N',
            )	
            chart.save(os.path.join(filtered_result_dir, f"{dataset}.html"))

main()









