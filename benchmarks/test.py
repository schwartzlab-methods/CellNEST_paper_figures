import os
import pandas as pd
import copy
import csv
import numpy as np
import sys
from collections import defaultdict
import pickle
import gzip
import matplotlib.pyplot as plt
import altairThemes
import altair as alt
import yaml
import pandas as pd
import glob
from collections import defaultdict
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

# filter by strength (default = positive ccc)
def filter_ccc_df(
        df, 
        strength_col, 
        strength_thresh = 0
):
    df_filtered = df[df[strength_col] > strength_thresh]
    return df_filtered 

# convert ground-truth dictionary from spot to cluster resolution
def convert_gt_dict(
    gt_cell_dict,
    cluster_df: pd.DataFrame
):
    gt_cell_dict = dict(gt_cell_dict)
    spot_to_cluster = cluster_df.set_index("spot")["cluster"].to_dict() 
    gt_cluster_dict = {}
    for sender_cell, ccc_dict in gt_cell_dict.items():
        sender_cluster = spot_to_cluster[sender_cell]
        print(f"sender Cell: {sender_cell}, cluster: {sender_cluster}")
        if sender_cluster not in gt_cluster_dict:
            gt_cluster_dict[sender_cluster] = {}
        for receiver_cell, lr_list in ccc_dict.items():
            receiver_cluster = spot_to_cluster[receiver_cell]
            print(f"receiver Cell: {receiver_cell}, cluster: {receiver_cluster}")
            if receiver_cluster not in gt_cluster_dict[sender_cluster]:
                gt_cluster_dict[sender_cluster][receiver_cluster] = []
            lr_tuple = tuple(lr_list)
            print(f"LR pair: {lr_tuple}")
            if lr_tuple not in gt_cluster_dict[sender_cluster][receiver_cluster]:
                print(f"adding LR pair {lr_tuple} to sender cluster {sender_cluster}, receiver cluster {receiver_cluster}")
                gt_cluster_dict[sender_cluster][receiver_cluster].append(lr_tuple)
            print(f"current dictionary: {gt_cluster_dict}")
    print(gt_cluster_dict)
    return gt_cluster_dict


def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    datasets = config["params"]["datasets"]
    tools = config["params"]["tools"]
    ##### until created ######
    tools.remove("TWCOM")
    datasets.remove("lymph_node")
    data_dir = config["directories"]["data"]
    raw_result_dir = config["directories"]["raw_result"]

    for dataset in datasets:
        coordinate_file = glob.glob(os.path.join(data_dir, dataset, "*coordinate"))[0]
        ground_truth_file = glob.glob(os.path.join(data_dir, dataset, "*ground_truth_ccc"))[0]
        with gzip.open(coordinate_file, 'rb') as fp:
            x_index, y_index, no_need = pickle.load(fp)
        with gzip.open(ground_truth_file, 'rb') as fp:  
            lr_database, lig_rec_dict_TP_cell, random_activation = pickle.load(fp)
        print(lig_rec_dict_TP_cell)
        cluster_df = pd.read_csv(os.path.join(data_dir, dataset, "leiden.csv"))
        # convert ground truth ligand-receptor dictionary from spot to cluster resolution 
        lig_rec_dict_TP = convert_gt_dict(lig_rec_dict_TP_cell, cluster_df)
        # print(len(lig_rec_dict_TP))
        # datapoint_size = x_index.shape[0] ### total number of cells or datapoints 
        datapoint_size = len(cluster_df["cluster"].unique())
        print(datapoint_size)
        tp = 0 # true positives
        for i in lig_rec_dict_TP.keys():
            for j in lig_rec_dict_TP[i].keys():
                tp = tp + len(lig_rec_dict_TP[i][j]) # list in fatema's case
        positive_class = tp # WE NEED THIS TO CALCULATE 'TRUE POSITIVE RATE'
        print(positive_class)
        for tool in tools:
            filtered_result_dir = os.path.join(config["directories"]["filtered_result"], tool)
            os.makedirs(filtered_result_dir, exist_ok = True)
            ccc_result = pd.read_csv(os.path.join(raw_result_dir, tool, f"{dataset}.csv"))
            col_dict = tool_params(tool)
            ccc_result = filter_ccc_df(ccc_result, col_dict["strength"])
            lig_rec_dict = {}
            attention_scores = {}
            distribution = []
            for _, row in ccc_result.iterrows():
                sender = row[col_dict["sender"]]
                receiver = row[col_dict["receiver"]]
                ligand = row[col_dict["ligand"]]
                receptor = row[col_dict["receptor"]]
                strength = row[col_dict["strength"]]
                if sender not in lig_rec_dict:
                    lig_rec_dict[sender] = {}
                if receiver not in lig_rec_dict[sender]:
                    lig_rec_dict[sender][receiver] = []
                lr_pair = (ligand, receptor)
                if lr_pair not in lig_rec_dict[sender][receiver]:
                    lig_rec_dict[sender][receiver].append(lr_pair)
                if sender not in attention_scores:
                    attention_scores[sender] = {}
                if receiver not in attention_scores[sender]:
                    attention_scores[sender][receiver] = []
                attention_scores[sender][receiver].append(strength)
                distribution.append(strength)
            print(len(lig_rec_dict))
            print(len(attention_scores))
            print(len(distribution)) 
            confusion_matrix = np.zeros((2,2))
            print(lig_rec_dict.keys())
            print(lig_rec_dict)
            for i in range (0, datapoint_size):
                for j in range (0, datapoint_size):
                    # lr_pair_list = lig_rec_dict[i][j]
                    lr_pair_list = lig_rec_dict.get(i, {}).get(j, [])
                    ground_truth = lig_rec_dict_TP.get(i, {}).get(j, [])
                    # if len(lr_pair_list)>0:
                    for k in lr_pair_list: # tuple
                        # print(f"Tool detected: {k}, Ground truth: {lig_rec_dict_TP[i][j]}") 
                        # if i in lig_rec_dict_TP and j in lig_rec_dict_TP[i] and k in lig_rec_dict_TP[i][j]:
                        if k in ground_truth:
                            confusion_matrix[0][0] += 1
                            #print("i=%d j=%d k=%d"%(i, j, k))
                            # confusion_matrix[0][0] = confusion_matrix[0][0] + 1 # detected true positive by COMMOT
                        else:
                            confusion_matrix[1][0] += 1
                            # confusion_matrix[1][0] = confusion_matrix[1][0] + 1 #  detected false positive by COMMOT            
            # print(confusion_matrix[0][0]) 
            negative_class = len(distribution) - confusion_matrix[0][0] # WE NEED THIS TO CALCULATE 'FALSE POSITIVE RATE'
            print(negative_class)
            distribution = sorted(distribution, reverse=True) 
            # start roc plot here. select top 10% (90th), 20% (80th), 30% (70th), ... ccc and calculate TPR and FPR 
            plot_dict = defaultdict(list)
            for percentile_value in [90, 80, 70, 60, 50, 40, 30, 20, 10, 0]:
                threshold_percentile =  np.percentile(distribution, percentile_value)
                existing_lig_rec_dict = [] # record COMMOT detected edges that are above the threshold percentile attention score
                for i in range (0, datapoint_size):
                    existing_lig_rec_dict.append([])   
                    for j in range (0, datapoint_size):	
                        existing_lig_rec_dict[i].append([])   
                        existing_lig_rec_dict[i][j] = []

                total_edges_count = 0
                for i in range (0, datapoint_size):
                    for j in range (0, datapoint_size):
                        if i in attention_scores and j in attention_scores[i]: # DP added check 
                            atn_score_list = attention_scores[i][j]
                            for k in range (0, len(atn_score_list)):
                                if attention_scores[i][j][k] >= threshold_percentile: 
                                # connecting_edges[i][j] = 1
                                    existing_lig_rec_dict[i][j].append(lig_rec_dict[i][j][k])
                                    total_edges_count = total_edges_count + 1
                ############# 
                #print('total edges %d'%total_edges_count)
                #negative_class = 0
                confusion_matrix = np.zeros((2,2))
                for i in range (0, datapoint_size):
                    for j in range (0, datapoint_size):
                        # if len(existing_lig_rec_dict[i][j])>0:
                        for k in existing_lig_rec_dict[i][j]:   
                            if i in lig_rec_dict_TP and j in lig_rec_dict_TP[i] and k in lig_rec_dict_TP[i][j]:
                                    #print("i=%d j=%d k=%d"%(i, j, k))
                                confusion_matrix[0][0] = confusion_matrix[0][0] + 1
                            else:
                                confusion_matrix[1][0] = confusion_matrix[1][0] + 1                 
                print(confusion_matrix[0][0]) 
                print('%d, %g, %g'%(percentile_value,  (confusion_matrix[1][0]/negative_class)*100, (confusion_matrix[0][0]/positive_class)*100))    
                FPR_value = (confusion_matrix[1][0]/negative_class)#*100
                TPR_value = (confusion_matrix[0][0]/positive_class)#*100
                plot_dict['FPR'].append(FPR_value)
                plot_dict['TPR'].append(TPR_value)
                plot_dict['Type'].append('COMMOT') # no noise

            # with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/synthetic_data/type_gaussian_distribution/no_noise/" + options +'_COMMOT_roc', 'wb') as fp: #b, b_1, a  11to20runs
                # pickle.dump(plot_dict, fp) #a - [0:5]

            data_list_pd = pd.DataFrame(plot_dict)    
            chart = alt.Chart(data_list_pd).mark_line().encode(
                x='FPR:Q',
                y='TPR:Q',
                color='Type:N',
            )	
            chart.save(os.path.join(filtered_result_dir, f"{dataset}_roc.html"))

        #with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/synthetic_data/type_gaussian_distribution/no_noise/uniform_distribution_input_graph" , 'rb') as fp:           
        #    row_col, edge_weight, lig_rec  = pickle.load(fp) 

        ######################### COMMOT ###############################################################################################################
        # options = 'uniform_distribution_no_noise'
        # with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/synthetic_data/type_gaussian_distribution/no_noise/" + options + '_commot_result', 'rb') as fp:
        #     attention_scores, lig_rec_dict, distribution = pickle.load(fp)            

        # lig_rec_dict[i][j]=[...] # is a list of lr pairs (edges) between cell i and cell j 
        # attention_scores[i][j]=[...] # is a list of attention scores of the lr pairs (edges) between cell i and cell j
        # distribution=[...] is a combined list of attention scores of all edges. 

        # let's first find the total count of negative classes in the COMMOT result. 
        # We calculate the tp detected by COMMOT and then deduct it from total detection by COMMOT to get the negative classes. 

main()
