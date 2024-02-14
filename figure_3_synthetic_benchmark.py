import os
import pandas as pd
import copy
import csv
import numpy as np
import sys
import altair as alt
from collections import defaultdict
import scipy
import pickle
import gzip
import matplotlib.pyplot as plt


 


import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--data_path', type=str, default='NEST_figures_input_synthetic/' , help='The path to dataset') 
parser.add_argument( '--output_path', type=str, default='NEST_figures_output/' , help='The path to dataset') 

args = parser.parse_args()
save_path = args.output_path
############################## Equidistant #################################################


sample_name = ["equidistant_withCCCpattern_lr1467_cellCount3000_noise0", 
              "equidistant_withCCCpattern_lowNoise_lr1467_cellCount3000_noise30",
              "equidistant_withCCCpattern_highNoise_lr1467_cellCount3000_noise30"]

sample_type = ["", "_LowNoise", "_HighNoise"]

sample_name_alt = ["equidistant_withCCCpattern_knn_lr1467_cellCount3000_noise0", 
              "equidistant_withCCCpattern_knn_lowNoise_lr1467_cellCount3000_noise30",
              "equidistant_withCCCpattern_knn_highNoise_lr1467_cellCount3000_noise30"]

 
output_name = ["c", "d", "e"]

for t in range (0, 3):
    plot_dict = defaultdict(list)
    with gzip.open(args.data_path + sample_name[t] +'_'+'naive_model', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("Naive"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("Naive"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    with gzip.open(args.data_path + sample_name[t]  +'_'+'rank_product_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]

    plot_dict_temp['FPR'].append(1.0)
    plot_dict_temp['TPR'].append(1.0)
    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST"+sample_type[t]) #(plot_dict_temp['Type'][i])
    

    with gzip.open(args.data_path + sample_name_alt[t]  +'_'+'rank_product_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]

    plot_dict_temp['FPR'].append(1.0)
    plot_dict_temp['TPR'].append(1.0)
    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST_alternate_cutOff"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST_alternate_cutOff"+sample_type[t]) #(plot_dict_temp['Type'][i])
 

######################
    with gzip.open(args.data_path + sample_name[t] +'_'+'rank_product_relu_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
#    plot_dict_temp['FPR'].append(1.0)
#    plot_dict_temp['TPR'].append(1.0)
#    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    

    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST_ReLU"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST_ReLU"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    with gzip.open(args.data_path + sample_name[t]  +'_'+'COMMOT', 'rb') as fp: # t = 0,1
        plot_dict_temp = pickle.load(fp) #a - [0:5]
        
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append('COMMOT'+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append('COMMOT'+sample_type[t]) #(plot_dict_temp['Type'][i])
    
        


###################

	
    data_list_pd = pd.DataFrame(plot_dict)    
    chart = alt.Chart(data_list_pd).mark_line().encode(
        x='FPR:Q',
        y='TPR:Q',
        color='Type:N',
    )	

    chart.save(save_path+'plot_equidistant_figure_'+output_name[t]+'.html')


#############################################  Uniform   ##########################################################
sample_type = ["", "_LowNoise", "_HighNoise"]

sample_name = ["uniform_distribution_withCCCpattern_lr112_cellCount5000_noise0", 
              "uniform_distribution_withCCCpattern_lowNoise_lr112_cellCount5000_noise30",
              "uniform_distribution_withCCCpattern_highNoise_lr112_cellCount5000_noise30"]
sample_name_alt = ["uniform_distribution_withCCCpattern_knn_lr112_cellCount5000_noise0", 
              "uniform_distribution_withCCCpattern_knn_lowNoise_lr112_cellCount5000_noise30",
              "uniform_distribution_withCCCpattern_knn_highNoise_lr112_cellCount5000_noise30"]

output_name = ["h", "i", "j"]

for t in range (0, 3): #len(sample_name)):
    plot_dict = defaultdict(list)
    with gzip.open(args.data_path +  sample_name[t] +'_'+'naive_model', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("Naive"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("Naive"+sample_type[t]) #(plot_dict_temp['Type'][i])
    ###
    
    with gzip.open(args.data_path + sample_name[t] +'_'+'rank_product_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict_temp['FPR'].append(1.0)
    plot_dict_temp['TPR'].append(1.0)
    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######
    with gzip.open(args.data_path + sample_name_alt[t] +'_'+'rank_product_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict_temp['FPR'].append(1.0)
    plot_dict_temp['TPR'].append(1.0)
    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST_alternate_cutOff"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST_alternate_cutOff"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######
    with gzip.open(args.data_path + sample_name[t] +'_'+'rank_product_relu_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
#    plot_dict_temp['FPR'].append(1.0)
#    plot_dict_temp['TPR'].append(1.0)
#    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST_ReLU"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST_ReLU"+sample_type[t]) #(plot_dict_temp['Type'][i])

 
    
    ######
    with gzip.open(args.data_path + sample_name[t]  +'_'+'Niches', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
        
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append('Niches'+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append('Niches'+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######
    if t!=2:
        with gzip.open(args.data_path + sample_name[t]  +'_'+'COMMOT', 'rb') as fp:
            plot_dict_temp = pickle.load(fp) #a - [0:5]
        # t = 2 -- did not work
            
        plot_dict['FPR'].append(0)
        plot_dict['TPR'].append(0)
        plot_dict['Type'].append('COMMOT'+sample_type[t]) #(plot_dict_temp['Type'][0])
        for i in range (0, len(plot_dict_temp['Type'])):
            plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
            plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
            plot_dict['Type'].append('COMMOT'+sample_type[t]) #(plot_dict_temp['Type'][i])
    
        

    
    data_list_pd = pd.DataFrame(plot_dict)    
    chart = alt.Chart(data_list_pd).mark_line().encode(
        x='FPR:Q',
        y='TPR:Q',
        color='Type:N',
    )	

    chart.save(save_path+'plot_uniform_figure_'+output_name[t]+'.html')

################################################################################# mixture of distribution #########################
sample_type = ["", "_LowNoise", "_HighNoise"]
sample_name = ["mixture_of_distribution_withCCCpattern_lr112_cellCount5000_noise0", 
              "mixture_of_distribution_withCCCpattern_lowNoise_lr112_cellCount5000_noise30",
              "mixture_of_distribution_withCCCpattern_highNoise_lr112_cellCount5000_noise30"]

sample_name_alt = ["mixture_of_distribution_withCCCpattern_thresholdDistance_lr112_cellCount5000_noise0", 
              "mixture_of_distribution_withCCCpattern_thresholdDistance_lowNoise_lr112_cellCount5000_noise30",
              "mixture_of_distribution_withCCCpattern_thresholdDistance_highNoise_lr112_cellCount5000_noise30"]

output_name = ["m", "n", "o"]


for t in range (0, 3): #len(sample_name)):
    plot_dict = defaultdict(list)
    with gzip.open(args.data_path +  sample_name[t] +'_'+'naive_model', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("Naive"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("Naive"+sample_type[t]) #(plot_dict_temp['Type'][i])
    ###

    with gzip.open(args.data_path +  sample_name[t] +'_'+'rank_product_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict_temp['FPR'].append(1.0)
    plot_dict_temp['TPR'].append(1.0)
    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######
    with gzip.open(args.data_path + sample_name_alt[t] +'_'+'rank_product_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict_temp['FPR'].append(1.0)
    plot_dict_temp['TPR'].append(1.0)
    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST_alternate_cutOff"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST_alternate_cutOff"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######
    with gzip.open(args.data_path +  sample_name[t] +'_'+'rank_product_relu_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    #plot_dict_temp['FPR'].append(1.0)
    #plot_dict_temp['TPR'].append(1.0)
    #plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST_relu"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST_relu"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######
    
    with gzip.open(args.data_path +  sample_name[t]  +'_'+'Niches', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
        
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append('Niches'+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append('Niches'+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######
    with gzip.open(args.data_path +  sample_name[t]  +'_'+'COMMOT', 'rb') as fp: #
    #with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" + sample_name[t]  +'_'+'COMMOT', 'rb') as fp: #
        plot_dict_temp = pickle.load(fp) #a - [0:5]
        
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append('COMMOT'+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append('COMMOT'+sample_type[t]) #(plot_dict_temp['Type'][i])
    
       
    
    data_list_pd = pd.DataFrame(plot_dict)    
    chart = alt.Chart(data_list_pd).mark_line().encode(
        x='FPR:Q',
        y='TPR:Q',
        color='Type:N',
    )	

    chart.save(save_path+'plot_mixture_of_distributions_figure_'+output_name[t]+'.html')

############################################################## random ccc ############################
sample_type = ["", "", ""]
sample_name = ["equidistant_randomCCC_lr105_cellCount3000", 
              "uniform_distribution_randomCCC_lr105_cellCount5000",
              "mixture_of_distribution_randomCCC_lr105_cellCount5000"]

sample_name_alternate = ["equidistant_randomCCC_knn_lr105_cellCount3000",
                "uniform_distribution_randomCCC_knn_lr105_cellCount5000",
                "mixture_of_distribution_randomCCC_threshold_distance_lr105_cellCount5000"]

output_name = ['plot_equidistant_randomCCC_figure_b', 'plot_uniform_randomCCC_figure_g', 'plot_mix_randomCCC_figure_l' ]
for t in range (0, 3):
    plot_dict = defaultdict(list)
    with gzip.open(args.data_path +  sample_name[t] +'_'+'naive_model', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("Naive"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("Naive"+sample_type[t]) #(plot_dict_temp['Type'][i])
    ###
    
    ######
    with gzip.open(args.data_path +  sample_name[t] +'_'+'rank_product_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict_temp['FPR'].append(1.0)
    plot_dict_temp['TPR'].append(1.0)
    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######
    with gzip.open(args.data_path +  sample_name_alternate[t] +'_'+'rank_product_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
    plot_dict_temp['FPR'].append(1.0)
    plot_dict_temp['TPR'].append(1.0)
    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST_alternate_cutOff"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST_alternate_cutOff"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######

    with gzip.open(args.data_path +   sample_name[t] +'_'+'rank_product_relu_10runs', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
    
#    plot_dict_temp['FPR'].append(1.0)
#    plot_dict_temp['TPR'].append(1.0)
#    plot_dict_temp['Type'].append(plot_dict_temp['Type'][1])
    
    
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append("NEST_ReLU"+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append("NEST_ReLU"+sample_type[t]) #(plot_dict_temp['Type'][i])
    
    ######
    with gzip.open(args.data_path +   sample_name[t]  +'_'+'COMMOT', 'rb') as fp: #b, b_1, a
        plot_dict_temp = pickle.load(fp) #a - [0:5]
        
    plot_dict['FPR'].append(0)
    plot_dict['TPR'].append(0)
    plot_dict['Type'].append('COMMOT'+sample_type[t]) #(plot_dict_temp['Type'][0])
    for i in range (0, len(plot_dict_temp['Type'])):
        plot_dict['FPR'].append(plot_dict_temp['FPR'][i])
        plot_dict['TPR'].append(plot_dict_temp['TPR'][i])
        plot_dict['Type'].append('COMMOT'+sample_type[t]) #(plot_dict_temp['Type'][i])
    
       

    data_list_pd = pd.DataFrame(plot_dict)    
    chart = alt.Chart(data_list_pd).mark_line().encode(
        x='FPR:Q',
        y='TPR:Q',
        color='Type:N',
    )	
    chart.save(save_path+output_name[t]+'.html')


