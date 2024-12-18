import pandas as pd
import os 
import glob 
import pickle
import gzip
import yaml

def checkDistance(dataDir):
    coords = pd.read_csv(os.path.join(dataDir, "coords.csv"), index_col = 0)
    groundTruth = glob.glob(os.path.join(dataDir, "*ground_truth_ccc"))[0]
    with gzip.open(groundTruth, 'rb') as fp:  
        _, TP, random_activation = pickle.load(fp)
    max_dist = 0
    for i, j in TP.items():
        for k in j:
            if k in coords.index:
                i_coords = coords.loc[i].values
                k_coords = coords.loc[k].values 
                distance = ((i_coords[0] - k_coords[0]) ** 2 + (i_coords[1] - k_coords[1]) ** 2) ** 0.5
                max_dist = max(max_dist, distance)
    return max_dist 

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    datasets = config["params"]["datasets"]
    data_dir = config["directories"]["data"]
    for ds in datasets:
        dataDir = os.path.join(data_dir, ds)
        print(ds)
        print(checkDistance(dataDir))

if __name__ == "__main__":
    main()



