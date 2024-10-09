# filter raw communication results
import pandas as pd
import yaml
import os
from itertools import product

def filter_results(
    result_csv: str,
    annot_csv: str,
    sender_col: str,
    receiver_col: str,
    lri_col: str,
    strength_col: str,
    filtered_csv: str, 
    deconv: bool
):
    result_df = pd.read_csv(
        result_csv, 
        index_col = 0
    )
    annot_df = pd.read_csv(
        annot_csv
    )
    if deconv: 
        annot_df.set_index("SpotID", inplace = True)
        annot_df['max_value'] = annot_df.max(axis = 1)
        annot_df['max_count'] = (annot_df.iloc[:, :-1] == annot_df['max_value'].values.reshape(-1, 1)).sum(axis=1)
        # filter out rows where max value occurs in more than one column
        # implictly handles 0s 
        annot_df = annot_df[annot_df['max_count'] == 1]
        annot_df['cell_type'] = annot_df.idxmax(axis=1)
        annot_df.reset_index(inplace = True)
        annot_df.rename(columns = {"SpotID": "spot", "cell_type": "cluster"}, inplace = True)
        annot_df = annot_df[["cluster", "spot"]]

    result_df = result_df[[sender_col, receiver_col, lri_col, strength_col]]
    # dict mapping clusters to spots/cells 
    cluster_to_spots = annot_df.groupby('cluster')['spot'].apply(list).to_dict()
    print(cluster_to_spots)
    filtered = []
    for _, row in result_df.iterrows():
        sender_spots = cluster_to_spots.get(row[sender_col], [])
        print(len(sender_spots))
        receiver_spots = cluster_to_spots.get(row[receiver_col], [])
        print(len(receiver_spots))
        # handle missing clusters in annot_df due to filtering deconvolution results 
        if not sender_spots or not receiver_spots:
            continue
        # cartesian product of sender spots and receiver spots
        combinations = product(sender_spots, receiver_spots)
        for sender_spot, receiver_spot in combinations:
            filtered.append({
                'sender': sender_spot,
                'receiver': receiver_spot,
                'lri': row[lri_col],
                'strength': row[strength_col]
            })
    filtered_df = pd.DataFrame(filtered)
    filtered_df.to_csv(filtered_csv)

def tool_params(tool, column):
    tool_dict = {
        "cellchat": {
            "sender": "source",
            "receiver": "target",
            "lri": "interaction_name",
            "strength": "prob"
        },
        "giotto": {
            "sender": "lig_cell_type",
            "receiver": "rec_cell_type",
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

    return tool_dict[tool][column]

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    tools = ["giotto", "cellchat"]
    datasets = config["params"]["datasets"]
    datasets.reverse()
    for tool in tools:
        filtered_result_dir = os.path.join(config["directories"]["filtered_result"], tool)
        os.makedirs(filtered_result_dir, exist_ok = True)
        for dataset in datasets:
            data_dir = os.path.join(config["directories"]["data"], dataset)
            if dataset == "lymph_node" and tool == "TWCOM":
                deconv = True
                annot_csv = os.path.join(config["directories"]["data"], "cytospace", dataset, "fractional_abundances_by_spot.csv")
            else: 
                deconv = False
                annot_csv = os.path.join(data_dir, "leiden.csv")
            print(f"filtering results for {tool} and {dataset}")
            filter_results(
                result_csv = os.path.join(config["directories"]["raw_result"], tool, f"{dataset}.csv"),
                annot_csv = annot_csv,
                sender_col = tool_params(tool, "sender"),
                receiver_col = tool_params(tool, "receiver"),
                lri_col = tool_params(tool, "lri"),
                strength_col = tool_params(tool, "strength"),
                filtered_csv = os.path.join(filtered_result_dir, f"{dataset}.csv"),
                deconv = deconv
            )
if __name__ == "__main__": 
    main()
