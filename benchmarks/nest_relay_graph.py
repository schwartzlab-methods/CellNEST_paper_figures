import numpy as np
import pandas as pd
from typing import List, Literal, Dict, Tuple, Union
import networkx as nx
import yaml 
import os

def query_tf(
        grn_db: pd.DataFrame, 
        target_gene: str,
        activation_only: bool
    ) -> Dict[str, float]:
    tf_scores = {}
    target_df = grn_db[grn_db["target"] == target_gene]
    
    if activation_only:
        target_df = target_df[target_df["mode"] == 1]
    
    if target_df.empty:
        return {}
    
    for _, row in target_df.iterrows():
        tf = row["source"]
        tf_scores[tf] = row["confidence_score"]
    
    return tf_scores

def find_paths(
        grn_db: pd.DataFrame, 
        ppi_db: pd.DataFrame,
        receptor_1: str, 
        ligand_2: str,
        min_edge_weight: float, 
        return_first_path: bool,
        activation_only: bool = True
    ) -> Union[List[Dict[str, Tuple[float, float]]], None]:

    ppi_db = ppi_db[ppi_db["experimental_score"] >= min_edge_weight]
    l2_tf_dict = query_tf(grn_db, ligand_2, activation_only)
    
    if not l2_tf_dict:
        print(f"No transcription factors found for ligand 2 ({ligand_2})")
        return None

    l2_tfs = set(l2_tf_dict.keys())
    G = nx.DiGraph()

    for _, row in ppi_db.iterrows():
        G.add_edge(row["source"], row["target"], weight = row["experimental_score"])

    all_paths = []

    for tf in l2_tfs:
        queue = [(receptor_1, [receptor_1])]
        while queue:
            current_node, path = queue.pop(0)
            if current_node == tf:
                avg_score = sum(G[u][v]['weight'] for u, v in zip(path[:-1], path[1:])) / (len(path) - 1)
                grn_score = l2_tf_dict[tf]
                path_string = " -> ".join(path + [ligand_2])
                all_paths.append({path_string: (avg_score, grn_score)})
                if return_first_path:
                    print(all_paths)
                    return all_paths
            for neighbour in G[current_node]:
                if neighbour not in path:  # prevent cycles
                    queue.append((neighbour, path + [neighbour]))

    if all_paths:
        print(all_paths)
        return all_paths 
    else:
        print(f"No paths found from receptor '{receptor_1}' to transcription factors of ligand '{ligand_2}'.")
        return None 
    
def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    database_dir = os.path.join(config["directories"]["data"], "protein_database")
    grn_db = pd.read_csv(os.path.join(database_dir, "dorothea_filtered.csv"))
    ppi_db = pd.read_csv(os.path.join(database_dir, "signaling_network_human_21122021_stringdb_scored.csv"))
    find_paths(grn_db, ppi_db, "CXCR4", "CCL21", 0.3, False)

if __name__ == "__main__":
    main()


     
