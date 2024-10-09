import requests
from time import sleep
import pandas as pd
import decoupler as dc
from typing import List, Literal
import yaml 
import os

def query_stringdb(
        stringdb: str,
        string_names: str,
        escore_min: float = 0
):
    db = pd.read_csv(stringdb, sep = "\s+")
    names = pd.read_csv(string_names, sep = "\t")
    mapping = {row["#string_protein_id"]: row["preferred_name"] for idx, row in names.iterrows()}
    db["protein_1_name"] = db["protein1"].map(mapping)
    db["protein_2_name"] = db["protein2"].map(mapping)
    db["experimental"] = pd.to_numeric(db["experimental"], errors="coerce").fillna(0)
    # score is x1000 in STRING DB's text files compared to GUI
    db["experimental"] /= 1000
    # retain only interactions with experimental score > escore_min
    db = db[db["experimental"] > escore_min]
    db = db[["protein_1_name", "protein_2_name", "experimental"]]

    return db 
    
def filter_nichenet(
        nichenet_db: str,
        stringdb_df: pd.DataFrame
):
    nichenet_df = pd.read_csv(nichenet_db)
    nichenet_df = nichenet_df.drop_duplicates(subset = ['from', 'to'])
    # ensure stringdb has bidirectional relationships
    forward_pairs = set(zip(stringdb_df['protein_1_name'], stringdb_df['protein_2_name']))
    reverse_pairs = set(zip(stringdb_df['protein_2_name'], stringdb_df['protein_1_name']))
    if forward_pairs == reverse_pairs:
        print("bidirectional relationships")
        merged_df = nichenet_df.merge(
            stringdb_df,
            left_on=["from", "to"],
            right_on=["protein_1_name", "protein_2_name"],
            how="inner"
        )
        merged_df = merged_df[['from', 'to', 'experimental']].rename(
            columns={'from': 'source', 'to': 'target', 'experimental': 'experimental_score'}
        )
        duplicate_count = merged_df.duplicated(subset=['source', 'target']).sum()
        print(f"duplicate count: {duplicate_count}")

        return merged_df

def filter_dorothea(organism = "human"):
    confidence_levels = ['A','B','C','D'] # confidence level E: interactions only supported in literature 
    grn = dc.get_dorothea(
        organism = organism, 
        levels = confidence_levels
    )
    grn['mode'] = grn['weight'].apply(lambda x: 1 if x > 0 else -1) # +1: activation, -1: inhibition
    grn["confidence_score"] = grn["weight"].abs()
    grn.drop(columns = ["confidence", "weight"], inplace = True) # drop alphabetical confidence scores 
 
    return grn 
    
def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    database_dir = os.path.join(config["directories"]["data"], "protein_database")
    stringdb = os.path.join(database_dir, "stringdb", "9606.protein.physical.links.detailed.v12.0.txt")
    string_names = os.path.join(database_dir, "stringdb", "9606.protein.info.v12.0.txt")
    stringdb_df = query_stringdb(stringdb, string_names)
    nichenet_db = os.path.join(database_dir, "nichenet_v2", "signaling_network_human_21122021.csv")
    nichenet_filtered = filter_nichenet(nichenet_db, stringdb_df)
    nichenet_out = os.path.join(database_dir, f"{os.path.splitext(os.path.basename(nichenet_db))[0]}_stringdb_scored.csv")
    nichenet_filtered.to_csv(nichenet_out, index = False)
    print(f"Filtered PPI database written to: {nichenet_out}")
    dorothea_filtered = filter_dorothea()
    dorothea_out = os.path.join(database_dir, "dorothea_filtered.csv")
    dorothea_filtered.to_csv(dorothea_out, index = False)
    print(f"Filtered GRN database written to: {dorothea_out}")

if __name__ == "__main__":
    main()








