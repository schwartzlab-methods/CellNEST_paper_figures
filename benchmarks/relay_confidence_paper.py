import os
import subprocess
import glob 
import yaml
import pandas as pd
import argparse
import numpy as np

def subsample_lrdb(db: str, num_combos: int = 100) -> pd.DataFrame:
    db_df = pd.read_csv(db)
    np.random.seed(52)
    combinations = []
    pairs = list(zip(db_df["Ligand"], db_df["Receptor"]))
    for _ in range(num_combos):
        pair_1 = pairs[np.random.randint(len(pairs))]
        pair_2 = pairs[np.random.randint(len(pairs))]
        combinations.append({
            "ligand_1": pair_1[0],
            "receptor_1": pair_1[1],
            "ligand_2": pair_2[0],
            "receptor_2": pair_2[1],
            "Relay Patterns": f"{pair_1[0]}-{pair_1[1]} to {pair_2[0]}-{pair_2[1]}"
        })
    df = pd.DataFrame(combinations)

    return df

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    data_dir = config["directories"]["data"]
    nest_database = os.path.join(data_dir, "nest_database.csv")
    synthetic = subsample_lrdb(nest_database)
    synthetic_outdir = os.path.join(data_dir, "NEST_relay_validation", "synthetic")
    os.makedirs(synthetic_outdir, exist_ok = True)
    synthetic.to_csv(os.path.join(synthetic_outdir, "100_relay_count.csv"))
    result_dir = os.path.join(config["directories"]["filtered_result"], "two_hop") 
    samples = {"lymph_node": "human", "lung": "human", "exp2_D1": "human", "exp1_C1": "human", "merfish": "mouse", "synthetic": "human"}
    # samples = {"lymph_node": "human"}
    for dataset, organism in samples.items():
        try:
            input_path = glob.glob(os.path.join(data_dir, "NEST_relay_validation", dataset, "*relay_count.csv"))
            if not input_path:
                print(f"No input file found for {dataset}")
                continue
            subprocess.run(
                [
                    "python", "relay_confidence.py", 
                    "--input_path", input_path[0],
                    "--database_dir", os.path.join(data_dir, "protein_database"), 
                    "--organism", organism, 
                    "--output_path", os.path.join(result_dir, f"{dataset}_confidence_2.csv"),
                ]
            )
        except Exception as e:
            print(f"Error {e} for sample: {dataset}")


if __name__ == "__main__":
    main()

