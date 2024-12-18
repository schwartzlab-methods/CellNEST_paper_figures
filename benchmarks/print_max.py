import os
import pandas as pd

def tool_params(ccc_tool):
    tool_dict = {
        "cellchat": {
            "sender": "source",
            "receiver": "target",
            "ligand": "ligand",
            "receptor": "receptor",
            # "lri": "interaction_name",
            "strength": "prob"
        },
        "giotto": {
            "sender": "lig_cell_type",
            "receiver": "rec_cell_type",
            "ligand": "ligand",
            "receptor": "receptor",
            # "lri": "LR_comb",
            "strength": "PI"  # PI: significance score: log2fc * -log10(p.adj)
        },
        "TWCOM": {
            "sender": "CellType_Sender",
            "receiver": "CellType_Receiver",
            "ligand": "Ligand",
            "receptor": "Receptor",
            # "lri": "lri",
            "strength": "Estimate"
        }
    }

    return tool_dict[ccc_tool]

for root, dirs, files in os.walk(os.path.join(os.getcwd(), "raw_result")):
    for dir_name in dirs:  # loop through each subdirectory
        dir_path = os.path.join(root, dir_name)
        for file in os.listdir(dir_path):  # list files in the subdirectory
            read_file = os.path.join(dir_path, file)
            if file.endswith('.csv'):
                read_df = pd.read_csv(read_file)
                tool = dir_name
                strength_col = tool_params(tool)["strength"]
                max_value = read_df[strength_col].max()
                print(f"Maximum value for {tool} in {file}: {max_value}")
