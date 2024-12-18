import yaml
import pandas as pd
import os
import altairThemes
import altair as alt
alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")


def visualize_lri(
        score_dir: str, 
        tool: str, 
        col_dict: dict
    ):
    tool_output = pd.read_csv(os.path.join(score_dir, tool, "lymph_node.csv"), index_col = 0)
    # CCL19-CCR7 in the T-cell zone (cluster 0-cluster 0) 
    sender_cluster_oi = 0
    receiver_cluster_oi = 0 
    ligand_oi = "CCL19"
    receptor_oi = "CCR7"
    # cols in tool_output 
    strength_col = col_dict["strength"][tool]
    ligand_col = col_dict["ligand"][tool]
    receptor_col = col_dict["receptor"][tool]
    sender_cluster_col = col_dict["sender_cluster"][tool]
    receiver_cluster_col = col_dict["receiver_cluster"][tool]
    tool_title = col_dict["tool_title"][tool]
    # filter for lri of interest 
    df_oi = tool_output[(tool_output[ligand_col] == ligand_oi) & (tool_output[receptor_col] == receptor_oi)]
    df_oi = df_oi.drop_duplicates()
    df_oi = df_oi.sort_values(by = strength_col, ascending = False)
    print(df_oi)
    if tool == "TWCOM":
        domain = [-1,1] # domain of tool 
        angle = -45
        width, height = 150, 150
        append = "cell type"
    else:
        domain = [tool_output[strength_col].min(), tool_output[strength_col].max()]
        print(domain)
        angle = 0
        width, height = 100, 100
        append = "cluster"
 

    heatmap = alt.Chart(df_oi).mark_rect().encode(
        x=alt.X(f'{receiver_cluster_col}:N', title=f"Receiver {append}", axis=alt.Axis(labelAngle=angle)),
        y=alt.Y(f'{sender_cluster_col}:N', title=f"Sender {append}"),
        color=alt.Color(
            f'{strength_col}:Q', 
            scale=alt.Scale(domain = domain, scheme='viridis', zero = True), 
            title=f'{strength_col}'
        )
    ).properties(
        title = f'{tool_title}',
        width = width,  
        height = height  
    )

    heatmap.save(os.path.join(score_dir, tool, "lymph_node_heatmap.html")) 

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    score_dir = config["directories"]["raw_result_user"]
    tools = config["params"]["cluster_tools"]
    col_dict = {
        "strength": {"giotto": "PI", "cellchat": "prob", "TWCOM": "Estimate"},
        "ligand": {"giotto": "ligand", "cellchat": "ligand", "TWCOM": "Ligand"},
        "receptor": {"giotto": "receptor", "cellchat": "receptor", "TWCOM": "Receptor"},
        "sender_cluster": {"giotto": "lig_cell_type", "cellchat": "source", "TWCOM": "CellType_Sender"},
        "receiver_cluster": {"giotto": "rec_cell_type", "cellchat": "target", "TWCOM": "CellType_Receiver"},
        "tool_title": {"cellchat": "CellChat v2", "giotto": "Giotto", "TWCOM": "TWCOM"}
    }
    for t in tools:
        visualize_lri(score_dir, t, col_dict)

if __name__ == "__main__":
    main()
