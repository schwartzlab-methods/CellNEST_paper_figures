import os
import pandas as pd
import numpy as np
from collections import defaultdict
from typing import Dict, Tuple, List
import pickle
import gzip
import yaml
import matplotlib.pyplot as plt
import altairThemes
import altair as alt
import glob
alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")

def toolParams(
        tool: str
    ) -> Dict[str, str]:
    paramsDict = {
        "cellchat": {
            "sender": "source",
            "receiver": "target",
            "ligand": "ligand",
            "receptor": "receptor",
            "strength": "prob"
        },
        "giotto": {
            "sender": "lig_cell_type",
            "receiver": "rec_cell_type",
            "ligand": "ligand",
            "receptor": "receptor",
            "strength": "PI"  # PI: significance score: log2fc * -log10(p.adj)
        },
        "TWCOM": {
            "sender": "CellType_Sender",
            "receiver": "CellType_Receiver",
            "ligand": "Ligand",
            "receptor": "Receptor",
            "strength": "Estimate"
        }
    }

    return paramsDict[tool]

def getCmat(
    cmat: np.ndarray, 
    stat: str
) -> float:
    tp = cmat[0][0]
    fp = cmat[0][1]
    fn = cmat[1][0]
    tn = cmat[1][1]
    tpr = tp / (tp + fn)
    fnr = fn / (fn + tp)
    precision = tp / (tp + fp)
    tnr = tn / (tn + fp) 
    f1 = (2 * precision * tpr) / (precision + tpr) if (precision + tpr) != 0 else np.nan
    balanced_accuracy = (tpr + tnr) / 2
    dict = {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tpr": tpr,
        "fnr": fnr,
        "precision": precision,
        "f1": f1,
        "balanced_accuracy": balanced_accuracy
    }

    return dict[stat]

def barChart(df: pd.DataFrame, metric: str):
    df_long = pd.melt(
        df, 
        id_vars = ["Tool", "Dataset"], 
        value_vars = [metric], 
        var_name = "Metric", 
        value_name = "Value"
    )
    df_long['Value'].fillna(0, inplace = True) ### filling NA with 0 (fix)
    df = df.stack()
    chart = alt.Chart(df_long).mark_bar().encode(
        x=alt.X('Tool:N', title='Tool', axis=alt.Axis(labelAngle=-45)), 
        y=alt.Y('Value:Q', title=metric),
        color='Tool:N', 
        column=alt.Column('Dataset:N',title = None, spacing = 40),
        tooltip=['Tool:N', 'Dataset:N', 'Value:Q'] 
    )

    return chart

def heatMap(df: pd.DataFrame, metric: str, mechanistic: bool):
    x_order = [
        "Equidistant without relay", "Equidistant no noise", "Equidistant low noise", "Equidistant high noise",
        "Uniform without relay", "Uniform no noise", "Uniform low noise", "Uniform high noise",
        "Mixed without relay", "Mixed no noise", "Mixed low noise", "Mixed high noise"
    ]
    if mechanistic: 
        filler = "mechanistic" 
        x_order = [f"{condition.split()[0]} {filler} {' '.join(condition.split()[1:])}" for condition in x_order]
        x_order = [condition.replace("mechanistic", "").strip() for condition in x_order]
    y_order = ["CellNEST", "CellNEST_Alternate-Cutoff", "CellNEST_ReLU", "Naive", "COMMOT", "CytoSignal", "NICHES", "CellChat v2", "Giotto", "TWCOM"]
    metric_map = {"balanced_accuracy": "Balanced accuracy", "f1": "F1 score"}
    x_order = [x for x in x_order if x in df['Dataset'].unique()]
    y_order = [y for y in y_order if y in df['Tool'].unique()]
    x_label = "Diffusion-based model" if mechanistic else "Synthetic setup"
    heatmap = alt.Chart(df).mark_rect().encode(
        x=alt.X('Dataset:N', sort = x_order, title = "Synthetic setup", axis=alt.Axis(labelAngle=-45)),  
        y=alt.Y('Tool:N', sort = y_order, title = "Method"),
        color=alt.Color(f'{metric}:Q', title = metric_map[metric], scale=alt.Scale(scheme='viridis', domain=[0, 1])), 
        tooltip=[alt.Tooltip(f'{metric}:Q', title=metric)]
    ).properties(
        width=400, 
        height=300, 
    )

    return heatmap

# format the GT
def convertGT(
        cellGT: dict, 
        clustDf: pd.DataFrame
    ):
    clustGT = defaultdict(int)
    modifCellGT = defaultdict(int)
    cellMemb = clustDf.set_index("spot")["cluster"].to_dict()
    for i, j in cellGT.items():
        if i not in cellMemb:
            continue
        i_clust = cellMemb[i]
        for l, m in j.items():
            if l not in cellMemb:
                continue
            l_clust = cellMemb[l]
            for idx in m:
                clustGT[(i_clust, l_clust, idx)] += 1
                modifCellGT[(i, l, idx)] = 1
    
    return dict(clustGT), dict(modifCellGT)

def filterToolResult(
        toolResult: pd.DataFrame, 
        strengthCol: str, 
        strengthThresh: int = 0
) -> pd.DataFrame:
    toolResult = toolResult[toolResult[strengthCol] > strengthThresh]

    return toolResult 

# load single-cell tool results 
def loadToolSc(tool: str, toolResult: pd.DataFrame, clustDf: pd.DataFrame) -> dict:
    cellSet = set(clustDf["spot"])
    toolDict = defaultdict(float)
    for _, row in toolResult.iterrows():
        sender = row["from"]
        receiver = row["to"]
        if sender and receiver in cellSet:
            if "lr" in toolResult.columns:
                idx = row["lr"]
            elif "lr pair" in toolResult.columns: # NEST and variants 
                idx = row["lr pair"]
            if "score" in toolResult.columns:
                score = row["score"]
            elif "rank" in toolResult.columns: # NEST and variants 
                rank = row["rank"]
                score = 1 / rank  # convert rank to score
            toolDict[(sender, receiver, idx)] = score
    
    return dict(toolDict) 

# load cluster tool results 
def loadToolCluster(tool: str, toolResult: pd.DataFrame, lriDb: pd.DataFrame) -> dict:
    toolDict = defaultdict(float)
    colDict = toolParams(tool)
    sender_col, receiver_col, ligand, receptor, strength = (
        colDict['sender'], colDict['receiver'], colDict['ligand'], colDict['receptor'], colDict['strength']
    )
    lriPairs = lriDb.apply(lambda row: f"{row['ligand']}-{row['receptor']}", axis = 1).tolist()
    lriMap = {pair: idx for idx, pair in enumerate(lriPairs)}
    toolResult = filterToolResult(toolResult, strength)
    for _, row in toolResult.iterrows():
        sender = row[sender_col]
        receiver = row[receiver_col]
        idx = lriMap[f"{row[ligand]}-{row[receptor]}"]
        score = row[strength]
        toolDict[(sender, receiver, idx)] = score
    
    return dict(toolDict)

# retain only top % of tool results 
def threshToolResults(toolDict: dict, thresh: int = 80) -> dict:
    allScores = list(toolDict.values())
    cutoff = np.percentile(allScores, 80) # 80th percentile
    filteredDict = {k: v for k, v in toolDict.items() if v >= cutoff}

    return filteredDict 

def confusionMatrix(modifCellGT: dict, clustGT: dict, toolDict: dict, toolType: str, method: str, clustDf: pd.DataFrame, lriDb: pd.DataFrame):
    spotMemb = clustDf.set_index("spot")["cluster"].to_dict()
    all_possible_cluster = (len(clustDf["cluster"].unique())**2)*len(lriDb)
    all_possible_sc = (len(clustDf)**2)*len(lriDb) 
    cmat = np.zeros((2, 2)) 
    if toolType == "sc":
        if method == "sc":
            for key in modifCellGT:
                if key in toolDict:
                    cmat[0][0] += 1  # TP
                else:
                    cmat[1][0] += 1  # FN
            for key in toolDict:
                if key not in modifCellGT:
                    cmat[0][1] += 1  # FP
            cmat[1,1] = all_possible_sc - cmat[0,0] - cmat[0,1] - cmat[1,0]
        elif method == "cluster":
            toolSet = {
                (spotMemb.get(s), spotMemb.get(r), idx) 
                for (s, r, idx) in toolDict
            }
            for clustKey in toolSet:
                if clustKey in clustGT:
                    cmat[0][0] += 1  # TP 
                else:
                    cmat[0][1] += 1  # FN 
            for clustKey in clustGT:
                if clustKey not in toolSet:
                    cmat[1][0] += 1  # FN
            cmat[1,1] = all_possible_cluster - cmat[0,0] - cmat[0,1] - cmat[1,0] 
    elif toolType == "cluster":
        if method == "cluster":
            for key in clustGT:
                if key in toolDict:
                    cmat[0][0] += 1  # TP
                else:
                    cmat[1][0] += 1  # FN
            for key in toolDict:
                if key not in clustGT:
                    cmat[0][1] += 1  # FP
            cmat[1,1] = all_possible_cluster - cmat[0,0] - cmat[0,1] - cmat[1,0]
        elif method == "sc":
            clustSize = clustDf['cluster'].value_counts().to_dict()
            productDict = {}
            for i in clustSize:
                for j in clustSize:
                    productDict[(i, j)] = clustSize[i] * clustSize[j]
            for (s, r, idx), predScore in toolDict.items():
                gt = clustGT.get((s, r, idx), None)
                allPairs = productDict[(s, r)]
                if gt is not None: 
                    cmat[0][0] += gt  # TP
                    cmat[0][1] += (allPairs - gt)  # FP 
                else:
                    cmat[0][1] += allPairs  # FP
            for (s, r, idx), gt in clustGT.items():
                if (s, r, idx) not in toolDict:
                    cmat[1][0] += gt  # FN
            cmat[1,1] = all_possible_sc - cmat[0,0] - cmat[0,1] - cmat[1,0] 

    return cmat

def confusionMatrixThresh(modifCellGT: dict, clustGT: dict, toolDict: dict, toolType: str, method: str, clustDf: pd.DataFrame, lriDb: pd.DataFrame, threshold: float = 0.0):
    cmat = np.zeros((2, 2))  # confusion matrix [TP, FP], [FN, TN]
    all_possible_cluster = (len(clustDf["cluster"].unique())**2)*len(lriDb)
    all_possible_sc = (len(clustDf)**2)*len(lriDb)
    if toolType == "sc":
        if method == "sc":
            # print("length of modifCellGT", len(modifCellGT))
            # print("length of tool dict", len(toolDict))
            for key in modifCellGT:
                if key not in toolDict:
                    toolDict[key] = 0  # FN score = 0 for membership check
            # print("length of tool dict after adding empty lists", len(toolDict))  
            for key in toolDict:
                toolPred = toolDict[key]
                gt = modifCellGT.get(key, None)
                predicted_class = 1 if toolPred >= threshold else 0
                if predicted_class == 1 and gt is not None:
                    cmat[0,0] += 1 # TP 
                elif predicted_class == 1 and gt is None:
                    cmat[0,1] += 1 # FP 
                elif predicted_class == 0 and gt is not None:
                    cmat[1,0] += 1 # FN
                tn = all_possible_sc - cmat[0,0] - cmat[0,1] - cmat[1,0]
                if threshold == 0:
                    cmat[0,1] += tn #FP 
                else:
                    cmat[1,1] += tn 
        elif method == "cluster":
            # print("length of clustGT", len(clustGT))
            spotMemb = clustDf.set_index("spot")["cluster"].to_dict()
            clustScore = {}
            for (s, r, idx), score in toolDict.items():
                clustS, clustR = spotMemb[s], spotMemb[r]
                clustScore.setdefault((clustS, clustR, idx), []).append(score)
            # print("length of clustScore", len(clustScore))
            for key in clustGT:
                if key not in clustScore:
                    clustScore[key] = []
            # print("length of clustScore after adding empty lists", len(clustScore))
            for key in clustScore:
                toolPred = clustScore[key]
                gt = clustGT.get(key, None)
                predicted_class = 1 if any(score >= threshold for score in toolPred) else 0
                if predicted_class == 1 and gt is not None:
                    cmat[0,0] += 1 # TP 
                elif predicted_class == 1 and gt is None:
                    cmat[0,1] += 1 # FP 
                elif predicted_class == 0 and gt is not None:
                    cmat[1,0] += 1 # FN 
                tn = all_possible_cluster - cmat[0,0] - cmat[0,1] - cmat[1,0]
                if threshold == 0:
                    cmat[0,1] += tn 
                else:
                    cmat[1,1] += tn 
    elif toolType == "cluster":
        if method == "cluster":
            for key in clustGT:
                if key not in toolDict:
                    toolDict[key] = 0  # FN score = 0 for membership check  
            for key in toolDict:
                toolPred = toolDict[key]
                gt = clustGT.get(key, None)
                predicted_class = 1 if toolPred >= threshold else 0
                if predicted_class == 1 and gt is not None:
                    cmat[0,0] += 1 # TP 
                elif predicted_class == 1 and gt is None:
                    cmat[0,1] += 1 # FP 
                elif predicted_class == 0 and gt is not None:
                    cmat[1,0] += 1 # FN 
                tn = all_possible_cluster - cmat[0,0] - cmat[0,1] - cmat[1,0]
                if threshold == 0:
                    cmat[0,1] += tn 
                else:
                    cmat[1,1] += tn 
        elif method == "sc":
            clustSize = clustDf['cluster'].value_counts().to_dict()
            productDict = {}
            for i in clustSize:
                for j in clustSize:
                    productDict[(i, j)] = clustSize[i] * clustSize[j]
            for key in clustGT:
                if key not in toolDict:
                    toolDict[key] = 0  # FN score = 0 for membership check  
            for (s, r, idx), predScore in toolDict.items():
                toolPred = toolDict[(s, r, idx)]
                gt = clustGT.get((s, r, idx), None)
                predicted_class = 1 if toolPred >= threshold else 0
                allPairs = productDict[(s,r)]
                if predicted_class == 1 and gt is not None:
                    cmat[0,0] += gt # TP
                    cmat[0][1] += (allPairs - gt)  # FP 
                elif predicted_class == 1 and gt is None:
                    cmat[0,1] += allPairs # FP 
                elif predicted_class == 0 and gt is not None:
                    cmat[1,0] += gt # FN 
                tn = all_possible_sc - cmat[0,0] - cmat[0,1] - cmat[1,0]
                if threshold == 0:
                    cmat[0,1] += tn #FP 
                else:
                    cmat[1,1] += tn 


    return cmat

def analyze(method: str, mechanistic: bool):
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    # params 
    datasets = config["params"]["datasets"]
    datasets.remove("lymph_node")
    sc_tools = config["params"]["sc_tools"]
    cluster_tools = config["params"]["cluster_tools"]
    tools = sc_tools + cluster_tools
    print(tools)
    if mechanistic == True:
        mech_title = "mechanistic"
        datasets = [dataset for dataset in datasets if "mechanistic" in dataset]
    else:
        mech_title = "non_mechanistic"
        datasets = [dataset for dataset in datasets if "mechanistic" not in dataset]
    # dirs 
    data_dir = config["directories"]["data"] 
    raw_result_dir = config["directories"]["raw_result_user"]
    filtered_result_dir = os.path.join(config["directories"]["filtered_result"], "metrics", mech_title, method)
    os.makedirs(filtered_result_dir, exist_ok = True)
    metric_data = []
    for dataset in datasets:
        LRI_gt_cells = glob.glob(os.path.join(data_dir, dataset, "*ground_truth_ccc"))[0]
        with gzip.open(LRI_gt_cells, 'rb') as fp:  
            _, cellGT, _ = pickle.load(fp)
        lriDb = pd.read_csv(glob.glob(os.path.join(data_dir, dataset, "*lrdb.csv"))[0])
        clustDf = pd.read_csv(os.path.join(data_dir, dataset, "leiden.csv"))
        clustGT, modifCellGT = convertGT(cellGT, clustDf)
        dataset_list = []
        for tool in tools:
            result_file = os.path.join(raw_result_dir, tool, f"{dataset}.csv")
            if os.path.exists(result_file):
                toolResult = pd.read_csv(result_file)
                if tool in sc_tools:
                    toolType = "sc"
                    toolDict = loadToolSc(tool, toolResult, clustDf)
                elif tool in cluster_tools:
                    toolType = "cluster"
                    toolDict = loadToolCluster(tool, toolResult, lriDb)
                cmatTotal = confusionMatrix(modifCellGT, clustGT, toolDict, toolType, method, clustDf, lriDb)
                print(method, dataset, tool)
                print(cmatTotal) 
                metrics = ["f1", "balanced_accuracy"] 
                metric_values = {metric: getCmat(cmatTotal, metric) for metric in metrics}
                metric_data.append({'Tool': tool, 'Dataset': dataset, **metric_values})
                # distribution = list(toolDict.values())
                # min, max = 0, np.max(distribution) + 0.01
                # thresholds = np.linspace(min, max, 15)
                # for threshold in thresholds:
                    # cmatThresh = confusionMatrixThresh(modifCellGT, clustGT, toolDict, toolType, method, clustDf, lriDb, threshold)
                    # print(f"threshold: {threshold}")
                    # print(cmatThresh)
                    # recall = getCmat(cmatThresh, "tpr")
                    # if threshold == max: # at max: tp + fp = 0, so precision = nan but should be 1 
                        # precision = 1
                    # else:
                        # precision = getCmat(cmatThresh, "precision")
                    # dataset_list.append({'Recall': recall, 'Precision': precision, 'Tool': tool})
        # if dataset_list:
            # dataset_df = pd.DataFrame(dataset_list)
            # print(dataset_df)
            # dataset_df["Tool"] = dataset_df["Tool"].replace({"cellchat": "CellChat v2","giotto": "Giotto"})
            # pr_chart = alt.Chart(dataset_df).mark_line().encode(
                # x=alt.X('Recall:Q', title='Recall', scale=alt.Scale(domain=[0, 1])),
                # y=alt.Y('Precision:Q', scale=alt.Scale(domain=[0, 1])),
                # color='Tool:N',
            # )
            # pr_chart.save(os.path.join(filtered_result_dir, f"pr_{dataset}.html"))
    if metric_data:
        metric_df = pd.DataFrame(metric_data)
        metric_df["Tool"] = metric_df["Tool"].replace({"cellchat": "CellChat","giotto": "Giotto"})
        print(metric_df)
        metric_df["Dataset"] = (
            metric_df["Dataset"]
            .str.replace('_', ' ', regex=False)
            .str.replace('wo', 'without ', regex=False)
            .str.replace('distribution', '', regex=False)
            .str.replace('ccc', '', regex=False)
            .str.replace('random', '', regex=False)
            .str.replace(r'\s+', ' ', regex=True)
            .str.strip()
            .str.capitalize()
        )
        # metric_df = metric_df.replace({np.nan: None}) # metric_df['F1'].fillna(0, inplace=True)
        metrics = ["f1", "balanced_accuracy"]
        for metric in metrics:
            metric_heatmap = heatMap(metric_df, metric, mechanistic)
            metric_heatmap.save(os.path.join(filtered_result_dir, f"{metric}_heatmap.html"))

if __name__ == "__main__":
    analyze(method = "cluster", mechanistic = True)
    analyze(method = "sc", mechanistic = True)
    analyze(method = "cluster", mechanistic = False)
    analyze(method = "sc", mechanistic = False)
