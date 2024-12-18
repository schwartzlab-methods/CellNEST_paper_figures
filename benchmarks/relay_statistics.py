# mann-whitney u test to compare actual scores to null scores 
# replace None with 0? 
import numpy as np
from scipy import stats
import pandas as pd
import scikit_posthocs as sp
import yaml
import os
from statsmodels.stats.multitest import multipletests
import altair as alt
import altairThemes
alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")

def mwu(
        df_actual: pd.DataFrame,
        df_null: pd.DataFrame,
        dataset: str
    ) -> pd.DataFrame:
    actual_scores = df_actual.dropna(subset=['combined_score'])["combined_score"].to_numpy().astype(float)
    null_scores = df_null.dropna(subset=['combined_score'])["combined_score"].to_numpy().astype(float)
    u_stat, p_val = stats.mannwhitneyu(actual_scores, null_scores, alternative='greater')
    # test_results = pg.mannwhitneyu(actual_scores, null_scores, alternative='greater')
    n1 = len(actual_scores)
    n2 = len(null_scores)
    median_actual = np.median(actual_scores)
    median_null = np.median(null_scores)
    effect_size = u_stat / (n1 * n2) if n1 * n2 > 0 else 0
    df = pd.DataFrame(
        [[dataset, u_stat, p_val, effect_size, median_actual, median_null, n1, n2]],
        columns = ["dataset", "U_stat", "p_value", "effect_size", "median_actual", "median_null", "size_actual", "size_null"]
    )

    return df

def kW_Dunn(
        df: pd.DataFrame, 
        pAdj: str = 'fdr_bh'
    ):
    categories = df['dataset'].unique()
    groups = [df[df['dataset'] == category]['score'].dropna().tolist() for category in categories]
    kwResults = stats.kruskal(*groups)
    dunnResults = sp.posthoc_dunn(groups, pAdj)
    print(categories)
    print("Kruskal-Wallis")
    print(kwResults)
    print("Dunn")
    print(dunnResults)
    print(dunnResults.to_string(formatters={'p-value': '{:.3g}'.format}))

    return categories, kwResults, dunnResults


def prepare_combined_data(datasets, result_dir):
    combined_data = []
    df_null = pd.read_csv(f"{result_dir}/synthetic_confidence.csv")
    null_scores = df_null.dropna(subset=['combined_score'])["combined_score"].to_numpy().astype(float)
    null_labels = ['null'] * len(null_scores)
    combined_data.append(pd.DataFrame({
        'score': null_scores,
        'dataset': null_labels
    }))
    for dataset in datasets:
        df_actual = pd.read_csv(f"{result_dir}/{dataset}_confidence.csv")
        actual_scores = df_actual.dropna(subset=['combined_score'])["combined_score"].to_numpy().astype(float)
        actual_labels = [dataset] * len(actual_scores) 
        combined_data.append(pd.DataFrame({
            'score': actual_scores,
            'dataset': actual_labels
        }))
    combined_df = pd.concat(combined_data, ignore_index=True)
    kW_Dunn(combined_df)

    return combined_df



def create_boxplot(df):
    dataset_rename = {
        "lymph_node": "Lymph node",
        "lung": "LUAD",
        "exp1_C1": "PDAC_140694",
        "exp2_D1": "PDAC_64630",
        # "merfish": "Female mouse brain",
        "null": "Null"
    }
    df["dataset"] = df["dataset"].map(dataset_rename)
    dataset_order = ["Null", "Lymph node", "LUAD", "PDAC_140694", "PDAC_64630"]

    set1 = altairThemes.get_colour_scheme("Set1", len(df["dataset"].unique()))

    box_plot = alt.Chart(df).mark_boxplot().encode(
        x = alt.X("dataset:N", title = None, axis = alt.Axis(labelAngle=-45), sort = dataset_order),
        y = alt.Y("score:Q", title = "Confidence score", axis = alt.Axis(values = [1, 0.1, 0.01, 0.001]), scale = alt.Scale(type = 'log')),   
        color=alt.Color('dataset:N', scale = alt.Scale(range = set1), legend = None), 
    ).properties(
        # width = 160, 
        # height = 200
    )

    return box_plot


def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    result_dir = os.path.join(config["directories"]["filtered_result"], "two_hop", "paths")
    datasets = ["lymph_node", "lung", "exp1_C1", "exp2_D1"] # , "merfish"]
    n_subsample = 100
    all_results = []
    for d in datasets:
        df_actual = pd.read_csv(os.path.join(result_dir, f"{d}_confidence.csv"))
        df_null = pd.read_csv(os.path.join(result_dir, f"synthetic_confidence.csv"))
        result_df = mwu(df_actual, df_null, d)
        all_results.append(result_df)
    all_results_df = pd.concat(all_results, ignore_index = True)
    all_results_df["p_value"] = all_results_df["p_value"].astype(float)
    all_results_df["p_value_fdr_bh"] = multipletests(all_results_df["p_value"], method="fdr_bh")[1] 
    cols_sci_not = ["p_value", "p_value_fdr_bh", "effect_size", "U_stat"]
    for col in cols_sci_not:
        all_results_df[col] = all_results_df[col].apply(lambda x: format(x, '.3g'))
    output_file = os.path.join(result_dir, f"actual_vs_shuffled_stats_{n_subsample}.csv")
    all_results_df.to_csv(output_file, index = False)
    # boxplot 
    combined_df = prepare_combined_data(datasets, result_dir)
    box_plot = create_boxplot(combined_df)
    box_plot.save(os.path.join(result_dir, f"actual_vs_shuffled_stats_{n_subsample}.html"))

if __name__ == "__main__":
    main()

