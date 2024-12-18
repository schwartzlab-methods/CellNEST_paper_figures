import scanpy as sc
import pandas as pd
import anndata as ad 
from scipy.sparse import csr_matrix
import numpy as np
import altair as alt
import yaml
import os
from typing import List, Optional
import glob
import gzip 
import pickle

# for the merfish data, 
# determine the mean gene expression across all cell types

def network_cell_ids(
        relay_cell_info, 
        network: str, 
    ) -> List[str]:
    with gzip.open(relay_cell_info, "rb") as fp:
        relay = pickle.load(fp)
    instances = relay[network]
    id_set = set() 
    for instance in instances:
        for item in instance:
            id = item[0] 
            id_set.add(id)
    print(len(list(id_set)))
    return list(id_set)

def prep_adata(
        gex_csv: str,
        ct_csv: str,
        network_ids = None
) -> ad.AnnData:
    ct_df = pd.read_csv(ct_csv, index_col = 0)
    cols_to_rmv = ['Animal_ID', 'Animal_sex', 'Behavior', 'Bregma', 'Centroid_X', 'Centroid_Y', 'Cell_class', 'Neuron_cluster_ID']
    gex_df = pd.read_csv(gex_csv, index_col = 0)
    gex_df.drop(cols_to_rmv, inplace = True, axis = 1)
    gex_df.columns = gex_df.columns.str.upper() 
    adata = ad.AnnData(gex_df) # cell ids as rows, genes as cols
    common_indices = ct_df.index.intersection(gex_df.index)
    if network_ids is not None:
        print("retaining cell ids in network")
        common_indices = common_indices.intersection(network_ids)
        adata = adata[adata.obs_names.isin(common_indices), :]
        print("printing adata:", adata)
    ct_df = ct_df[ct_df.index.isin(common_indices)]
    adata.obs["cell_type"] = ct_df.loc[adata.obs_names, "annotation"]

    return adata

def mean_gex(
        adata: ad.AnnData,
        genes: list,
        ct_col: str, # cell type column in adata
    ):
    sc.pp.filter_genes(adata, min_cells = 1)
    adata = adata[:, adata.var_names.isin(genes)]
    # sc.pp.normalize_total(adata, target_sum = 1e4) 
    # sc.pp.log1p(adata)
    adata_df = adata.to_df()
    adata_df[ct_col] = adata.obs[ct_col]
    mean_gex_df = pd.DataFrame(index = adata.var_names, columns = adata.obs[ct_col].unique())
    for g in genes:
        mean_gex_df.loc[g] = adata_df.groupby(ct_col)[g].mean()
    print(mean_gex_df)

    return mean_gex_df

def plot_heatmap(
        adata: ad.AnnData,
        genes: list,
        ct_col: str, # cell type column in adata,
        network_ids: Optional[List[str]] = None
    ):
    # re-format 
    gex_df = mean_gex(adata, genes, ct_col)
    data = gex_df.reset_index().melt(id_vars=["index"], var_name="cell_type", value_name="normalized_exp") # divided by cell volume and scaled by 1000
    print(data)
    data.rename(columns={"index": "gene"}, inplace=True)
    plot = (
        alt.Chart(data).mark_rect().encode(
            x=alt.X('cell_type:N', title='Cell type', axis = alt.Axis(labelAngle=-45)),
            y=alt.Y('gene:N', title='Gene'),
            color=alt.Color(
                'normalized_exp:Q', 
                title='Normalized expression', 
                scale=alt.Scale(scheme='viridis')
            ),
            tooltip=[
            alt.Tooltip('normalized_exp:Q', title='Normalized expression')
            ]
        ).properties(
            width=400,
            height=300,
        )
    )

    return plot 

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    data_dir = os.path.join(config["directories"]["data"], "merfish")
    relay_dir = os.path.join(config["directories"]["data"], "NEST_relay_validation", "merfish")
    out_dir = os.path.join(config["directories"]["filtered_result"], "two_hop", "cell_type", "merfish")
    network_ids = network_cell_ids(glob.glob(os.path.join(relay_dir, "*pattern_distribution_cell_info"))[0], network = "PNOC-OPRD1 to PNOC-LPAR1")
    adata_network = prep_adata(os.path.join(data_dir, "merfish_animal19_bregma011.csv"), os.path.join(data_dir, "merfish_meta_Parenting_Female.csv"), network_ids)
    heatmap_network = plot_heatmap(
        adata_network,
        genes = ["PNOC", "OPRD1", "LPAR1"], 
        ct_col = "cell_type"
    )
    heatmap_network.save(os.path.join(out_dir, "gene_expression_heatmap_network_ids.html"))
    adata_all = prep_adata(os.path.join(data_dir, "merfish_animal19_bregma011.csv"), os.path.join(data_dir, "merfish_meta_Parenting_Female.csv"))
    heatmap_all = plot_heatmap(
        adata_all,
        genes = ["PNOC", "OPRD1", "LPAR1"], 
        ct_col = "cell_type"
    )
    heatmap_all.save(os.path.join(out_dir, "gene_expression_heatmap.html"))


if __name__ == "__main__":
    main()
