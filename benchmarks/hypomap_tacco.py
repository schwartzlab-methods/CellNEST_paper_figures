import tacco as tc
import scanpy as sc
import pandas as pd
import yaml
import numpy as np

def sample_adata(
    adata, 
    ct_col, 
    num_cells
    ):
    cell_types = adata.obs[ct_col].value_counts()
    num_cell_types = len(cell_types)
    sampled_indices = []
    for cell_type, count in cell_types.items():
        cell_type_indices = adata.obs[adata.obs[ct_col] == cell_type].index
        if count <= num_cells:
            sampled_indices.extend(cell_type_indices)
        else:
            sampled_indices.extend(np.random.choice(cell_type_indices, num_cells, replace = False))
    adata_subset = adata[sampled_indices, :].copy()

    return adata_subset

def annotate(
    merfish_path: str,
    sc_path: str, 
    sc_annot_col: str,
    filter_ref = False,
    subsample_ref = False,
    dataset = None,
    num_cells = 500
    ):
    merfish_adata = sc.read(merfish_path)
    sc_adata = sc.read(sc_path)
    sc_adata.X = sc_adata.raw.X.copy()
    annots_to_exclude = ["Unassigned", "NA", "Unassigned", "Doublets", "Excluded"]
    sc_adata = sc_adata[~sc_adata.obs[sc_annot_col].isin(annots_to_exclude)].copy()
    if subsample_ref == True:
        sc_adata = sample_adata(sc_adata, sc_annot_col, num_cells)
        print(sc_adata) 
    if filter_ref == True:
        sc_adata = sc_adata[sc_adata.obs["Dataset"] == dataset]
    annot_key = "tacco"
    tc.tl.annotate(
        adata = merfish_adata,
        reference = sc_adata,
        annotation_key = sc_annot_col,
        result_key = annot_key,
        assume_valid_counts = True # disable integer checking in merfish adata
    )
    probabilities = merfish_adata.obsm[annot_key].max(axis=1)
    annots = merfish_adata.obsm[annot_key].idxmax(axis = 1)
    merfish_adata.obs[annot_key] = annots 

    annot_df = pd.DataFrame({
        'annotation': annots,  # predicted cell type (from idxmax)
        'probability': probabilities  # probability (max value)
    })

    return annot_df, merfish_adata

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    data_dir = f"{config['directories']['data']}/merfish"
    merfish_path = f"{data_dir}/animal_19_merfish.h5ad"
    sc_path = f"{data_dir}/hypomap.h5ad"
    # annot_df, annotated_adata = annotate(merfish_path, sc_path, "Author_Class", subsample_ref = True)
    annot_df, annotated_adata = annotate(merfish_path, sc_path, "Author_Class", filter_ref = True, dataset = "Moffit10x")
    annot_df.to_csv(f"{data_dir}/tacco_probabilities.csv")
    annot_df['annotation'].to_csv(f"{data_dir}/tacco_annotation.csv")
    annotated_adata.write(f"{data_dir}/merfish_tacco.h5ad")

if __name__ == "__main__":
    main()
