import pandas as pd 
import scanpy as sc
import os
import matplotlib
import matplotlib.pyplot as plt
import glob
import gzip
import pickle
import qnorm
import anndata as ad
import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix


# 1. https://www.nature.com/articles/s42255-022-00657-y
# 2. https://www.science.org/doi/10.1126/sciadv.adf6251 
# 3. https://www.jneurosci.org/content/43/24/4541 
# 4. https://www.science.org/doi/full/10.1126/science.aau5324?casa_token=eO1ty8GfCAgAAAAA%3AhW3A_CTsImgtmHBktanJ67EmteDFb7H-cgDdoNBFUB7kkx1YA0kPctTdUbEHK9p2IkGg3PwsuvVn_r4 
# 5. https://pmc.ncbi.nlm.nih.gov/articles/PMC3665473/
# 6. https://scbrainmap.sysneuro.net/#/celltypes 
# 7. https://www.proteinatlas.org/ENSG00000134853-PDGFRA/brain 
# 8. https://www.proteinatlas.org/ENSG00000197430-OPALIN/brain
# 9. https://www.proteinatlas.org/ENSG00000136541-ERMN/summary/rna
# 10 https://www.nature.com/articles/s41467-022-34590-1 
# 11 https://www.proteinatlas.org/ENSG00000179520-SLC17A8
# 12 https://www.nature.com/articles/s42255-022-00657-y   
# Endothelial: Igfbp7 (1), Itm2a (2), Cldn5 (2), FLT1 (3), CLDN5 (3), ABCB1 (3), VWF (3), Cd31 (5), Cd105 (5), Cd146 (5), Slc2a1 (6), Bdg (6), Pdgfra (7), Ermn (9) 
# Oligodendrocyte: Mbp (1), Mog (1), Plp1 (2), Mag (2), Slc24a2 (6), Opalin (8)  
# Inhibitory: Gad1 (4), Gad2 (4), Slc32a1 (4), Cntnap2 (6), Epha3 (6),  
# Excitatory: Slc17a6 (4), Kitl (6), Lmx1a (6), Dnaaf3 (6), Slc17a7 (10), Slc17a8 (11), Gja1 (12)
# Aqp4: marker of astrocytes and ependymal cells
# Pnoc and Oprd1: according to HPA, expression in inhibitory neurons   

cell_markers = {
    "Endothelial": [
        "Igfbp7", "Itm2a", "Cldn5", "Cldn5", "Abcb1", "Vwf", 
        "Cd31", "Cd105", "Cd146", "Slc2a1", "Bdg", "Fn1", "Pdgfra"
    ],
    "Ependymal": ["Gja1", "Aqp4", "Mlc1"], 
    "Oligodendrocyte": [
        "Mbp", "Mog", "Plp1", "Mag", "Slc24a2", "Opalin", "Ermn"
    ],
    "Inhibitory": [
        "Gad1", "Gad2", "Slc32a1", "Cntnap2", "Epha3", "Pnoc", "Oprd1",
    ],
    "Excitatory": [
        "Slc17a6", "Kitl", "Lmx1a", "Dnaaf3", "Slc17a7"
    ]
    # "NEST": ["Lpar1", "Bdnf", "Pnoc", "Oprd1"]
} 

def quantile(adata):
    # qnormed = qnorm.quantile_normalize(np.transpose(sparse.csr_matrix.toarray(adata_filtered.X)))
    if not isinstance(adata.X, np.ndarray):
        dense_data = adata.X.toarray()  # Convert sparse to dense
    else:
        dense_data = adata.X  # Use directly if already dense
    qnormed = qnorm.quantile_normalize(np.transpose(dense_data))
    qnormed = np.transpose(qnormed)
    adata.X = csr_matrix(qnormed)

    return adata 

def load_anndata(input_dir, save_dir, relay_network, quantile_normalize = False):
    matplotlib.rcParams[ 'svg.fonttype'] = 'none'
    matplotlib.rcParams['figure.figsize'] = [6,10]
    gene_exp = os.path.join(input_dir, "merfish_animal19_bregma011.csv") 
    df = pd.read_csv(gene_exp, index_col = 0)
    print(df.columns)
    metadata = df.columns[0:8]
    expression = df.drop(metadata, axis = 1)
    cell_types = df[["Cell_class"]]
    adata = sc.AnnData(expression)
    if quantile_normalize == True:
        adata = quantile(adata) # modify adata.X 
    print(adata)
    adata.obs = cell_types
    relay_cell_info = glob.glob(os.path.join(input_dir, "*pattern_distribution_cell_info"))[0]
    with gzip.open(relay_cell_info, "rb") as fp:
        relay = pickle.load(fp)
    instances = relay[relay_network]
    print(len(instances))
    cell_ids = []
    for instance in instances:
        cell_id = instance[2][0]
        cell_ids.append(cell_id)
    print(len(cell_ids))
    unique_counts = pd.Series(cell_ids).value_counts()
    for element, count in unique_counts.items():
         print(f"Element {element}: {count}")
    n_cells = len(set(cell_ids))
    adata.obs["New_Cell_class"] = adata.obs["Cell_class"]
    adata.obs.loc[cell_ids, "New_Cell_class"] = "NEST_" + adata.obs.loc[cell_ids, "Cell_class"]
    print(f"Number of unique cells: {n_cells}")
    for cell_type in list(cell_markers.keys()):
        markers = cell_markers[cell_type]
        filtered_markers = [m for m in markers if m in adata.var_names]
        if filtered_markers:
            cell_markers[cell_type] = filtered_markers
        else:
            del cell_markers[cell_type]
    adata_filtered = adata[adata.obs_names.isin(cell_ids), :].copy()
    print(adata_filtered)
    print(adata_filtered.obs["New_Cell_class"])
    matplotlib.rcParams['font.family'] = 'Arial'
    matplotlib.rcParams['font.size'] = 8
    # fig, ax = plt.subplots(figsize=(4, 4))
    sc.pl.dotplot(adata_filtered, cell_markers, groupby = 'Cell_class', show = False)
    plt.savefig(os.path.join(save_dir, "merfish_dotplot.svg"), bbox_inches = "tight") #, dpi = 300)
    adata.write_h5ad(os.path.join(input_dir, "adata.h5ad"))
    print(adata)
    return adata

def toomanycells_pp(adata, obs_cols, out_dir):
    print("adata.X shape:", adata.X.shape)
    print("adata.var_names length:", len(adata.var_names))
    print("adata.obs_names length:", len(adata.obs_names))
    print("adata.X.T shape:", adata.X.T.shape)

    expression_df = pd.DataFrame(
        adata.X.toarray().T,  
        index = adata.var_names,  
        columns = adata.obs_names 
    )
    for col in obs_cols:
        label_df = adata.obs[col].reset_index()
        label_df.columns = ["item", "label"]
        label_df.to_csv(os.path.join(out_dir, f"{col}_label.csv"), index = False)

    expression_df.to_csv(os.path.join(out_dir, "matrix.csv"), index = False)

def process(relay_network):
    # input_dir = "/mnt/data0/dpaliwal/software/NEST_paper_figures/benchmarks/relay_validation_sample_data/merfish"
    input_dir = "/mnt/data0/dpaliwal/software/NEST_paper_figures/benchmarks/data/NEST_relay_validation/merfish"
    save_dir = f"/mnt/data0/dpaliwal/software/NEST_paper_figures/benchmarks/raw_result/misc/{relay_network}"
    os.makedirs(save_dir, exist_ok = True)
    adata_merfish = load_anndata(input_dir, save_dir, relay_network)
    # toomanycells_pp(adata_merfish, ["New_Cell_class"], save_dir)

if __name__ == "__main__":
    process("PNOC-OPRD1 to PNOC-LPAR1")
    process("PNOC-LPAR1 to BDNF-ESR1")
