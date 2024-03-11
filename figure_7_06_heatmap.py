import pandas as pd
import os
import argparse
import numpy as np
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.use14corefonts"] = True

def gene_exp_prep(data_dir, pdo_file, outlier, genes_to_plot):
    """
    Log2 transform and z-score scale the DESeq2 normalized counts of PDOs.

    Parameters
    ----------
    data_dir : str
        Data directory containing deseq2 normalized PDO data 
    pdo_file : str
        File containing deseq2 normalized pdo data 
    outlier: str
        Sample ID of outlier
    genes_to_plot: List
        List of gene symbols to plot in the heatmap 

    Returns
    -------
    scaled : pd.Dataframe
        Dataframe containing normalized, log2 transformed, and z-score scaled PDO data
    """

    # read in the deseq2 normalized data of the pure pdos
    deseq2_data = pd.read_csv(os.path.join(data_dir, pdo_file), index_col=0)
    deseq2_data.index.name = "sample"
    # exclude the myc-amplified outlier
    deseq2_data = deseq2_data.drop(outlier)
    # log2 transform
    log2_transformed = np.log2(deseq2_data + 1)
    # z score scale between samples 
    mean_vals = log2_transformed.mean(axis=0)
    std_vals = log2_transformed.std(axis=0)
    scaled = (log2_transformed - mean_vals) / std_vals
    # only retain the genes of interest
    scaled = scaled.loc[:, genes_to_plot]
    # delete the deseq2 normalized data of the pure pdos (intermediate file) 
    os.remove(os.path.join(data_dir, pdo_file))

    return scaled

def metadata_prep(data_dir, meta_file, gene_exp_df):
    """
    Prepare the metadata file containing subtype assignments for the heatmap 

    Parameters
    ----------
    data_dir : str
        Data directory containing PDO metadata file 
    meta_file : str
        Metadata file with Notta & Moffit subtype assignments 
    gene_exp_df: pd.Dataframe
        Dataframe containing normalized, log2 transformed, and z-score scaled PDO data
    Returns
    -------
    meta_subset : pd.Dataframe
        Dataframe containing metadata of only pure classical and basal-like PDOs 
    """
    meta_df = pd.read_csv(os.path.join(data_dir, meta_file))
    samples = gene_exp_df.index
    meta_subset = meta_df[meta_df["Full_name"].isin(samples)]
    cols_to_retain = [
        "Full_name",
        "Biobank_ID",
        "Tumor_Subtype_Moffit",
        "Tumor_PDO_Subtype_Notta",
    ]
    meta_subset = meta_subset.loc[:, cols_to_retain]

    return meta_subset

def plot_heatmap(fig_dir, gene_exp_df, meta_df, fig_file_name):
    """
    Plot the heatmap 

    Parameters
    ----------
    fig_dir : str
        Output figure directory
    gene_exp_df : pd.Dataframe
        Dataframe containing normalized, log2 transformed, and z-score scaled PDO data
    meta_file : str
        Metadata file with Notta & Moffit subtype assignments 
    fig_file_name : str
        File name of heatmap 
    """

    meta_df.rename(columns={"Full_name": "sample"}, inplace=True)
    meta_df.set_index("sample", inplace=True)
    merged = gene_exp_df.join(meta_df)
    merged.set_index("Biobank_ID", inplace=True)

    moffit = merged.pop("Tumor_Subtype_Moffit")
    notta = merged.pop("Tumor_PDO_Subtype_Notta")

    moffit_palette = {"Basal-like": "#BEC100", "Classical": "#B0ABFF"}

    notta_palette = {
        "Basal-likeA": "#814C42",
        "Basal-likeB": "#FF7311",
        "ClassicalA": "#1C6CAB",
    }

    moffit_colours = moffit.map(moffit_palette)
    notta_colours = notta.map(notta_palette)

    notta_legend_patches = [
        mpatches.Patch(color=color, label=label)
        for label, color in notta_palette.items()
    ]
    moffit_legend_patches = [
        mpatches.Patch(color=color, label=label)
        for label, color in moffit_palette.items()
    ]

    g = sns.clustermap(
        merged,
        row_cluster=True,
        col_cluster=False,
        cmap="coolwarm",
        row_colors=[moffit_colours, notta_colours],
        annot=False,
        fmt="g",
        xticklabels=True,
        yticklabels=True,
        figsize=(10, 9),
    )

    notta_legend = g.ax_heatmap.legend(
        loc="center left",
        bbox_to_anchor=(-0.25, 0.3),
        handles=notta_legend_patches,
        frameon=False,
    )
    notta_legend.set_title(title="Notta Subtype", prop={"size": 10})
    moffit_legend = g.ax_heatmap.legend(
        loc="center left",
        bbox_to_anchor=(-0.25, 0.5),
        handles=moffit_legend_patches,
        frameon=False,
    )
    moffit_legend.set_title(title="Moffit Subtype", prop={"size": 10})
    g.ax_heatmap.add_artist(notta_legend)
    g.ax_heatmap.add_artist(moffit_legend)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")

    plt.savefig(os.path.join(fig_dir, fig_file_name))


def main():
    parser = argparse.ArgumentParser(
        "Plot a heatmap to display expression of MET, classical-associated genes, and basal-like-associated genes"
    )
    parser.add_argument(
        "--data_dir", type=str, help="Data directory containing input file"
    )
    parser.add_argument("--fig_dir", type=str, help="Directory to store heatmap")
    parser.add_argument(
        "--input_pdo_meta_file", type=str, help="input metadata for PDOs"
    )
    parser.add_argument("--outlier_id", type=str, help="Sample ID containing outlier")
    args = parser.parse_args()

    genes = [
        "MET",
        "TFF1",
        "TFF2",
        "TFF3",
        "CEACAM6",
        "LGALS4",
        "ST6GALNAC1",
        "PLA2G10",
        "TSPAN8",
        "LYZ",
        "MYO1A",
        "VSIG2",
        "CLRN3",
        "CDH17",
        "AGR3",
        "AGR2",
        "BTNL8",
        "ANXA10",
        "FAM3D",
        "CTSE",
        "REG4",
        "SPRR3",
        "SERPINB3",
        "SERPINB4",
        "VGLL1",
        "DHRS9",
        "SPRR1B",
        "KRT17",
        "KRT15",
        "TNS4",
        "SCEL",
        "KRT6A",
        "KRT7",
        "CST6",
        "LY6D",
        "FAM83A",
        "AREG",
        "FGFBP1",
        "GPR87",
        "LEMD1",
        "S100A2",
        "SLC2A1",
    ]

    scaled = gene_exp_prep(
        data_dir=args.data_dir,
        pdo_file="pdo_deseq2_normalized.csv",
        outlier=args.outlier_id,
        genes_to_plot=genes,
    )

    meta = metadata_prep(
        data_dir=args.data_dir, meta_file=args.input_pdo_meta_file, gene_exp_df=scaled
    )

    plot_heatmap(
        fig_dir=args.fig_dir,
        gene_exp_df=scaled,
        meta_df=meta,
        fig_file_name="heatmap.svg",
    )

if __name__ == "__main__":
    main()
