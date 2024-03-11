import argparse
import altair as alt
from vega_datasets import data
import pandas as pd
import os
from altair_themes import altairThemes

alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")


def gene_box_plot(data_dir, normalized_file, outlier, genes, cat, colours, fig_dir):
    """
    Create a box plot to compare DESEQ2 normalized data

    Parameters
    ----------
    data_dir : str
        Data directory
    normalized_file : str
        .csv file with normalized data
    outlier : str
        Sample ID of outlier
    genes : str
        Genes to compare in the box plot
    cat : str
        Comparison of interest: subtype or type (fibroblast vs. PDO)
    colours : list
        List of colour ids for each category
    """
    genes = genes.split(",")
    df = pd.read_csv(os.path.join(data_dir, normalized_file))

    # exclude myc-amplified outlier
    df = df[df["sample"] != outlier]

    for gene in genes:
        box_plot_df = pd.DataFrame({cat: df[cat], "expression": df[gene]})

        box_plot = (
            alt.Chart(box_plot_df)
            .mark_boxplot()
            .encode(
                alt.X(f"{cat}:N").axis(labelAngle=45, title=None),
                alt.Y("expression:Q", title=f"{gene} expression"),
                alt.Color(f"{cat}:N", scale=alt.Scale(range=colours), legend=None),
            )
            .properties(
                width=70,
                height=100,
            )
        )

        fig_file_name = f"{cat}_{gene}_boxplot.html"
        box_plot.save(os.path.join(fig_dir, fig_file_name))

def main():
    parser = argparse.ArgumentParser("Create boxplots for the CCC genes of interest")
    parser.add_argument("--data_dir", type = str, help = "Data directory")
    parser.add_argument("--outlier_id", type = str, help = "Sample ID of outlier")
    parser.add_argument("--ccc_genes", type = str, help = "CCC genes of interest")
    parser.add_argument("--fig_dir", type = str, help = "Figure directory")
    args = parser.parse_args()

    # Classical vs. Basal-like PDO analysis
    gene_box_plot(
        data_dir=args.data_dir,
        normalized_file="pdo_goi_deseq2_normalized.csv",
        outlier=args.outlier_id,
        genes=args.ccc_genes,
        cat="subtype",
        colours=["#BEC100", "#B0ABFF"],
        fig_dir=args.fig_dir,
    )

if __name__ == "__main__":
    main()
