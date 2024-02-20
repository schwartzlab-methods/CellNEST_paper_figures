import altair as alt
from vega_datasets import data
import pandas as pd
import os
from altair_themes import altairThemes

alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")

## Create a box plot to compare gene expression of DESeq2-normalized data
def gene_box_plot(data_dir, normalized_file, genes, cat, colours, fig_dir):
    df = pd.read_csv(os.path.join(data_dir, normalized_file))

    # exclude myc-amplified outlier
    outlier = "PDA_107407_Lv_M_OBp9M3_pr3"
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

        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)

        box_plot.save(os.path.join(fig_dir, fig_file_name))

def main():
    # Classical vs. Basal-like PDO analysis - compare MET gene expression
    gene_box_plot(
        data_dir="/mnt/data0/dpaliwal/software/nest-analyses/data",
        normalized_file="pdo_goi_deseq2_normalized.csv",
        genes=["MET"],
        cat="subtype",
        colours=["#BEC100", "#B0ABFF"],
        fig_dir="/mnt/data0/dpaliwal/software/nest-analyses/figures",
    )

    # Fibroblast vs. PDO analysis - compare MET & PLXNB2 epression
    gene_box_plot(
        data_dir="/mnt/data0/dpaliwal/software/nest-analyses/data",
        normalized_file="fibro.pdo_goi_deseq2_normalized.csv",
        genes=["MET", "PLXNB2"],
        cat="type",
        colours=["#393b79", "#e7ba52"],
        fig_dir="/mnt/data0/dpaliwal/software/nest-analyses/figures",
    )


if __name__ == "__main__":
    main()
