import pandas as pd
import os

data_dir = "/mnt/data0/dpaliwal/software/nest-analyses/data"

## Subset the PDO metadata file to only contain pdos with pure Classical and Basal-like subtype signatures
def subset_pdos(input_file):
    # load the pdo metadata file
    pdo_meta = pd.read_csv(os.path.join(data_dir, input_file))

    # pure Classical and Basal-like samples defined by the Notta lab
    pure_samples = [
        "PDA_101053_Lv_M_OBp11M3_pr3",
        "PDA_097484_Lv_M_OBp10M3_pr3",
        "PDA_093429_Lv_M_OBp11M3_pr3",
        "PDA_091168_Lv_M_OBp11M3_pr3",
        "PDA_102454_Lv_M_OBp7_pr3",
        "PDA_106651_Lv_M_OBp9M3_pr3",
        "PDA_104760_Lv_M_OBp14M3_pr3",
        "PDA_100809_Lv_M_OBp8M3_pr3",
        "PDA_107872_Lv_M_OBp6M3_pr3",
        "PDA_101781_Lv_M_OBp7M3_pr3",
        "PDA_091416_Lv_M_OBp11M3_pr3",
        "PDA_091416_Lv_M_OBp11M3_pr3",
        "PDA_107407_Lv_M_OBp9M3_pr3",
        "PDA_097768_Lv_M_OBp11M3_pr3",
    ]

    pdo_pure = pdo_meta[pdo_meta["Full_name"].isin(pure_samples)]

    # only consider Basal-like & Classical samples
    pdo_pure = pdo_pure[pdo_pure["Tumor_PDO_Subtype_Notta"] != "Hybrid"]

    # convert from the Classical and Basal-like A/B subtype scheme to Classical and Basal-like
    basal_mask = pdo_pure["Tumor_PDO_Subtype_Notta"].str.contains(
        "Basal-like", case=False
    )
    classical_mask = pdo_pure["Tumor_PDO_Subtype_Notta"].str.contains(
        "Classical", case=False
    )
    pdo_pure.loc[basal_mask, "Modified_Tumor_PDO_Subtype_Notta"] = "Basal-like"
    pdo_pure.loc[classical_mask, "Modified_Tumor_PDO_Subtype_Notta"] = "Classical"

    return pdo_pure

## Subset PDO & fibroblast metadata file to only contain fibroblast samples
def subset_fibroblasts(input_file):
    fibro_pdo_meta = pd.read_csv(os.path.join(data_dir, input_file))
    fibro_meta = fibro_pdo_meta[fibro_pdo_meta["Sample Type"] == "Fibroblast"]

    return fibro_meta

## Prepare metadata files for DESeq2 input
def deseq2_meta_prep(
    input_fibro_meta, input_pdo_meta, output_pdo_file, output_fibro_pdo_file
):
    deseq2_fibro_meta = pd.DataFrame(
        {
            "sample": input_fibro_meta["Sample Name/ID"],
            "type": "fibroblast",
            "subtype": "none",
        }
    )

    deseq2_pdo_meta = pd.DataFrame(
        {
            "sample": input_pdo_meta["Full_name"],
            "type": "organoid",
            "subtype": input_pdo_meta["Modified_Tumor_PDO_Subtype_Notta"],
        }
    )

    # write the pdo metadata to a csv file for the Classical vs. Basal-like MET analysis
    deseq2_pdo_meta.to_csv(os.path.join(data_dir, output_pdo_file), index=False)

    print(deseq2_pdo_meta)
    print(len(deseq2_pdo_meta))

    # concatenate the pdo and fibroblast metadata
    deseq2_fibro_pdo_meta = pd.concat([deseq2_fibro_meta, deseq2_pdo_meta])
    # write the fibroblast & pdo metadata to a csv file for the MET and PLXNB2 analyses
    deseq2_fibro_pdo_meta.to_csv(
        os.path.join(data_dir, output_fibro_pdo_file), index=False
    )

    print(deseq2_fibro_pdo_meta)
    print(len(deseq2_fibro_pdo_meta))

    return deseq2_pdo_meta, deseq2_fibro_pdo_meta


## Prepare bulk RNA-seq counts files for DESeq2 input
def deseq2_counts_prep(input_counts_file, input_meta_df, output_counts_file):
    # load the bulk RNA-seq counts csv file
    pdo_raw = pd.read_csv(os.path.join(data_dir, input_counts_file), index_col=0)

    input_meta_df.set_index("sample", inplace=True)

    common_indices = pdo_raw.index.intersection(input_meta_df.index)
    subset_raw = pdo_raw.loc[common_indices]
    subset_raw = subset_raw.T
    subset_raw.index.name = "genes"

    subset_raw.to_csv(os.path.join(data_dir, output_counts_file))

    return subset_raw


def main():
    pdo = subset_pdos(input_file="RNA_pdo.csv")

    fibro = subset_fibroblasts(input_file="fibroblast_organoid.csv")

    pdo_df, fibro_pdo_df = deseq2_meta_prep(
        input_fibro_meta=fibro,
        input_pdo_meta=pdo,
        output_pdo_file="pdo_metadata.csv",
        output_fibro_pdo_file="fibro_pdo_metadata.csv",
    )

    ## Prepare the counts file for running DESeq2 for the Classical vs. Basal-like PDO analysis
    deseq2_counts_prep(
        input_counts_file="transposed_raw_counts.csv",
        input_meta_df=pdo_df,
        output_counts_file="deseq2_pdo_raw_counts.csv",
    )

    ## Prepare the counts file for running DESeq2 for the fibroblast vs. PDO analysis
    deseq2_counts_prep(
        input_counts_file="transposed_raw_counts.csv",
        input_meta_df=fibro_pdo_df,
        output_counts_file="deseq2_fibro_pdo_raw_counts.csv",
    )


if __name__ == "__main__":
    main()
