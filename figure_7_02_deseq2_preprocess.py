import os
import argparse
import pandas as pd

def subset_pdos(data_dir, input_file, pure_ids):
    """
    Subset the PDO metadata file to only contain pdos with pure Classical and Basal-like subtype signatures.

    Parameters
    ----------
    data_dir : str
        Data directory containing input PDO metadata file
    input_file : str
        PDO metadata .csv file
    pure_ids: str
        Sample IDs of pure Classical and Basal-like PDOs defined by the Notta Lab

    Returns
    -------
    df_pure : pd.DataFrame
        Dataframe containing pure PDO metadata
    """
    df = pd.read_csv(os.path.join(data_dir, input_file))
    # subset pure organoids defined by the Notta lab
    pure_ids_list = pure_ids.split(",")
    df_pure = df[df["Full_name"].isin(pure_ids_list)]
    # only consider Basal-like & Classical samples
    df_pure = df_pure[df_pure["Tumor_PDO_Subtype_Notta"] != "Hybrid"]
    # convert from the Classical and Basal-like A/B subtype scheme to Classical and Basal-like subtype scheme
    
    basal_mask = df_pure["Tumor_PDO_Subtype_Notta"].str.contains(
        "Basal-like", case=False
    )
    classical_mask = df_pure["Tumor_PDO_Subtype_Notta"].str.contains(
        "Classical", case=False
    )
    df_pure.loc[basal_mask, "Modified_Tumor_PDO_Subtype_Notta"] = "Basal-like"
    df_pure.loc[classical_mask, "Modified_Tumor_PDO_Subtype_Notta"] = "Classical"
    

    return df_pure

def deseq2_meta_prep(pdo_meta):
    """
    Prepare metadata files in the format required by DESEQ2.

    Parameters
    ----------
    pdo_meta : pd.DataFrame
        Dataframe containing PDO metadata

    Returns
    -------
    pre_deseq2_pdo_meta : pd.DataFrame
        Dataframe containing PDO metadata for subtype comparison
    """
    pre_deseq2_pdo_meta = pd.DataFrame(
        {
            "sample": pdo_meta["Full_name"],
            "type": "organoid",
            "subtype": pdo_meta["Tumor_Subtype_Moffit"],
        }
    )
    print(pre_deseq2_pdo_meta)
    return pre_deseq2_pdo_meta


def deseq2_counts_prep(data_dir, input_counts_file, input_meta_df):
    """
    Prepare bulk RNA-seq counts files for DESeq2 input

    Parameters
    ----------
    data_dir : str
        Data directory containing input raw counts .csv file
    input_counts_file : str
        Input .csv file containing raw counts
    input_meta_df: pd.DataFrame
        Dataframe containing metadata for samples of interest

    Returns
    -------
    subset_raw :  pd.DataFrame
        Dataframe containing subsetted raw counts for samples of interest
    """
    raw_counts = pd.read_csv(os.path.join(data_dir, input_counts_file), index_col=0)
    input_meta_df.set_index("sample", inplace=True)
    common_indices = raw_counts.index.intersection(input_meta_df.index)
    subset_raw = raw_counts.loc[common_indices]
    subset_raw = subset_raw.T
    subset_raw.index.name = "genes"

    return subset_raw


def main():
    parser = argparse.ArgumentParser(
        "Preprocess metadata and counts files for running DESEQ2"
    )
    parser.add_argument(
        "--data_dir", type=str, help="Data directory for input and output files"
    )
    parser.add_argument(
        "--input_pdo_meta_file", type=str, help="input metadata for PDOs"
    )
    parser.add_argument(
        "--pure_ids",
        type=str,
        help="list of pure PDO sample ids defined by the Notta lab",
    )
    parser.add_argument(
        "--input_counts_file", type=str, help="input pdo raw counts .csv file"
    )
    args = parser.parse_args()

    pdo_meta = subset_pdos(
        data_dir=args.data_dir,
        input_file=args.input_pdo_meta_file,
        pure_ids=args.pure_ids,
    )

    # metadata for the Classical vs. Basal-like PDO analysis
    pre_deseq2_pdo_meta = deseq2_meta_prep(
        pdo_meta=pdo_meta
    )

    pre_deseq2_pdo_meta.to_csv(os.path.join(args.data_dir, "pre_deseq2_pdo_meta.csv"), index = False)

    # counts for running DESeq2 for the Classical vs. Basal-like PDO analysis
    pre_deseq2_pdo_raw_subset = deseq2_counts_prep(
        data_dir=args.data_dir,
        input_counts_file=args.input_counts_file,
        input_meta_df=pre_deseq2_pdo_meta,
    )
    pre_deseq2_pdo_raw_subset.to_csv(
        os.path.join(args.data_dir, "pre_deseq2_pdo_raw_subset.csv")
    )

if __name__ == "__main__":
    main()
