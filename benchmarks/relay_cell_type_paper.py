import os
import subprocess
import yaml
import pandas as pd

# convert annotation file from granular cell types to broader cell types 
def change_annots(
        annotation_file, 
        spot_sample
    ):
    new_cts = {
    "lymph_node": {
        "B cell": {"B GC DZ", "B GC LZ", "B Cycling", "B IFN", "B activated", "B mem", "B naive", "B plasma", "B preGC"},
        "CD4 T cell": {"T CD4", "T CD4  TfH", "T CD4  TfH GC", "T CD4  naive"},
        "CD8 T cell": {"T CD8  CD161", "T CD8  cytotoxic", "T CD8  naive"},
        "Macrophages": {"Macrophages M1", "Macrophages M2"},
        "DC": {"DC CCR7", "DC cDC1", "DC cDC2", "DC pDC"}
        }
    }

    if spot_sample in new_cts:
        df = pd.read_csv(annotation_file)
        summed_values = {}
        existing_cols_set = set()  # track columns that have been summed
        for new_key, old_vals in new_cts[spot_sample].items():
            existing_cols = [col for col in old_vals if col in df.columns]
            if existing_cols:
                summed_values[new_key] = df[existing_cols].sum(axis=1)
                existing_cols_set.update(existing_cols)
        summed_df = pd.DataFrame(summed_values)
        other_cols = df[[col for col in df.columns if col not in existing_cols_set]]
        final_df = pd.concat([other_cols, summed_df], axis=1)
        base, ext = os.path.splitext(annotation_file)
        new_annotation_file = f"{base}_broad_cell_types{ext}"
        final_df.to_csv(new_annotation_file, index=False)

        return new_annotation_file
    else:
        return annotation_file

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)

    data_dir = config["directories"]["data"]
    filtered_result_dir = os.path.join(config["directories"]["filtered_result"], "two_hop", "cell_type")
    spot_samples = ["lymph_node"]
    # Visium samples
    datasets = [
        {
            'input_dir': f'{data_dir}/NEST_relay_validation/{spot_sample}',
            'output_dir': f'{filtered_result_dir}/{spot_sample}',
            'annotation_file': f'{data_dir}/cytospace/{spot_sample}/fractional_abundances_by_spot.csv',
            'modality': 'spot'
        }
        for spot_sample in spot_samples
    ]
    # sc MERFISH brain sample
    datasets.append({
        'input_dir': f'{data_dir}/NEST_relay_validation/merfish',
        'output_dir': f'{filtered_result_dir}/merfish',
        'annotation_file': f'{data_dir}/merfish/tacco_annotation.csv',
        'modality': 'sc',
        'additional_network': 'PNOC-OPRD1 to PNOC-LPAR1'
    })
    for dataset in datasets:
        new_annotation_file = change_annots(dataset['annotation_file'], dataset['annotation_file'].split('/')[-2])
        command = [
            'python', 'relay_cell_type.py',
            '--input_dir', dataset['input_dir'],
            '--output_dir', dataset['output_dir'],
            '--annotation_file', new_annotation_file,
            '--modality', dataset['modality']
        ]
        if 'additional_network' in dataset:
            command += ['--additional_network', dataset['additional_network']]
        try:
            subprocess.run(command, check=True)
            print(f"Executed: {dataset['input_dir']} -> {dataset['output_dir']} successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {dataset['input_dir']}: {e}")

if __name__ == "__main__":
    main()
