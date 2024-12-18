# move results for the paper to a designated folder 
import os
import shutil 
import yaml

# spot-based deconvolution results 
def deconvolution(
        input_dir: str, 
        output_dir: str,
        samples: list
    ):
    os.makedirs(output_dir, exist_ok = True)
    files = [
        "cell_type_assignments_by_spot.png", 
        "cell_type_assignments_by_spot_jitter.png"
    ]
    for s in samples:
        for f in files:
            img = os.path.join(input_dir, s, f)
            dest = os.path.join(output_dir, f"{s}_{f}")
            if os.path.isfile(img):
                shutil.copy(img, dest)
                print(f"copied: {img} to {dest}")
            else:
                print(f"file not found: {img}")

# roc plots for synthetic benchmarks
def roc(
        input_dir: str, 
        output_dir: str
    ):
    os.makedirs(output_dir, exist_ok = True)
    for f in os.listdir(input_dir):
        if f.endswith("roc.html"):
            img = os.path.join(input_dir, f)
            dest = os.path.join(output_dir, f)
            shutil.copy(img, dest)
            print(f"copied: {img} to {dest}")

# cell type relay analysis
def cell_type(
        input_dir: str, 
        output_dir: str,
    ):
    os.makedirs(output_dir, exist_ok = True)
    # specific interactions requested by reviewer - represented by pie charts
    tissues_network = {
        "lymph_node": {
            "CCL19-CCR7 to CCL21-CXCR4",
            "CCL21-CXCR4 to CCL21-CXCR4",
            "CCL21-CCR7 to CCL21-CXCR4",
            "CCL21-CXCR4 to CCL21-CCR7",
            "CCL21-CCR7 to CCL21-CCR7",
        },
        "merfish": {
            "PNOC-OPRD1 to PNOC-LPAR1",
            "PNOC-LPAR1 to BDNF-ESR1"
        },
        "lung": {
            "PSAP-LRP1 to APOE-LRP1"
        },
        "exp2_D1": {
            "FN1-RPSA to FN1-RPSA"
        }
    }
    # tissue interactions requested by reviewer (fig 2f) - represented by bar chart 
    bar_tissues = ["lymph_node"]
    # copy pie charts
    for t, networks in tissues_network.items():
        for n in networks:
            scrubbed_n = n.replace(" ", "_").replace('"', '')
            # pie chart 
            img = os.path.join(input_dir, t, f"{scrubbed_n}_pie.html")
            dest = os.path.join(output_dir, f"{t}_{scrubbed_n}_pie.html")
            shutil.copy(img, dest) 
            print(f"copied: {img} to {dest}")
    # copy bar charts 
    for t in bar_tissues:
        img = os.path.join(input_dir, t, "bar.html")
        dest = os.path.join(output_dir, f"{t}_bar.html")
        shutil.copy(img, dest)
        print(f"copied: {img} to {dest}")

# relay network validation files
def relay_paths(
        input_dir: str, 
        output_dir: str,
    ):
    os.makedirs(output_dir, exist_ok = True)
    shutil.copytree(input_dir, output_dir, dirs_exist_ok = True)

def main():
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    base_output_dir = config["directories"]["paper_figs"]
    samples = [
        "exp1_C1", 
        "exp2_D1", 
        "lung", 
        "lymph_node",
        "merfish"
    ]
    deconvolution(
        input_dir = os.path.join(config["directories"]["data"], "cytospace"),
        output_dir = os.path.join(base_output_dir, "deconvolution"),
        samples = samples 
    )
    roc(
        input_dir = os.path.join(config["directories"]["filtered_result"], "roc"),
        output_dir = os.path.join(base_output_dir, "roc")
    )
    cell_type(
        input_dir = os.path.join(config["directories"]["filtered_result"], "two_hop", "cell_type"),
        output_dir = os.path.join(base_output_dir, "cell_type")
    )
    relay_paths(
        input_dir = os.path.join(config["directories"]["filtered_result"], "two_hop", "paths"), 
        output_dir = os.path.join(base_output_dir, "relay_paths")
    )

if __name__ == "__main__":
    main()
