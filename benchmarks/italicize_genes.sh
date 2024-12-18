base_dir="/mnt/data0/dpaliwal/software/NEST_paper_figures/benchmarks/paper_figures"
orig_dir="$base_dir/svg_orig"
out_dir="$base_dir/svg_fixed"
mkdir -p "$out_dir"
for file in "$orig_dir"/*; do
    if [[ -f "$file" ]]; then
        filename=$(basename "$file")
        sed -r 's/(font-style="normal">)([A-Z0-9]+)/\1<font-style="italic">\2<\/font-style>/g' "$file" > "$out_dir/$filename"

    fi
done
