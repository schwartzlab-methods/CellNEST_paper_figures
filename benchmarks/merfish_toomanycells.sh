base_dir="raw_result/misc"

for dir in "$base_dir"/*; do
    if [ -d "$dir" ]; then
        matrix_path="$dir/matrix.csv"
        labels_path="$dir/New_Cell_class_label.csv"
        output_path="$dir/toomanycells_clusters.csv"
        too-many-cells make-tree \
            --matrix-path="$matrix_path" \
            --labels-file="$labels_path" \
            --output="$dir" \
            --draw-node-number \
            --normalization=NoneNorm > "$output_path"
        echo "Processed: $dir"
    fi
done
