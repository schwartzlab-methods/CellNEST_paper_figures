#!/bin/bash
#SBATCH -p himem
#SBATCH -t 72:00:00
#SBATCH --mem=50GB
#SBATCH -J cytospace
#SBATCH -c 4
#SBATCH -N 1

echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""
# ---------------------------------------------------------------------
cd /cluster/projects/schwartzgroup/deisha/GCN_clustering/cytospace
eval "$(conda shell.bash hook)"
conda activate /cluster/home/t114340uhn/miniconda3/envs/cytospace

cytospaceDir="/cluster/projects/schwartzgroup/deisha/software/NEST_paper_figures/benchmarks/data/cytospace"
for dir in "$cytospaceDir"/*/; do
    base_name=$(basename "$dir")
    cytospace \
        --scRNA-path "$dir/sc_counts.csv"  \
        --cell-type-path "$dir/sc_cell_types.csv"  \
        --st-path "$dir/spatial_counts.csv" \
        --coordinates-path "$dir/spatial_coordinates.csv" \
        -o "$dir"
done

echo "Job finished with exit code $? at: `date`"

