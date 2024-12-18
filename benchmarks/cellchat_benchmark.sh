singularity exec --bind $(pwd):/app \
    /mnt/data0/dpaliwal/software/NEST_paper_figures/derivations/cellchat/cellchat.sif \
    bash -c "cd /app && time Rscript cellchat_benchmark.R" \
    >> $(pwd)/cellchat_benchmark.log 2>&1 &


