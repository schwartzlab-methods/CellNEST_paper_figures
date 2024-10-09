#docker run --rm \
#    -v $(pwd):/app \
#    -w /app \
#    cellchat-image:latest \
#    Rscript cellchat_benchmark.R

singularity exec --bind $(pwd):/app \
    /mnt/data0/dpaliwal/software/NEST_paper_figures/derivations/cellchat/cellchat.sif \
    bash -c "cd /app && Rscript cellchat_benchmark.R"
