#docker run --rm \
#    -v $(pwd):/app \
#    -w /app \
#    josschavezf/giotto:latest \
#    Rscript giotto_benchmark.R

singularity exec --bind $(pwd):/app \
    /mnt/data0/dpaliwal/software/NEST_paper_figures/derivations/giotto/giotto.sif \
    bash -c "cd /app && Rscript giotto_benchmark.R"

