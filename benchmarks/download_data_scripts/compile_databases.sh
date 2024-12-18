singularity exec --bind $(pwd):/app \
    /mnt/data0/dpaliwal/software/NEST_paper_figures/derivations/omnipath/omnipath.sif \
    bash -c "cd /app && python compile_databases.py"
