singularity exec --bind $(pwd):/app \
    /mnt/data0/dpaliwal/software/NEST_paper_figures/derivations/tacco/tacco.sif \
    bash -c "cd /app && python hypomap_tacco.py" \
    # >> $(pwd)/hypomap_tacco.log 2>&1 &
