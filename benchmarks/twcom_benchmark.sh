docker run --rm \
    -v $(pwd):/app \
    -w /app \
    twcom-image:latest \
    Rscript twcom_benchmark.R
