#!/bin/bash
#SBATCH --job-name=twcom_en          # Job name
#SBATCH --output=twcom_en.log      # Standard output and error log
#SBATCH --account=def-gregorys
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=16              # Number of CPU cores per task
#SBATCH --time=100:00:00                 # Time limit (hh:mm:ss)
#SBATCH --mem=100000M                       # Memory limit

module load apptainer
apptainer exec --nv -B $(pwd) twcom.sif Rscript twcom_benchmark.R equidistant no

