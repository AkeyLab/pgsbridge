#!/bin/bash
#SBATCH --job-name=theory_validation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --output=theory_validation_%j.out

# Submit from the repository root:  sbatch simulation/03_theory_validation.slurm

module load R/4.5.1

Rscript simulation/03_theory_validation.R
