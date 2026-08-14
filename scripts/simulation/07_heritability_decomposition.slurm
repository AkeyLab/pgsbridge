#!/bin/bash
#SBATCH --job-name=heritability_decomp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48G
#SBATCH --time=24:00:00
#SBATCH --output=heritability_decomp_%j.out

# Submit from the repository root:  sbatch simulation/07_heritability_decomposition.slurm

module load R/4.5.1

Rscript simulation/07_heritability_decomposition.R
