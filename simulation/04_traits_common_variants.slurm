#!/bin/bash
#SBATCH --job-name=traits_common
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --output=traits_common_%j.out

# Submit from the repository root:  sbatch simulation/04_traits_common_variants.slurm

module load R/4.5.1

Rscript simulation/04_traits_common_variants.R
