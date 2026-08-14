#!/bin/bash
#SBATCH --job-name=traits_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=12:00:00
#SBATCH --output=traits_normal_%j.out

# Submit from the repository root:  sbatch simulation/05_traits_normal_effects.slurm

module load R/4.5.1

Rscript simulation/05_traits_normal_effects.R
