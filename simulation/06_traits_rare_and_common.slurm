#!/bin/bash
#SBATCH --job-name=traits_rare_common
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48G
#SBATCH --time=48:00:00
#SBATCH --output=traits_rare_common_%j.out

# Submit from the repository root:  sbatch simulation/06_traits_rare_and_common.slurm

module load R/4.5.1

Rscript simulation/06_traits_rare_and_common.R
