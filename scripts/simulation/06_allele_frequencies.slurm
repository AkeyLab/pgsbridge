#!/bin/bash
#SBATCH --job-name=allele_frequencies
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48G
#SBATCH --time=04:00:00
#SBATCH --output=allele_frequencies_%j.out

# Submit from the repository root:  sbatch scripts/simulation/06_allele_frequencies.slurm

module load R/4.5.1

Rscript scripts/simulation/06_allele_frequencies.R
