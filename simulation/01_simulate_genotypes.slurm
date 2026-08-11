#!/bin/bash
#SBATCH --job-name=simulate_genotypes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=simulate_genotypes_%j.out

# Submit from the repository root:  sbatch simulation/01_simulate_genotypes.slurm

module load anaconda3

python simulation/01_simulate_genotypes.py --out sim.vcf --seed 20220926
