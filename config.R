# config.R
#
# Every path used anywhere in this repository is defined here, and nowhere else.
# Change these six values to match your machine and all scripts will follow.
#
# Scripts are meant to be run from the repository root, for example
#
#     Rscript figures/figure_02.R
#
# and each one begins with `source("config.R")`.

# Root of this repository. Scripts run from the repository root, so the default
# is the working directory. Set the environment variable PRS_ROOT to override.
PRS_ROOT <- Sys.getenv("PRS_ROOT", unset = getwd())

# Pre-computed simulation results that the figure scripts read. These files are
# distributed with the repository under data/.
PRS_DATA <- file.path(PRS_ROOT, "data")

# Where figures and newly generated results are written.
PRS_OUTPUT <- file.path(PRS_ROOT, "output")

# Directory holding the PLINK genotype files produced by
# simulation/01_simulate_genotypes.py and simulation/02_prepare_genotypes.sh.
# Only the trait-simulation scripts in simulation/ need this; the figure
# scripts do not, because they read the pre-computed results in data/.
PRS_GENOTYPES <- Sys.getenv("PRS_GENOTYPES", unset = file.path(PRS_ROOT, "genotypes"))

# Additional R library path, for installations that keep packages outside the
# default location. Leave as "" if all packages are on the default path.
PRS_RLIB <- Sys.getenv("PRS_RLIB", unset = "")
if (nzchar(PRS_RLIB)) .libPaths(c(PRS_RLIB, .libPaths()))

# Number of cores for the parallel simulations. Under SLURM this follows the
# job allocation.
PRS_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))

dir.create(PRS_OUTPUT, showWarnings = FALSE, recursive = TRUE)
