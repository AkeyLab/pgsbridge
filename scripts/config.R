# Paths used by the analysis scripts.  Change these to match your machine and
# every script follows.  Scripts are run from the repository root, for example
#
#     Rscript scripts/figures/figure_02.R
#
# The package itself needs none of this; these are only the locations of the
# pre-computed results and of the genotypes.

PRS_ROOT <- Sys.getenv("PRS_ROOT", unset = getwd())

# Pre-computed simulation results, distributed with the repository.
PRS_DATA <- file.path(PRS_ROOT, "data")

# Where figures and newly generated results are written.
PRS_OUTPUT <- file.path(PRS_ROOT, "output")

# PLINK genotypes written by scripts/simulation/02_prepare_genotypes.sh.  Only
# the simulation scripts need these; the figure scripts do not.
PRS_GENOTYPES <- Sys.getenv("PRS_GENOTYPES",
                            unset = file.path(PRS_ROOT, "genotypes"))

# Cores for the parallel simulations; follows the SLURM allocation when present.
PRS_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))

dir.create(PRS_OUTPUT, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(pgsbridge)
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
})
