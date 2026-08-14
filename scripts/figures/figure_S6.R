# scripts/figures/figure_S6.R
#
# S6 Fig. Rare-variant sweeps, all three pairs, heritability 0.8.
#
# The figure itself is built by pgsbridge::figure_S6(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_S6.R

source("scripts/config.R")

figure <- figure_S6(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "S6_Fig.pdf"), figure, height = 10, width = 14)
ggsave(file.path(PRS_OUTPUT, "S6_Fig.png"), figure, height = 10, width = 14,
       dpi = 170)
cat("S6 Fig written to", PRS_OUTPUT, "\n")
