# scripts/figures/figure_S5.R
#
# S5 Fig. Rare-variant sweeps, all three pairs, heritability 0.4.
#
# The figure itself is built by pgsbridge::figure_S5(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_S5.R

source("scripts/config.R")

figure <- figure_S5(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "S5_Fig.pdf"), figure, height = 15, width = 7)
ggsave(file.path(PRS_OUTPUT, "S5_Fig.png"), figure, height = 15, width = 7,
       dpi = 300)
cat("S5 Fig written to", PRS_OUTPUT, "\n")
