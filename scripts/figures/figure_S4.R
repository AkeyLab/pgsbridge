# scripts/figures/figure_S4.R
#
# S4 Fig. Rare-variant sweeps, remaining pairs, heritability 0.1.
#
# The figure itself is built by pgsbridge::figure_S4(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_S4.R

source("scripts/config.R")

figure <- figure_S4(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "S4_Fig.pdf"), figure, height = 2.83, width = 7.93)
ggsave(file.path(PRS_OUTPUT, "S4_Fig.png"), figure, height = 2.83, width = 7.93,
       dpi = 300)
cat("S4 Fig written to", PRS_OUTPUT, "\n")
