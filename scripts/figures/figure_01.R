# scripts/figures/figure_01.R
#
# Figure 1. Theoretical determinants of intrinsic transferability.
#
# The figure itself is built by pgsbridge::figure_01(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_01.R

source("scripts/config.R")

figure <- figure_01(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "Fig1.pdf"), figure, height = 7.4, width = 9.2)
ggsave(file.path(PRS_OUTPUT, "Fig1.png"), figure, height = 7.4, width = 9.2,
       dpi = 300)
cat("Figure 1 written to", PRS_OUTPUT, "\n")
