# scripts/figures/figure_S7.R
#
# S7 Fig. Rare variants are lost outside Africa more often than common variants.
#
# The figure itself is built by pgsbridge::figure_S7(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_S7.R

source("scripts/config.R")

figure <- figure_S7(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "S7_Fig.pdf"), figure, height = 4.5, width = 9)
ggsave(file.path(PRS_OUTPUT, "S7_Fig.png"), figure, height = 4.5, width = 9,
       dpi = 300)
cat("S7 Fig written to", PRS_OUTPUT, "\n")
