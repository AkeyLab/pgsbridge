# scripts/figures/figure_S2.R
#
# S2 Fig. Transferability is robust to the effect-size distribution.
#
# The figure itself is built by pgsbridge::figure_S2(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_S2.R

source("scripts/config.R")

figure <- figure_S2(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "S2_Fig.pdf"), figure, height = 9, width = 13)
ggsave(file.path(PRS_OUTPUT, "S2_Fig.png"), figure, height = 9, width = 13,
       dpi = 300)
cat("S2 Fig written to", PRS_OUTPUT, "\n")
