# scripts/figures/figure_S1.R
#
# S1 Fig. Transferability is robust to the number of causal variants.
#
# The figure itself is built by pgsbridge::figure_S1(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_S1.R

source("scripts/config.R")

figure <- figure_S1(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "S1_Fig.pdf"), figure, height = 17.1, width = 13)
ggsave(file.path(PRS_OUTPUT, "S1_Fig.png"), figure, height = 17.1, width = 13,
       dpi = 300)
cat("S1 Fig written to", PRS_OUTPUT, "\n")
