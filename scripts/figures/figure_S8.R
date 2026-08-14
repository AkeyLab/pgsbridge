# scripts/figures/figure_S8.R
#
# S8 Fig. African variants that are monomorphic elsewhere.
#
# The figure itself is built by pgsbridge::figure_S8(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_S8.R

source("scripts/config.R")

figure <- figure_S8(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "S8_Fig.pdf"), figure, height = 5, width = 16)
ggsave(file.path(PRS_OUTPUT, "S8_Fig.png"), figure, height = 5, width = 16,
       dpi = 170)
cat("S8 Fig written to", PRS_OUTPUT, "\n")
