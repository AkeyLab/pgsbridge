# scripts/figures/figure_02.R
#
# Figure 2. Transferability for three population pairs in both directions.
#
# The figure itself is built by pgsbridge::figure_02(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_02.R

source("scripts/config.R")

figure <- figure_02(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "Fig2.pdf"), figure, height = 5.4, width = 12)
ggsave(file.path(PRS_OUTPUT, "Fig2.png"), figure, height = 5.4, width = 12,
       dpi = 300)
cat("Figure 2 written to", PRS_OUTPUT, "\n")
