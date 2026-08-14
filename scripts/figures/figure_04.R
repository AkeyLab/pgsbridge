# scripts/figures/figure_04.R
#
# Figure 4. Heritability split between common and rare causal variants.
#
# The figure itself is built by pgsbridge::figure_04(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_04.R

source("scripts/config.R")

figure <- figure_04(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "Fig4.pdf"), figure, height = 11.43, width = 10.28)
ggsave(file.path(PRS_OUTPUT, "Fig4.png"), figure, height = 11.43, width = 10.28,
       dpi = 300)
cat("Figure 4 written to", PRS_OUTPUT, "\n")
