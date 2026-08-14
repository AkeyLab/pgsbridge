# scripts/figures/figure_03.R
#
# Figure 3. Transferability against the rare-variant fraction.
#
# The figure itself is built by pgsbridge::figure_03(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_03.R

source("scripts/config.R")

figure <- figure_03(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "Fig3.pdf"), figure, height = 5, width = 7)
ggsave(file.path(PRS_OUTPUT, "Fig3.png"), figure, height = 5, width = 7,
       dpi = 300)
cat("Figure 3 written to", PRS_OUTPUT, "\n")
