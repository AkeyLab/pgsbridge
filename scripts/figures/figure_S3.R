# scripts/figures/figure_S3.R
#
# S3 Fig. Average minor allele frequency in the three populations.
#
# The figure itself is built by pgsbridge::figure_S3(); this script only supplies the
# paths and writes the files, so that the figure and the code that draws it stay
# together in the package.
#
# Usage:  Rscript scripts/figures/figure_S3.R

source("scripts/config.R")

figure <- figure_S3(PRS_DATA)

ggsave(file.path(PRS_OUTPUT, "S3_Fig.pdf"), figure, height = 5, width = 7)
ggsave(file.path(PRS_OUTPUT, "S3_Fig.png"), figure, height = 5, width = 7,
       dpi = 300)
cat("S3 Fig written to", PRS_OUTPUT, "\n")
