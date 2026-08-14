# scripts/figures/figure_S7.R
#
# S7 Fig. Allele frequency in African against non-African samples.
#
# The figure is drawn by pgsbridge::figure_S7(), which uses base graphics and so
# writes onto whichever device is open rather than returning an object.
#
# Usage:  Rscript scripts/figures/figure_S7.R

source("scripts/config.R")

pdf(file.path(PRS_OUTPUT, "S7_Fig.pdf"), height = 8, width = 8)
figure_S7(PRS_DATA)
invisible(dev.off())

png(file.path(PRS_OUTPUT, "S7_Fig.png"), height = 8, width = 8,
    units = "in", res = 300)
figure_S7(PRS_DATA)
invisible(dev.off())

cat("S7 Fig written to", PRS_OUTPUT, "\n")
