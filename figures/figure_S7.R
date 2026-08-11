# figures/figure_S7.R
#
# S7 Fig. Allele frequency in the African sample against the two non-African
# samples, separately for common and for rare variants.
#
#   a  African against European: common variants left, rare variants right
#   b  African against East Asian: common variants left, rare variants right
#
# The point of the figure is the mass of rare variants lying on the horizontal
# axis, which are variants that segregate in the African sample and have fallen
# to zero frequency in the non-African one. Rare variants are far more likely to
# be lost this way than common ones, which is what makes the heritability gap
# widen as more of the causal signal is placed on rare variants.
#
# The panels use a smoothed scatter because there are hundreds of thousands of
# variants, and plotting them as points would produce a solid block of ink.
# Shading is proportional to the local density of variants.
#
# The allele frequencies are produced by
# simulation/08_allele_frequencies.R.
#
# Usage:  Rscript figures/figure_S7.R

source("config.R")
source(file.path(PRS_ROOT, "R", "theme.R"))

DATA_DIR <- file.path(PRS_DATA, "allele_frequencies")

common <- readRDS(file.path(DATA_DIR, "common_af.RDS"))
rare <- readRDS(file.path(DATA_DIR, "rare_af.RDS"))

# Rare variants are confined to a small corner of the unit square, so their
# panels are drawn on their own axis range rather than the full one.
RARE_AXIS_LIMIT <- 0.02

density_ramp <- colorRampPalette(c("white", unname(okabe_ito["blue"])))

#' Draw the four panels of the figure into the device that is already open.
draw_figure <- function() {
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.4, 1.2))
  on.exit(par(op), add = TRUE)

  smoothScatter(common$pafr_common, common$peur_common,
                xlab = "AFR Common", ylab = "EUR Common",
                colramp = density_ramp)
  mtext("a", side = 3, adj = 0, line = 0.8, font = 2, cex = 1.2)

  smoothScatter(rare$pafr_rare, rare$peur_rare,
                xlim = c(0, RARE_AXIS_LIMIT), ylim = c(0, RARE_AXIS_LIMIT),
                xlab = "AFR Rare", ylab = "EUR Rare",
                colramp = density_ramp)

  smoothScatter(common$pafr_common, common$pasi_common,
                xlab = "AFR Common", ylab = "ASN Common",
                colramp = density_ramp)
  mtext("b", side = 3, adj = 0, line = 0.8, font = 2, cex = 1.2)

  smoothScatter(rare$pafr_rare, rare$pasi_rare,
                xlim = c(0, RARE_AXIS_LIMIT), ylim = c(0, RARE_AXIS_LIMIT),
                xlab = "AFR Rare", ylab = "ASN Rare",
                colramp = density_ramp)
}

pdf(file.path(PRS_OUTPUT, "S7_Fig.pdf"), height = 8, width = 8)
draw_figure()
invisible(dev.off())

png(file.path(PRS_OUTPUT, "S7_Fig.png"), height = 8, width = 8,
    units = "in", res = 300)
draw_figure()
invisible(dev.off())

cat("S7 Fig written to", PRS_OUTPUT, "\n")
