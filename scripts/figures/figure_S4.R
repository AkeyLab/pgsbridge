# figures/figure_S4.R
#
# S4 Fig. Transferability against the rare-variant fraction and the
# rare-to-common effect-size ratio, for the two population pairs not shown in
# Figure 3, at total heritability 0.1 in the European sample.
#
#   a  African and East Asian
#   b  European and East Asian
#
# Figure 3 shows the African and European pair at the same heritability.
#
# The underlying replicates are produced by
# simulation/06_traits_rare_and_common.R.
#
# Usage:  Rscript scripts/figures/figure_S4.R

source("scripts/config.R")

suppressMessages(library(patchwork))

DATA_DIR <- file.path(PRS_DATA, "transferability_rare_variants", "h2_0.1")
H2_TAG <- "h201"
Y_MAX <- 2.5           # shared by every panel of this figure

panel_a <- arrange_rare_variant_pair(DATA_DIR, H2_TAG, "AFR", "ASN", Y_MAX)
panel_b <- arrange_rare_variant_pair(DATA_DIR, H2_TAG, "EUR", "ASN", Y_MAX)

# wrap_elements keeps patchwork from tagging every leaf panel of the nested
# arrangements, so each population pair receives exactly one tag.
figure <- wrap_plots(lapply(list(panel_a, panel_b), wrap_elements), ncol = 1) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(PRS_OUTPUT, "S4_Fig.pdf"), figure, height = 10, width = 7)
ggsave(file.path(PRS_OUTPUT, "S4_Fig.png"), figure, height = 10, width = 7,
       dpi = 300)
cat("S4 Fig written to", PRS_OUTPUT, "\n")
