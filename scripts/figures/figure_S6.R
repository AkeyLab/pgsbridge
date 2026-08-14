# figures/figure_S6.R
#
# S6 Fig. Transferability against the rare-variant fraction and the
# rare-to-common effect-size ratio, for all three population pairs, at total
# heritability 0.8 in the European sample.
#
#   a  African and European
#   b  African and East Asian
#   c  European and East Asian
#
# Figure 3 and S4 Fig show the same three pairs at heritability 0.1.
#
# The underlying replicates are produced by
# simulation/06_traits_rare_and_common.R.
#
# Usage:  Rscript scripts/figures/figure_S6.R

source("scripts/config.R")

suppressMessages(library(patchwork))

DATA_DIR <- file.path(PRS_DATA, "transferability_rare_variants", "h2_0.8")
H2_TAG <- "h208"
Y_MAX <- 2.5           # shared by every panel of this figure

PAIRS <- list(c("AFR", "EUR"), c("AFR", "ASN"), c("EUR", "ASN"))

panels <- lapply(PAIRS, function(pair) {
  arrange_rare_variant_pair(DATA_DIR, H2_TAG, pair[1], pair[2], Y_MAX)
})

# wrap_elements keeps patchwork from tagging every leaf panel of the nested
# arrangements, so each population pair receives exactly one tag.
figure <- wrap_plots(lapply(panels, wrap_elements), ncol = 1) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(PRS_OUTPUT, "S6_Fig.pdf"), figure, height = 15, width = 7)
ggsave(file.path(PRS_OUTPUT, "S6_Fig.png"), figure, height = 15, width = 7,
       dpi = 300)
cat("S6 Fig written to", PRS_OUTPUT, "\n")
