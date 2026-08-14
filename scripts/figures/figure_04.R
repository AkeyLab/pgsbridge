# figures/figure_04.R
#
# Figure 4. Heritability split into the parts carried by common and by rare
# causal variants, in the three populations.
#
# Because transferability in the coalescent setting reduces to a ratio of
# heritabilities, this figure shows the quantity that drives Figures 2 and 3
# directly. Each bar is one population at one rare-variant fraction, split into
# the share of the trait variance explained by common causal variants and the
# share explained by rare ones.
#
# The layout is a grid: columns are the total heritability in the European
# sample, rows are the rare-to-common effect-size ratio, and within each cell
# the three populations are repeated across rare-variant fractions.
#
# The underlying values are produced by
# simulation/07_heritability_decomposition.R.
#
# Usage:  Rscript scripts/figures/figure_04.R

source("scripts/config.R")

suppressMessages({
  library(ggh4x)
  library(patchwork)
})

DATA_DIR <- file.path(PRS_DATA, "heritability_decomposition")

HERITABILITIES <- c("0.1", "0.8")
RATIOS <- c(2, 5, 10)

decomposition <- do.call(rbind, unlist(recursive = FALSE, lapply(
  HERITABILITIES, function(h2) lapply(RATIOS, function(ratio) {
    tag <- sub("0\\.", "0", h2)          # "0.1" -> "01", "0.8" -> "08"
    df <- readRDS(file.path(DATA_DIR,
                            sprintf("h2decomp_h2%s_s%d.RDS", tag, ratio)))
    df$h2_level <- h2
    df$ratio <- ratio
    df
  }))))

decomposition$pop <- factor(decomposition$pop, levels = c("ASN", "EUR", "AFR"))
decomposition$type <- factor(decomposition$type, levels = c("rare", "common"))
decomposition$ratio <- factor(decomposition$ratio, levels = RATIOS)
decomposition$h2_level <- factor(decomposition$h2_level, levels = HERITABILITIES)

# Constant columns that become the labels on the nested strips: fraction_title
# spans the rare-variant fractions along the top of each block, and ratio_title
# runs vertically down the right of it. The heritability heading is a plot title
# rather than a strip, added below to the top row only.
decomposition$ratio_title <- factor("Effect size ratio")
decomposition$fraction_title <- factor("Fraction of rare variants")

# The figure is a three by two grid of blocks: rows are the effect-size ratio,
# columns are the heritability in the European sample. Within a block the nine
# rare-variant fractions share one vertical scale, so the bars can be compared
# across fractions. Between blocks the scale is free, because the two
# heritability levels differ by nearly an order of magnitude and forcing one
# scale on both would flatten the 0.1 column into an unreadable strip. Each
# block is therefore built as its own plot and the six are assembled afterwards.

#' One block of the grid: three populations by nine rare-variant fractions
block <- function(heritability, ratio, show_column_title) {
  subset <- decomposition[decomposition$h2_level == heritability &
                            decomposition$ratio == ratio, ]

  plot <- ggplot(subset, aes(pop, h2, fill = type)) +
    geom_col(position = "stack", width = 0.8) +
    facet_nested(ratio_title + ratio ~ fraction_title + fraction) +
    scale_fill_manual("Variant type", values = decomp_pal,
                      breaks = c("common", "rare")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Population", y = "Heritability") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = 5.5),
          strip.text = element_text(size = 7),
          strip.background = element_rect(fill = "grey92", colour = NA),
          panel.spacing = unit(1.5, "pt"),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))

  if (show_column_title) {
    plot <- plot + ggtitle(sprintf("Heritability in EUR: %s", heritability))
  }
  plot
}

blocks <- list()
for (ratio in RATIOS) {
  for (heritability in HERITABILITIES) {
    blocks[[length(blocks) + 1]] <- block(
      heritability, ratio,
      show_column_title = identical(ratio, RATIOS[1]))
  }
}

figure <- wrap_plots(blocks, ncol = length(HERITABILITIES),
                     guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(PRS_OUTPUT, "Fig4.pdf"), figure, width = 15, height = 9)
ggsave(file.path(PRS_OUTPUT, "Fig4.png"), figure, width = 15, height = 9,
       dpi = 300)
cat("Figure 4 written to", PRS_OUTPUT, "\n")
