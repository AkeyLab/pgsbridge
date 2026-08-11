# figures/figure_S8.R
#
# S8 Fig. The percentage of variants segregating in the African sample that have
# zero frequency in a non-African sample, for common and for rare variants.
#
#   a  from African to European
#   b  from African to East Asian
#
# This puts a number on the pattern that S7 Fig shows graphically. A variant
# that is present in the African sample and absent from a population that passed
# through the out-of-Africa bottleneck contributes additive genetic variance in
# the first and none at all in the second, and rare variants are lost this way
# far more often than common ones.
#
# The allele frequencies are produced by
# simulation/08_allele_frequencies.R.
#
# Usage:  Rscript figures/figure_S8.R

source("config.R")
source(file.path(PRS_ROOT, "R", "theme.R"))

suppressMessages({
  library(ggplot2)
  library(patchwork)
})

DATA_DIR <- file.path(PRS_DATA, "allele_frequencies")

common <- readRDS(file.path(DATA_DIR, "common_af.RDS"))
rare <- readRDS(file.path(DATA_DIR, "rare_af.RDS"))

# The same two colours the heritability decomposition uses for the two variant
# classes, so a class keeps one colour across the paper.
type_pal <- unname(decomp_pal[c("common", "rare")])
names(type_pal) <- c("Common", "Rare")

#' Percentage of variants that segregate in the first population and have
#' become monomorphic in the second.
#'
#' The stored frequencies count the reference allele, which is held fixed across
#' populations, so a variant that is monomorphic in the target population has
#' frequency zero or one there. Both contribute exactly zero additive genetic
#' variance, so both count as lost.
#'
#' The denominator is every variant in the class, so the figure reports the
#' share of the whole class that is lost, not the share of those segregating in
#' the African sample.
percent_lost <- function(frequency_african, frequency_other) {
  segregating <- frequency_african > 0 & frequency_african < 1
  monomorphic <- frequency_other == 0 | frequency_other == 1
  round(100 * sum(segregating & monomorphic) / length(frequency_african), 2)
}

#' One panel: a two-bar comparison of common against rare.
loss_panel <- function(percent_common, percent_rare, target, tag) {
  df <- data.frame(Type = factor(c("Common", "Rare"),
                                 levels = c("Common", "Rare")),
                   percent = c(percent_common, percent_rare))
  ggplot(df, aes(x = Type, y = percent, fill = Type)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = percent), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = type_pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Type",
         y = sprintf("Percentage of variants reaching zero frequency in %s",
                     target),
         title = sprintf("From AFR to %s", target),
         tag = tag) +
    theme_minimal(base_size = PRS_BASE_SIZE) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.tag = element_text(face = "bold", size = 14))
}

panel_a <- loss_panel(
  percent_lost(common$pafr_common, common$peur_common),
  percent_lost(rare$pafr_rare, rare$peur_rare),
  target = "EUR", tag = "a")

panel_b <- loss_panel(
  percent_lost(common$pafr_common, common$pasi_common),
  percent_lost(rare$pafr_rare, rare$pasi_rare),
  target = "ASN", tag = "b")

figure <- panel_a | panel_b

ggsave(file.path(PRS_OUTPUT, "S8_Fig.pdf"), figure, height = 5, width = 9)
ggsave(file.path(PRS_OUTPUT, "S8_Fig.png"), figure, height = 5, width = 9,
       dpi = 300)
cat("S8 Fig written to", PRS_OUTPUT, "\n")
