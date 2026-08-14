# figures/figure_S2.R
#
# S2 Fig. Transferability for three population pairs in both directions, with
# fixed and with normally distributed effect sizes.
#
# The upper panel repeats Figure 2, in which every causal variant is given the
# same effect size. The lower panel repeats the same experiment with the effect
# sizes drawn from a normal distribution instead. The pattern is unchanged,
# which shows that the result does not depend on the effect sizes being equal.
#
# The underlying replicates are produced by
# simulation/04_traits_common_variants.R for the upper panel and
# simulation/05_traits_normal_effects.R for the lower one.
#
# Usage:  Rscript scripts/figures/figure_S2.R

source("scripts/config.R")

suppressMessages(library(patchwork))

HERITABILITIES <- c("0.1", "0.4", "0.8")

panel_fixed <- plot_direction_pairs(read_transferability_by_heritability(
  file.path(PRS_DATA, "transferability_common_variants",
            sprintf("transferability_m1000_h2%s.RDS", HERITABILITIES))))

panel_normal <- plot_direction_pairs(read_transferability_by_heritability(
  file.path(PRS_DATA, "transferability_normal_effects",
            sprintf("transferability_m1000_h2%s.RDS", HERITABILITIES))))

figure <- (panel_fixed / panel_normal) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(PRS_OUTPUT, "S2_Fig.pdf"), figure, height = 9, width = 13)
ggsave(file.path(PRS_OUTPUT, "S2_Fig.png"), figure, height = 9, width = 13,
       dpi = 300)
cat("S2 Fig written to", PRS_OUTPUT, "\n")
