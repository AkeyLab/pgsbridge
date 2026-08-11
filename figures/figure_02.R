# figures/figure_02.R
#
# Figure 2. Transferability for three population pairs in both directions.
#
# Each ordered pair of populations gives one transferability value, and the two
# orderings of a pair are joined by a line, so the slope shows how much the
# result depends on which population supplies the score and which receives it.
# The label "ASN-AFR" means a discovery population of East Asian and a target
# population of African. The three panels are the three levels of total
# heritability in the European sample.
#
# Traits carry 1,000 common causal variants. The underlying replicates are
# produced by simulation/04_traits_common_variants.R.
#
# Usage:  Rscript figures/figure_02.R

source("config.R")
source(file.path(PRS_ROOT, "R", "summaries.R"))
source(file.path(PRS_ROOT, "R", "plots.R"))

DATA_DIR <- file.path(PRS_DATA, "transferability_common_variants")

summary_df <- read_transferability_by_heritability(
  file.path(DATA_DIR, sprintf("transferability_m1000_h2%s.RDS",
                              c("0.1", "0.4", "0.8"))))

figure <- plot_direction_pairs(summary_df)

ggsave(file.path(PRS_OUTPUT, "Fig2.pdf"), figure, height = 5, width = 12)
ggsave(file.path(PRS_OUTPUT, "Fig2.png"), figure, height = 5, width = 12,
       dpi = 300)
cat("Figure 2 written to", PRS_OUTPUT, "\n")
