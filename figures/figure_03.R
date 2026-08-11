# figures/figure_03.R
#
# Figure 3. Transferability as a function of the rare-variant fraction and the
# rare-to-common effect-size ratio, for the African and European pair.
#
# Each trait carries 1,000 causal variants. The horizontal axis is the fraction
# of those that are rare, and the colour is the ratio of the effect size given
# to a rare causal variant to that given to a common one. Total heritability in
# the European sample is 0.1. The two panels are the two directions of transfer.
#
# The underlying replicates are produced by
# simulation/06_traits_rare_and_common.R.
#
# Usage:  Rscript figures/figure_03.R

source("config.R")
source(file.path(PRS_ROOT, "R", "summaries.R"))
source(file.path(PRS_ROOT, "R", "plots.R"))

DATA_DIR <- file.path(PRS_DATA, "transferability_rare_variants", "h2_0.1")
Y_MAX <- 2.25          # shared by both panels so they can be read against each other

figure <- arrange_rare_variant_pair(DATA_DIR, h2_tag = "h201",
                                    first = "AFR", second = "EUR",
                                    y_max = Y_MAX)

ggsave(file.path(PRS_OUTPUT, "Fig3.pdf"), figure, height = 5, width = 7)
ggsave(file.path(PRS_OUTPUT, "Fig3.png"), figure, height = 5, width = 7,
       dpi = 300)
cat("Figure 3 written to", PRS_OUTPUT, "\n")
