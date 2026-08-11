# figures/figure_S1.R
#
# S1 Fig. Average minor allele frequency in each of the three populations.
#
# The point of the figure is a single ordering: whether the average minor allele
# frequency at shared causal sites is higher in the African sample than in the
# other two. That ordering is what sets the additive genetic variance, and
# through it the heritability ceiling that drives every later figure.
#
# The averages are taken over so many variants that their confidence intervals
# are narrower than the plotting symbols, so only the averages are drawn.
#
# The variant definition matters
# ------------------------------
# Two definitions of a common variant are in use, and they do not give the same
# ordering, so the figure shows both rather than picking one.
#
#   pooled       minor allele frequency at least MAF_COMMON in the three
#                populations pooled. This is the definition the trait
#                simulations use, because it is the panel that
#                02_prepare_genotypes.sh writes.
#   per-population   minor allele frequency at least MAF_COMMON separately in
#                every one of the three populations. This is a strictly smaller
#                set, and requiring a variant to be common everywhere selects
#                against the variation that is private to any one population.
#
# Usage:  Rscript figures/figure_S1.R

source("config.R")
source(file.path(PRS_ROOT, "R", "theme.R"))

suppressMessages({
  library(ggplot2)
})

DATA_DIR <- file.path(PRS_DATA, "allele_frequencies")
MAF_COMMON <- 0.01
POPULATION_ORDER <- c("AFR", "EUR", "ASN")

# The stored names are af_eur_001, af_afr_001 and af_asi_001. Map them onto the
# population labels used throughout, in which the East Asian sample is ASN.
POPULATION_NAMES <- c(af_eur_001 = "EUR", af_afr_001 = "AFR", af_asi_001 = "ASN")

frequencies <- readRDS(file.path(DATA_DIR, "simulated_common_variants.RDS"))
stopifnot(setequal(names(frequencies), names(POPULATION_NAMES)))
names(frequencies) <- unname(POPULATION_NAMES[names(frequencies)])

minor <- function(p) pmin(p, 1 - p)

common_everywhere <- Reduce(`&`, lapply(frequencies,
                                        function(p) minor(p) >= MAF_COMMON))

summarise_under <- function(keep, definition) {
  data.frame(
    population = factor(POPULATION_ORDER, levels = POPULATION_ORDER),
    mean_maf = vapply(POPULATION_ORDER,
                      function(p) mean(minor(frequencies[[p]][keep])),
                      numeric(1)),
    definition = definition,
    n_variants = sum(keep),
    stringsAsFactors = FALSE)
}

summary_df <- rbind(
  summarise_under(rep(TRUE, length(common_everywhere)),
                  "Common in the pooled sample"),
  summarise_under(common_everywhere,
                  "Common in every population"))

figure <- ggplot(summary_df,
                 aes(x = population, y = mean_maf,
                     colour = definition, group = definition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = unname(okabe_ito[c("vermillion", "blue")]),
                      name = "Variant definition") +
  labs(x = "Population", y = "Average minor allele frequency") +
  expand_limits(y = 0) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_minimal(base_size = PRS_BASE_SIZE) +
  theme(legend.position = "bottom")

ggsave(file.path(PRS_OUTPUT, "S1_Fig.pdf"), figure, height = 5, width = 6)
ggsave(file.path(PRS_OUTPUT, "S1_Fig.png"), figure, height = 5, width = 6,
       dpi = 300)

print(summary_df)
cat("S1 Fig written to", PRS_OUTPUT, "\n")
