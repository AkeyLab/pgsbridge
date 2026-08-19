# simulation/04_traits_rare_and_common.R
#
# Simulates traits whose causal variants are part common and part rare, and
# measures transferability between every ordered pair of populations. Produces
# the data behind Figure 3 and S4 to S6 Figs.
#
# The experiment
# --------------
# Each trait carries `N_CAUSAL_VARIANTS` causal variants in total. Three things
# are swept:
#
#   the rare fraction, the share of those causal variants drawn from the rare
#     pool rather than the common pool;
#   the effect-size ratio s, the factor by which a rare causal variant's effect
#     exceeds a common one's;
#   the total heritability in the European sample, which is the outermost loop
#     and gives one output directory per level.
#
# Rare variants are given larger effects because a variant that is both common
# and strongly acting would be visible to selection; the ratio is set by hand
# here rather than derived from a selection model, so the size of its effect on
# transferability should be read as illustrative.
#
# Why this changes transferability
# --------------------------------
# The mechanism is the asymmetric loss of causal variants across the
# out-of-Africa bottleneck. A causal variant that segregates in the African
# sample can have drifted to zero frequency in a population that passed through
# the bottleneck, where it then contributes no additive genetic variance at all.
# Rare variants are the most likely to be lost this way, so the asymmetry grows
# as more of the causal signal is placed on them.
#
# Variant pools
# -------------
# Common variants are those with minor allele frequency at least
# `MAF_COMMON` in the three populations pooled, which is the panel that
# 02_prepare_genotypes.sh writes out. Rare variants are the remaining
# segregating variants, that is those whose pooled minor allele frequency lies
# strictly between zero and `MAF_COMMON`. Both pools are defined on the pooled
# sample, so a given column refers to the same variant in all three
# populations, which is what makes the cross-population comparison meaningful.
#
# Accuracy is measured on all `EVALUATION_SIZE` individuals of each population
# from each population rather than on all of them, as in
# 02_traits_common_variants.R.
#
# Runtime: many hours for the full sweep of 9 rare fractions by 4 effect-size
# ratios by 3 heritabilities. Use 04_traits_rare_and_common.slurm.
#
# Usage:  Rscript scripts/simulation/04_traits_rare_and_common.R

source("scripts/config.R")

set.seed(20230210)

CALIBRATION_POPULATION <- "EUR"

N_REPLICATES <- 200
EVALUATION_SIZE <- 1000          # every individual: nothing is fitted, so none is held out
INTERCEPT <- 0.1
NOISE_VARIANCE <- 0.1

N_CAUSAL_VARIANTS <- 1000
RARE_FRACTIONS <- seq(0.1, 0.9, by = 0.1)
EFFECT_SIZE_RATIOS <- c(1, 2, 5, 10)
HERITABILITIES <- c(0.1, 0.4, 0.8)

pools <- load_variant_pools(PRS_GENOTYPES)
n_individuals <- pools$n_individuals
cat(sprintf("%d individuals, %d common and %d rare variants\n",
            n_individuals, length(pools$common_index),
            length(pools$rare_index)))

ORDERED_PAIRS <- list(c("EUR", "AFR"), c("EUR", "ASN"),
                      c("AFR", "EUR"), c("AFR", "ASN"),
                      c("ASN", "EUR"), c("ASN", "AFR"))
PAIR_LABELS <- vapply(ORDERED_PAIRS, paste, character(1), collapse = "-")

#' Run every replicate for one combination of the swept parameters
run_experiment <- function(rare_fraction, effect_size_ratio, heritability) {
  n_rare <- round(N_CAUSAL_VARIANTS * rare_fraction)
  n_common <- N_CAUSAL_VARIANTS - n_rare

  values <- matrix(NA_real_, nrow = N_REPLICATES, ncol = length(ORDERED_PAIRS),
                   dimnames = list(NULL, PAIR_LABELS))

  for (replicate in seq_len(N_REPLICATES)) {
    common_genotypes <- draw_causal_genotypes(pools, "common", n_common)$genotypes
    rare_genotypes <- draw_causal_genotypes(pools, "rare", n_rare)$genotypes

    # The two effect sizes, chosen together so that the calibration population
    # reaches the target heritability at the requested ratio between them.
    effects <- split_effect_sizes(
      h2 = heritability, tau2 = NOISE_VARIANCE,
      p_common = allele_frequencies(common_genotypes[[CALIBRATION_POPULATION]]),
      p_rare = allele_frequencies(rare_genotypes[[CALIBRATION_POPULATION]]),
      s = effect_size_ratio)

    beta <- c(rep(effects$common, n_common), rep(effects$rare, n_rare))

    traits <- lapply(names(pools$genotypes), function(population) {
      X <- cbind(common_genotypes[[population]], rare_genotypes[[population]])
      simulate_trait(X, beta, alpha = INTERCEPT, tau2 = NOISE_VARIANCE)
    })
    names(traits) <- names(pools$genotypes)

    evaluated <- sample.int(n_individuals, EVALUATION_SIZE)

    for (k in seq_along(ORDERED_PAIRS)) {
      discovery <- ORDERED_PAIRS[[k]][1]
      target <- ORDERED_PAIRS[[k]][2]
      values[replicate, k] <- empirical_transferability(
        y_discovery = traits[[discovery]]$y[evaluated],
        y_hat_discovery = INTERCEPT + traits[[discovery]]$g[evaluated],
        y_target = traits[[target]]$y[evaluated],
        y_hat_target = INTERCEPT + traits[[target]]$g[evaluated])
    }
  }

  data.frame(emp_tau = as.vector(values),
             pop_comb = factor(rep(PAIR_LABELS, each = N_REPLICATES)))
}

for (heritability in HERITABILITIES) {
  h2_tag <- sub("0\\.", "0", format(heritability))     # 0.1 -> "01"
  output_dir <- file.path(PRS_DATA, "transferability_rare_variants",
                          sprintf("h2_%s", format(heritability)))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (rare_fraction in RARE_FRACTIONS) {
    fraction_tag <- sprintf("0%d", round(rare_fraction * 10))
    for (effect_size_ratio in EFFECT_SIZE_RATIOS) {
      started <- Sys.time()
      result <- run_experiment(rare_fraction, effect_size_ratio, heritability)
      saveRDS(result, file.path(
        output_dir,
        sprintf("transdf_mloci%d_h2%s_rare%s_s%d.RDS",
                N_CAUSAL_VARIANTS, h2_tag, fraction_tag, effect_size_ratio)))
      cat(sprintf("heritability %.1f, rare fraction %.1f, ratio %2d: %.0f s\n",
                  heritability, rare_fraction, effect_size_ratio,
                  as.numeric(difftime(Sys.time(), started, units = "secs"))))
      flush.console()
    }
  }
}
