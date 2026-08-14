# simulation/03_traits_normal_effects.R
#
# Repeats the experiment of 02_traits_common_variants.R with the effect sizes
# drawn from a normal distribution instead of being equal across causal
# variants. Produces the data behind the lower panel of S2 Fig.
#
# Everything else is unchanged: the same genotypes, the same common-variant
# pool, the same 1,000 causal variants, the same heritability targets, and
# the same accuracy measurement. The only difference is how the effect sizes are
# obtained.
#
# How the effect sizes are drawn
# ------------------------------
# Effect sizes are drawn once per replicate as beta_i ~ Normal(0, 1), shared by
# all three populations, and then rescaled by a single constant so that the
# trait's heritability in the European sample takes the target value. Rescaling
# a whole effect-size vector by a constant c multiplies the additive genetic
# variance by c^2, so the constant that achieves a target heritability h2 at
# non-genetic variance tau^2 is
#
#     c = sqrt( h2/(1-h2) * tau^2 / sigma^2_unscaled ),
#
# with sigma^2_unscaled the additive genetic variance of the unscaled draw. The
# relative sizes of the effects are therefore whatever the normal draw gave, and
# only their overall scale is fixed by the heritability target.
#
# Accuracy is measured on a random subsample of `EVALUATION_SIZE` individuals
# from each population rather than on all of them, as in
# 02_traits_common_variants.R.
#
# Runtime: comparable to 02_traits_common_variants.R for a single causal-variant
# count. Use 03_traits_normal_effects.slurm.
#
# Usage:  Rscript scripts/simulation/03_traits_normal_effects.R

source("scripts/config.R")

suppressMessages(library(BEDMatrix))

set.seed(20251020)

POPULATIONS <- c(AFR = "afr_common", EUR = "eur_common", ASN = "asn_common")
CALIBRATION_POPULATION <- "EUR"

N_REPLICATES <- 1000
EVALUATION_SIZE <- 180            # individuals used to measure accuracy
INTERCEPT <- 1
NOISE_VARIANCE <- 1

N_CAUSAL_VARIANTS <- 1000
HERITABILITIES <- c(0.1, 0.4, 0.8)

genotypes <- lapply(POPULATIONS, function(stem) {
  BEDMatrix(file.path(PRS_GENOTYPES, stem))
})
names(genotypes) <- names(POPULATIONS)

n_individuals <- nrow(genotypes[[1]])
n_variants <- ncol(genotypes[[1]])
cat(sprintf("%d individuals and %d common variants per population\n",
            n_individuals, n_variants))

ORDERED_PAIRS <- list(c("EUR", "AFR"), c("EUR", "ASN"),
                      c("AFR", "EUR"), c("AFR", "ASN"),
                      c("ASN", "EUR"), c("ASN", "AFR"))

run_experiment <- function(heritability) {
  pair_labels <- vapply(ORDERED_PAIRS, paste, character(1), collapse = "-")
  values <- matrix(NA_real_, nrow = N_REPLICATES, ncol = length(ORDERED_PAIRS),
                   dimnames = list(NULL, pair_labels))

  for (replicate in seq_len(N_REPLICATES)) {
    causal <- sample_causal_variants(n_variants, N_CAUSAL_VARIANTS)
    causal_genotypes <- lapply(genotypes, function(X) X[, causal])

    frequencies <- allele_frequencies(causal_genotypes[[CALIBRATION_POPULATION]])
    unscaled <- rnorm(N_CAUSAL_VARIANTS)
    scale_factor <- sqrt(
      heritability / (1 - heritability) * NOISE_VARIANCE /
        additive_genetic_variance(unscaled, frequencies))
    beta <- scale_factor * unscaled

    traits <- lapply(causal_genotypes, simulate_trait,
                     beta = beta, alpha = INTERCEPT, tau2 = NOISE_VARIANCE)

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
             pop_comb = factor(rep(pair_labels, each = N_REPLICATES)))
}

output_dir <- file.path(PRS_DATA, "transferability_normal_effects")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (heritability in HERITABILITIES) {
  started <- Sys.time()
  result <- run_experiment(heritability)
  saveRDS(result, file.path(output_dir,
                            sprintf("transferability_m%d_h2%s.RDS",
                                    N_CAUSAL_VARIANTS, format(heritability))))
  cat(sprintf("heritability = %.1f, done in %.0f seconds\n", heritability,
              as.numeric(difftime(Sys.time(), started, units = "secs"))))
  flush.console()
}
