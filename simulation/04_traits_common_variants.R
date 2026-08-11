# simulation/04_traits_common_variants.R
#
# Simulates traits with common causal variants only, on the coalescent
# genotypes, and measures transferability between every ordered pair of
# populations. Produces the data behind Figure 2 and S2 Fig.
#
# The experiment
# --------------
# Genotypes are fixed; each replicate redraws the causal variants and the
# non-genetic deviations. Within a replicate:
#
#   1. draw `m` causal variants, spread evenly across the common-variant pool;
#   2. choose the single effect size shared by all of them so that the trait's
#      heritability in the European sample takes the target value;
#   3. simulate the trait in all three populations using those same effect
#      sizes, so the only thing that differs between populations is their
#      allele-frequency spectrum;
#   4. measure the optimal predictor's accuracy in each population and form the
#      transferability ratio for each ordered pair.
#
# Effect sizes are the same in every population by construction. This is part of
# the best-case idealisation, and it is what isolates the allele-frequency
# channel: any difference in accuracy between populations must come from their
# allele frequencies, because nothing else about the genetic architecture
# differs.
#
# The heritability target is set in the European sample and the other two
# populations take whatever heritability the same effect sizes produce in them.
#
# Accuracy is measured on a random subsample of `EVALUATION_SIZE` individuals
# from each population rather than on all of them.
#
# Runtime: several hours for the full sweep. Use 04_traits_common_variants.slurm.
#
# Usage:  Rscript simulation/04_traits_common_variants.R

source("config.R")
source(file.path(PRS_ROOT, "R", "theory.R"))
source(file.path(PRS_ROOT, "R", "traits.R"))

suppressMessages(library(BEDMatrix))

set.seed(20220926)

POPULATIONS <- c(AFR = "afr_common", EUR = "eur_common", ASN = "asn_common")
CALIBRATION_POPULATION <- "EUR"   # the population whose heritability is fixed

N_REPLICATES <- 1000
EVALUATION_SIZE <- 180            # individuals used to measure accuracy
INTERCEPT <- 1
NOISE_VARIANCE <- 1               # tau^2; the heritability target fixes sigma^2

CAUSAL_VARIANT_COUNTS <- c(20, 50, 100, 1000)
HERITABILITIES <- c(0.1, 0.4, 0.8)

genotypes <- lapply(POPULATIONS, function(stem) {
  BEDMatrix(file.path(PRS_GENOTYPES, stem))
})
names(genotypes) <- names(POPULATIONS)

n_individuals <- nrow(genotypes[[1]])
n_variants <- ncol(genotypes[[1]])
cat(sprintf("%d individuals and %d common variants per population\n",
            n_individuals, n_variants))

# The six ordered pairs, written as discovery then target.
ORDERED_PAIRS <- list(c("EUR", "AFR"), c("EUR", "ASN"),
                      c("AFR", "EUR"), c("AFR", "ASN"),
                      c("ASN", "EUR"), c("ASN", "AFR"))

#' Run every replicate for one setting of the causal-variant count and the
#' heritability.
#'
#' @return A data frame with one row per replicate per ordered pair, holding the
#'   transferability value and the pair label.
run_experiment <- function(m_causal, heritability) {
  pair_labels <- vapply(ORDERED_PAIRS, paste, character(1), collapse = "-")
  values <- matrix(NA_real_, nrow = N_REPLICATES, ncol = length(ORDERED_PAIRS),
                   dimnames = list(NULL, pair_labels))

  for (replicate in seq_len(N_REPLICATES)) {
    causal <- sample_causal_variants(n_variants, m_causal)
    causal_genotypes <- lapply(genotypes, function(X) X[, causal])

    # One effect size, shared by every causal variant and every population,
    # chosen so that the calibration population reaches the target heritability.
    frequencies <- allele_frequencies(causal_genotypes[[CALIBRATION_POPULATION]])
    beta <- rep(shared_effect_size(heritability, NOISE_VARIANCE, frequencies),
                m_causal)

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

output_dir <- file.path(PRS_DATA, "transferability_common_variants")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (m_causal in CAUSAL_VARIANT_COUNTS) {
  for (heritability in HERITABILITIES) {
    started <- Sys.time()
    result <- run_experiment(m_causal, heritability)
    saveRDS(result, file.path(output_dir,
                              sprintf("transferability_m%d_h2%s.RDS",
                                      m_causal, format(heritability))))
    cat(sprintf("m = %4d, heritability = %.1f, done in %.0f seconds\n",
                m_causal, heritability,
                as.numeric(difftime(Sys.time(), started, units = "secs"))))
    flush.console()
  }
}
