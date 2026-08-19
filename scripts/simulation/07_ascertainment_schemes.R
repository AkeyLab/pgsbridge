# scripts/simulation/07_ascertainment_schemes.R
#
# Produces the data behind S8 Fig: intrinsic transferability under five schemes
# for deciding which causal variants count as common and which count as rare.
#
# Figure 3 classifies variants on the three populations pooled. That is one
# choice among several, and a variant that is rare in the pooled sample need not
# be rare in any single population. This script repeats the Figure 3 sweep under
# five rules. Writing m_S for the minor allele count in population S out of the
# 2,000 alleles it contributes, and m_global for the minor allele count out of
# all 6,000, the threshold of 0.01 corresponds to 20 within a population and 60
# globally, and the rules are
#
#   Model 1  all three   common: m_S >= 20 in all three populations
#                        rare:   m_S <  20 in all three populations
#                        A variant may be absent in one or two populations and
#                        still be rare. Variants common in one population and
#                        rare in another are in neither pool. Requiring a rare
#                        variant to segregate in all three as well leaves only
#                        1,089 variants, too few to draw from, so that stricter
#                        mirror is not used.
#   Model 2  European    common: 0 < n_EUR < 2000 and m_EUR >= 20
#                        rare:   0 < n_EUR < 2000 and m_EUR <  20
#                        Variants monomorphic in Europeans are in neither pool.
#   Model 3  African     the same rule applied to the African sample.
#   Model 4a pooled      common: m_global >= 60, rare: m_global < 60. This is
#                        the rule used for Figure 3, and it is also the global
#                        allele-frequency rule, since pooling the three samples
#                        and computing one frequency over 6,000 alleles is the
#                        same calculation.
#   Model 4b random      no ascertainment. Causal variants are drawn uniformly
#                        from the whole panel; the Model 4a split is applied
#                        afterwards only to decide which of them get the larger
#                        effect size.
#
# Everything else matches 04_traits_rare_and_common.R at the Figure 3 setting:
# 1,000 causal variants, total heritability 0.1 in the European sample, 200
# replicates, accuracy on all 1,000 individuals of each population.
#
# The pools are built from exact integer allele counts rather than from
# frequencies. Averaging the three per-population frequencies in floating point
# misclassifies the variants sitting exactly on the threshold, and how many
# depends on the order they are added in. The integer rule reproduces the
# `--maf 0.01` panel that coalescent/02_prepare_genotypes.sh writes, variant for variant.
#
# Runtime: several hours for the five schemes. The exploratory version of this
# analysis, which also varies the population the heritability is calibrated in
# and the level it is calibrated to, is kept outside this repository.
#
# Output:  data/ascertainment/transferability_by_scheme.csv
#          data/ascertainment/pool_sizes.csv
#
# Usage:  Rscript scripts/simulation/07_ascertainment_schemes.R

source("scripts/config.R")

suppressMessages({
  library(BEDMatrix)
  library(genio)
})

set.seed(20230210)

N_REPLICATES <- 200
INTERCEPT <- 0.1
NOISE_VARIANCE <- 0.1
HERITABILITY <- 0.1
CALIBRATION_POPULATION <- "EUR"

N_CAUSAL <- 1000L
RARE_FRACTIONS <- seq(0.1, 0.9, by = 0.1)
EFFECT_RATIOS <- c(1, 2, 5, 10)
RANDOM_N_CAUSAL <- c(1000L, 2000L, 5000L)

POPULATIONS <- c("AFR", "EUR", "ASN")
N_ALLELES <- 2000L                       # 1,000 diploid individuals
MAF_COUNT_WITHIN <- 20L                  # 0.01 of 2,000
MAF_COUNT_GLOBAL <- 60L                  # 0.01 of 6,000

# --- allele counts for every variant, in each population ---------------------
# Read from the PLINK panels in chunks rather than holding 3.2 million variants
# times three populations in memory. The panels are cut from one merged file, so
# a column index means the same variant in all three and the coded allele is the
# same; `load_variant_pools()` checks the first of those.

# Adding allele counts across populations is only meaningful if a variant's
# coded allele is the same in all three files. That is not automatic: PLINK will
# reset the coded allele to whichever is minor within each population unless the
# panels are cut with --keep-allele-order, and panels written without it disagree
# at tens of thousands of variants. Check before relying on it.
bim <- lapply(c(AFR = "afr", EUR = "eur", ASN = "asn"), function(stem) {
  read_bim(file.path(PRS_GENOTYPES, paste0(stem, "_all_variants.bim")))
})
for (population in names(bim)[-1]) {
  if (!identical(bim[[population]]$alt, bim[[1]]$alt) ||
      !identical(bim[[population]]$ref, bim[[1]]$ref)) {
    stop("the ", population, " panel does not share the coded allele of the ",
         names(bim)[1], " panel, so allele counts cannot be added across ",
         "populations. Rebuild the panels with --keep-allele-order, as ",
         "coalescent/02_prepare_genotypes.sh does.")
  }
}

pools_pooled <- load_variant_pools(PRS_GENOTYPES)
n_variants <- ncol(pools_pooled$genotypes[[1]])
n_individuals <- pools_pooled$n_individuals
cat(sprintf("%d individuals per population, %s variants\n",
            n_individuals, format(n_variants, big.mark = ",")))

allele_counts <- function(X, chunk_size = 20000L) {
  counts <- integer(ncol(X))
  for (start in seq(1L, ncol(X), by = chunk_size)) {
    stop_at <- min(start + chunk_size - 1L, ncol(X))
    block <- X[, start:stop_at, drop = FALSE]
    counts[start:stop_at] <- as.integer(colSums(block))
  }
  counts
}

counts <- vapply(POPULATIONS,
                 function(p) allele_counts(pools_pooled$genotypes[[p]]),
                 integer(n_variants))
colnames(counts) <- POPULATIONS
stopifnot(all(counts >= 0L), all(counts <= N_ALLELES))

minor <- pmin(counts, N_ALLELES - counts)
global <- rowSums(counts)
minor_global <- pmin(global, 3L * N_ALLELES - global)

# The panel holds segregating sites only, so nothing can be globally monomorphic.
if (any(global == 0L) || any(global == 3L * N_ALLELES)) {
  stop("a variant is monomorphic across all three populations, which means the ",
       "global allele count is being computed wrongly")
}

# --- the five schemes -------------------------------------------------------
index <- seq_len(n_variants)
segregating <- counts > 0L & counts < N_ALLELES

common_everywhere <- apply(minor >= MAF_COUNT_WITHIN, 1, all)
rare_everywhere   <- apply(minor <  MAF_COUNT_WITHIN, 1, all)

schemes <- list(
  model1_allthree = list(
    common = index[common_everywhere],
    rare   = index[rare_everywhere]),
  model2_european = list(
    common = index[segregating[, "EUR"] & minor[, "EUR"] >= MAF_COUNT_WITHIN],
    rare   = index[segregating[, "EUR"] & minor[, "EUR"] <  MAF_COUNT_WITHIN]),
  model3_african = list(
    common = index[segregating[, "AFR"] & minor[, "AFR"] >= MAF_COUNT_WITHIN],
    rare   = index[segregating[, "AFR"] & minor[, "AFR"] <  MAF_COUNT_WITHIN]),
  model4a_pooled = list(
    common = index[minor_global >= MAF_COUNT_GLOBAL],
    rare   = index[minor_global <  MAF_COUNT_GLOBAL]))
schemes$model4b_random <- list(
  all    = index,
  common = schemes$model4a_pooled$common,
  rare   = schemes$model4a_pooled$rare)

# Model 4a must reproduce the panel that 02_prepare_genotypes.sh wrote with
# --maf 0.01, variant for variant, since they are the same rule.
if (!identical(as.integer(sort(schemes$model4a_pooled$common)),
               as.integer(pools_pooled$common_index))) {
  stop("the pooled scheme does not reproduce the PLINK common panel")
}

pool_sizes <- do.call(rbind, lapply(names(schemes), function(scheme) {
  data.frame(model = scheme,
             common = length(schemes[[scheme]]$common),
             rare = length(schemes[[scheme]]$rare),
             stringsAsFactors = FALSE)
}))
print(pool_sizes)

# --- the sweep --------------------------------------------------------------
draw <- function(pool, m) {
  if (m == 0L) {
    empty <- lapply(pools_pooled$genotypes,
                    function(X) matrix(numeric(0), nrow = n_individuals))
    return(list(index = integer(0), genotypes = empty))
  }
  drawn <- pool[sample_causal_variants(length(pool), m)]
  list(index = drawn,
       genotypes = lapply(pools_pooled$genotypes,
                          function(X) X[, drawn, drop = FALSE]))
}

is_globally_rare <- logical(n_variants)
is_globally_rare[schemes$model4a_pooled$rare] <- TRUE

run_cell <- function(scheme, n_common, n_rare, effect_ratio, n_draw = 0L) {
  unascertained <- scheme == "model4b_random"
  ratios <- matrix(NA_real_, nrow = N_REPLICATES, ncol = length(ORDERED_PAIRS),
                   dimnames = list(NULL, PAIR_LABELS))
  rare_share <- numeric(N_REPLICATES)

  for (replicate in seq_len(N_REPLICATES)) {
    if (unascertained) {
      drawn <- draw(schemes[[scheme]]$all, n_draw)
      is_rare <- is_globally_rare[drawn$index]
      order_by <- c(which(!is_rare), which(is_rare))
      is_rare <- is_rare[order_by]
      genotypes <- lapply(drawn$genotypes,
                          function(X) X[, order_by, drop = FALSE])
    } else {
      common <- draw(schemes[[scheme]]$common, n_common)
      rare <- draw(schemes[[scheme]]$rare, n_rare)
      is_rare <- c(rep(FALSE, n_common), rep(TRUE, n_rare))
      genotypes <- lapply(POPULATIONS, function(p) {
        cbind(common$genotypes[[p]], rare$genotypes[[p]])
      })
      names(genotypes) <- POPULATIONS
    }
    rare_share[replicate] <- mean(is_rare)

    frequencies <- lapply(genotypes, allele_frequencies)
    calibration <- frequencies[[CALIBRATION_POPULATION]]
    effects <- split_effect_sizes(h2 = HERITABILITY, tau2 = NOISE_VARIANCE,
                                  p_common = calibration[!is_rare],
                                  p_rare = calibration[is_rare],
                                  s = effect_ratio)
    beta <- ifelse(is_rare, effects$rare, effects$common)

    traits <- lapply(POPULATIONS, function(p) {
      simulate_trait(genotypes[[p]], beta, alpha = INTERCEPT,
                     tau2 = NOISE_VARIANCE)
    })
    names(traits) <- POPULATIONS

    accuracy <- vapply(POPULATIONS, function(p) {
      empirical_accuracy(traits[[p]]$y, INTERCEPT + traits[[p]]$g)
    }, numeric(1))

    for (k in seq_along(ORDERED_PAIRS)) {
      ratios[replicate, k] <- accuracy[[ORDERED_PAIRS[[k]][2]]] /
        accuracy[[ORDERED_PAIRS[[k]][1]]]
    }
  }

  # The manuscript's estimator: trim each group once by the interquartile rule,
  # then average the per-replicate ratios.
  long <- data.frame(emp_tau = as.vector(ratios),
                     pop_comb = rep(PAIR_LABELS, each = N_REPLICATES),
                     stringsAsFactors = FALSE)
  summary <- group_summary(long, value = "emp_tau", groups = "pop_comb")
  summary$model <- scheme
  summary$ratio <- effect_ratio
  summary$n_causal <- if (unascertained) n_draw else n_common + n_rare
  summary$fraction <- if (unascertained) mean(rare_share) else n_rare / (n_common + n_rare)
  summary$rare_share <- mean(rare_share)
  summary
}

ORDERED_PAIRS <- list(c("EUR", "AFR"), c("EUR", "ASN"), c("AFR", "EUR"),
                      c("AFR", "ASN"), c("ASN", "EUR"), c("ASN", "AFR"))
PAIR_LABELS <- vapply(ORDERED_PAIRS, paste, character(1), collapse = "-")

rows <- list()
for (scheme in names(schemes)) {
  if (scheme == "model4b_random") {
    for (m in RANDOM_N_CAUSAL) for (s in EFFECT_RATIOS) {
      rows[[length(rows) + 1L]] <- run_cell(scheme, 0L, 0L, s, n_draw = m)
      cat(sprintf("%s m %5d ratio %2d done\n", scheme, m, s)); flush.console()
    }
  } else {
    for (fraction in RARE_FRACTIONS) {
      n_rare <- round(N_CAUSAL * fraction)
      for (s in EFFECT_RATIOS) {
        rows[[length(rows) + 1L]] <- run_cell(scheme, N_CAUSAL - n_rare,
                                              n_rare, s)
        cat(sprintf("%s rare %.1f ratio %2d done\n", scheme, fraction, s))
        flush.console()
      }
    }
  }
}

results <- do.call(rbind, rows)
output_dir <- file.path(PRS_DATA, "ascertainment")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(output_dir, "transferability_by_scheme.csv"),
          row.names = FALSE)
write.csv(pool_sizes, file.path(output_dir, "pool_sizes.csv"), row.names = FALSE)
cat("wrote", file.path(output_dir, "transferability_by_scheme.csv"), "\n")
