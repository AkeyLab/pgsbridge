# simulation/07_heritability_decomposition.R
#
# Splits each population's heritability into the part carried by common causal
# variants and the part carried by rare ones. Produces the data behind Figure 4.
#
# Why this figure exists
# ----------------------
# In the coalescent setting every population is in Hardy-Weinberg equilibrium,
# so transferability reduces exactly to the ratio of heritabilities. Figure 4
# therefore shows the quantity that drives Figures 2 and 3 directly, rather than
# through the ratio. It makes visible where the African sample's advantage comes
# from: as the rare fraction rises, rare variants account for a growing share of
# heritability in the African sample, and the gap between it and the two
# non-African samples widens.
#
# What is computed
# ----------------
# For each combination of the rare fraction and the effect-size ratio, traits
# are simulated exactly as in 06_traits_rare_and_common.R, and for each
# population the genetic value is split into its common-variant and
# rare-variant parts. Each part's variance is divided by the total trait
# variance, giving two heritability shares that sum to the total heritability up
# to the sampling covariance between the two parts.
#
# The shares are averaged over replicates and written as one tidy table per
# combination of heritability and effect-size ratio.
#
# Runtime: hours. Use 07_heritability_decomposition.slurm.
#
# Usage:  Rscript simulation/07_heritability_decomposition.R

source("config.R")
source(file.path(PRS_ROOT, "R", "theory.R"))
source(file.path(PRS_ROOT, "R", "traits.R"))
source(file.path(PRS_ROOT, "R", "variant_pools.R"))

set.seed(20230309)

CALIBRATION_POPULATION <- "EUR"
POPULATION_ORDER <- c("EUR", "AFR", "ASN")

N_REPLICATES <- 200
INTERCEPT <- 0.1
NOISE_VARIANCE <- 0.1

N_CAUSAL_VARIANTS <- 1000
RARE_FRACTIONS <- seq(0.1, 0.9, by = 0.1)
EFFECT_SIZE_RATIOS <- c(2, 5, 10)
HERITABILITIES <- c(0.1, 0.8)

pools <- load_variant_pools(PRS_GENOTYPES)
cat(sprintf("%d individuals, %d common and %d rare variants\n",
            pools$n_individuals, length(pools$common_index),
            length(pools$rare_index)))

#' Average heritability shares over replicates, for one parameter combination
#'
#' @return A data frame with columns `h2`, `pop`, `type` and `fraction`, holding
#'   one row per population per variant class.
decompose <- function(rare_fraction, effect_size_ratio, heritability) {
  n_rare <- round(N_CAUSAL_VARIANTS * rare_fraction)
  n_common <- N_CAUSAL_VARIANTS - n_rare

  totals <- matrix(0, nrow = length(POPULATION_ORDER), ncol = 2,
                   dimnames = list(POPULATION_ORDER, c("common", "rare")))

  for (replicate in seq_len(N_REPLICATES)) {
    common_genotypes <- draw_causal_genotypes(pools, "common", n_common)$genotypes
    rare_genotypes <- draw_causal_genotypes(pools, "rare", n_rare)$genotypes

    effects <- split_effect_sizes(
      h2 = heritability, tau2 = NOISE_VARIANCE,
      p_common = allele_frequencies(common_genotypes[[CALIBRATION_POPULATION]]),
      p_rare = allele_frequencies(rare_genotypes[[CALIBRATION_POPULATION]]),
      s = effect_size_ratio)

    beta <- c(rep(effects$common, n_common), rep(effects$rare, n_rare))

    for (population in POPULATION_ORDER) {
      X <- cbind(common_genotypes[[population]], rare_genotypes[[population]])
      trait <- simulate_trait(X, beta, alpha = INTERCEPT, tau2 = NOISE_VARIANCE)
      shares <- heritability_split(common_genotypes[[population]],
                                   rare_genotypes[[population]],
                                   effects$common, effects$rare, trait$y)
      totals[population, "common"] <- totals[population, "common"] + shares$common
      totals[population, "rare"] <- totals[population, "rare"] + shares$rare
    }
  }

  totals <- totals / N_REPLICATES

  data.frame(
    h2 = as.vector(t(totals)),
    pop = rep(POPULATION_ORDER, each = 2),
    type = rep(c("common", "rare"), times = length(POPULATION_ORDER)),
    fraction = rare_fraction,
    stringsAsFactors = FALSE)
}

output_dir <- file.path(PRS_DATA, "heritability_decomposition")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (heritability in HERITABILITIES) {
  h2_tag <- sub("0\\.", "0", format(heritability))
  for (effect_size_ratio in EFFECT_SIZE_RATIOS) {
    started <- Sys.time()
    table <- do.call(rbind, lapply(RARE_FRACTIONS, decompose,
                                   effect_size_ratio = effect_size_ratio,
                                   heritability = heritability))
    saveRDS(table, file.path(output_dir,
                             sprintf("h2decomp_h2%s_s%d.RDS",
                                     h2_tag, effect_size_ratio)))
    cat(sprintf("heritability %.1f, ratio %2d: %.0f seconds\n",
                heritability, effect_size_ratio,
                as.numeric(difftime(Sys.time(), started, units = "secs"))))
    flush.console()
  }
}
