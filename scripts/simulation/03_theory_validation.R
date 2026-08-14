# simulation/03_theory_validation.R
#
# Validates the closed form of transferability against forward simulation.
# Produces the data behind Figure 1a.
#
# What is being tested
# --------------------
# The closed form gives transferability as a ratio of two accuracy ceilings
# evaluated at the populations' determinants. This script asks whether a literal
# simulation of the best case reproduces that ratio. For each of many simulated
# pairs of populations it
#
#   1. draws the two populations' allele frequencies, within-population
#      structure and non-genetic variance independently, so that all three
#      determinants differ between them;
#   2. evaluates the closed form at the realised determinants, giving
#      T_theoretical;
#   3. simulates genotypes and traits in each population, measures the optimal
#      predictor's accuracy in each, and takes their ratio, giving T_empirical.
#
# Because the two populations differ in every determinant, agreement between the
# two shows that the full closed form governs transferability, and not merely
# the ratio of heritabilities that it reduces to under Hardy-Weinberg
# equilibrium.
#
# The simulation engine
# ---------------------
# Genotypes come from the bnpsd package, using a one-dimensional linear
# admixture over two subpopulations, which gives individuals within a population
# a non-trivial kinship structure and therefore an average self-kinship above
# the Hardy-Weinberg value of one half. Effect sizes are drawn once and shared
# between the two populations. The optimal predictor is the genetic value
# itself, so it carries no estimation error.
#
# The accuracy estimator
# ----------------------
# The constant-predictor baseline is taken about the population's true mean
# mu_S = 2 sum_i beta_i p_i^S rather than about the sample mean of one draw.
# Because mu_S is known exactly here, this makes the quantity the simulation
# estimates match the closed form rather than approximate it. Individuals in
# these populations are related, so a sample-mean baseline would carry a term in
# the average off-diagonal kinship that the closed form does not have.
#
# Accuracy is then a ratio of two averages, each estimated with error, and a
# ratio of two noisy averages is biased. To remove that bias the population is
# redrawn R times and the numerator and denominator losses are averaged across
# redraws before dividing:
#
#     PAhat_S = 1 -  ( (1/R) sum_r (1/n) sum_j (Y_j^{S,r} - g_j^{S,r})^2 )
#                  / ( (1/R) sum_r (1/n) sum_j (Y_j^{S,r} - mu_S)^2 ).
#
# Parameter ranges
# ----------------
# The differentiation coefficient of each population is drawn from [0.01, 0.18],
# the upper end being the widest fixation index observed among human continental
# groups, and the within-population inbreeding from [0.01, 0.10], matching the
# nearly unrelated samples that genome-wide association studies retain after
# removing relatives. Scenarios are then kept only if their realised
# heritabilities satisfy five conditions: h2_A >= 0.30, h2_B > 0.30,
# |h2_A - h2_B| <= 0.30, and h2_A and h2_B both at most 0.90. The first three
# keep each heritability inside the range reported for complex traits and the
# two populations close to one another; the two upper caps keep the realised
# heritabilities inside the range the parameters were drawn from. A sixth
# condition discards any scenario whose empirical transferability exceeds
# MAX_ABS_TRANSFERABILITY in absolute value, as a numerical failure rather than
# a result.
#
# Runtime: about two minutes on 28 cores. Use 03_theory_validation.slurm.
#
# Usage:  Rscript scripts/simulation/03_theory_validation.R

source("scripts/config.R")

suppressMessages({
  library(bnpsd)
  library(parallel)
})

# L'Ecuyer's generator gives independent streams across the parallel workers.
RNGkind("L'Ecuyer-CMRG")
set.seed(20260630)

N_INDIVIDUALS <- 600L      # individuals per population
N_VARIANTS <- 400L         # causal variants
N_REDRAWS <- 50L           # redraws averaged over in the accuracy estimator
N_SUBPOPULATIONS <- 2L     # subpopulations in the admixture model
N_SCENARIOS <- 400L        # simulated pairs of populations

# Ranges of the drawn parameters.
DIFFERENTIATION_RANGE <- c(0.01, 0.18)   # per-population fixation index
INBREEDING_RANGE <- c(0.01, 0.10)        # within-population structure
ACCURACY_RANGE <- c(0.30, 0.90)          # target accuracy before realisation
ACCURACY_GAP <- 0.30                     # bound on |h2_A - h2_B|

# Allele frequencies are rejection-sampled into this interval, so that no causal
# variant is nearly monomorphic in either population.
FREQUENCY_RANGE <- c(0.1, 0.9)
MAX_REJECTION_ATTEMPTS <- 200L

# Outlier guard on the simulated ratio. A scenario whose empirical
# transferability falls outside this range is discarded as a numerical failure
# rather than a result.
MAX_ABS_TRANSFERABILITY <- 5

n_cores <- PRS_CORES
seeds <- 20000L + seq_len(N_SCENARIOS)

# The admixture proportions are fixed across scenarios; what varies between two
# populations is the inbreeding of their subpopulations.
admixture <- admix_prop_1d_linear(N_INDIVIDUALS, k_subpops = N_SUBPOPULATIONS,
                                  sigma = 1)

#' Draw one allele frequency from the Balding-Nichols distribution
#'
#' The Balding-Nichols distribution is a Beta distribution with mean equal to
#' the ancestral frequency and dispersion set by the differentiation
#' coefficient: with ancestral frequency a and coefficient f, the shape
#' parameters are (1-f)/f * a and (1-f)/f * (1-a). A larger f makes two
#' independent draws about the same ancestral frequency more different, which is
#' how the two populations come to have different additive genetic variances.
#'
#' Draws are rejected until they fall inside `FREQUENCY_RANGE`. If no draw
#' succeeds within `MAX_REJECTION_ATTEMPTS` the last draw is clamped, which
#' keeps the simulation running rather than failing on a rare extreme.
draw_balding_nichols <- function(ancestral_frequency, differentiation) {
  scale <- (1 - differentiation) / differentiation
  shape1 <- scale * ancestral_frequency
  shape2 <- scale * (1 - ancestral_frequency)
  for (attempt in seq_len(MAX_REJECTION_ATTEMPTS)) {
    value <- rbeta(1, shape1, shape2)
    if (value >= FREQUENCY_RANGE[1] && value <= FREQUENCY_RANGE[2]) return(value)
  }
  min(max(value, FREQUENCY_RANGE[1]), FREQUENCY_RANGE[2])
}

#' Average self-kinship implied by a set of subpopulation inbreeding values
mean_self_kinship <- function(inbreeding) {
  mean(diag(coanc_to_kinship(coanc_admix(admixture, inbreeding))))
}

#' Draw all parameters of one scenario
#'
#' The two populations share their ancestral allele frequencies and their effect
#' sizes, and differ in everything else. `target_accuracy` is the accuracy each
#' population would have if its genetic variance equalled the reference value.
#' Constraining the second draw around the first keeps the two targets close,
#' but it does not bound the realised gap, because each population's own genetic
#' variance shifts its realised heritability away from its target. The realised
#' gap is enforced afterwards by the |h2_A - h2_B| <= ACCURACY_GAP filter.
draw_scenario <- function(seed) {
  set.seed(seed)
  ancestral <- runif(N_VARIANTS, 0.01, 0.99)
  beta <- rnorm(N_VARIANTS, 0, sqrt(1 / N_VARIANTS))

  differentiation_A <- runif(1, DIFFERENTIATION_RANGE[1], DIFFERENTIATION_RANGE[2])
  differentiation_B <- runif(1, DIFFERENTIATION_RANGE[1], DIFFERENTIATION_RANGE[2])
  inbreeding_A <- runif(N_SUBPOPULATIONS, INBREEDING_RANGE[1], INBREEDING_RANGE[2])
  inbreeding_B <- runif(N_SUBPOPULATIONS, INBREEDING_RANGE[1], INBREEDING_RANGE[2])

  target_A <- runif(1, ACCURACY_RANGE[1], ACCURACY_RANGE[2])
  target_B <- runif(1,
                    max(ACCURACY_RANGE[1], target_A - ACCURACY_GAP),
                    min(ACCURACY_RANGE[2], target_A + ACCURACY_GAP))

  frequencies_A <- vapply(ancestral, draw_balding_nichols, numeric(1),
                          differentiation = differentiation_A)
  frequencies_B <- vapply(ancestral, draw_balding_nichols, numeric(1),
                          differentiation = differentiation_B)

  phi_A <- mean_self_kinship(inbreeding_A)
  phi_B <- mean_self_kinship(inbreeding_B)

  sigma2_A <- additive_genetic_variance(beta, frequencies_A)
  sigma2_B <- additive_genetic_variance(beta, frequencies_B)

  list(beta = beta,
       frequencies_A = frequencies_A, frequencies_B = frequencies_B,
       inbreeding_A = inbreeding_A, inbreeding_B = inbreeding_B,
       target_A = target_A, target_B = target_B,
       mean_A = population_mean(0, beta, frequencies_A),
       mean_B = population_mean(0, beta, frequencies_B),
       genetic_variance_A = 2 * sigma2_A * phi_A,
       genetic_variance_B = 2 * sigma2_B * phi_B)
}

#' Measure the optimal predictor's accuracy in one simulated population
#'
#' Redraws the population `N_REDRAWS` times, accumulating the numerator and
#' denominator losses, and divides only at the end.
measure_accuracy <- function(frequencies, beta, inbreeding, tau2, true_mean) {
  numerator <- 0
  denominator <- 0
  for (redraw in seq_len(N_REDRAWS)) {
    subpopulation_frequencies <- draw_p_subpops(frequencies, inbreeding)
    individual_frequencies <- make_p_ind_admix(subpopulation_frequencies,
                                               admixture)
    genotypes <- draw_genotypes_admix(individual_frequencies)
    genetic_value <- as.vector(beta %*% genotypes)
    trait <- genetic_value + rnorm(N_INDIVIDUALS, 0, sqrt(tau2))
    numerator <- numerator + mean((trait - genetic_value)^2)
    denominator <- denominator + mean((trait - true_mean)^2)
  }
  1 - (numerator / N_REDRAWS) / (denominator / N_REDRAWS)
}

started <- Sys.time()

# The non-genetic variance is put on an absolute scale, tau2 = (1-q)/q * Gbar
# with Gbar a fixed reference genetic variance shared by every scenario. Were
# tau2 instead set relative to each population's own genetic variance, the
# realised heritability would be q by construction and the additive genetic
# variance and self-kinship would cancel out of the closed form, so the test
# would no longer exercise them.
reference_genetic_variance <- mean(unlist(mclapply(seeds, function(seed) {
  scenario <- draw_scenario(seed)
  c(scenario$genetic_variance_A, scenario$genetic_variance_B)
}, mc.cores = n_cores)))

cat(sprintf("reference genetic variance = %.4f\n", reference_genetic_variance))
flush.console()

run_scenario <- function(seed) {
  scenario <- draw_scenario(seed)

  tau2_A <- (1 - scenario$target_A) / scenario$target_A * reference_genetic_variance
  tau2_B <- (1 - scenario$target_B) / scenario$target_B * reference_genetic_variance

  # Realised accuracies from the closed form. Because the genetic variance is
  # already 2 sigma^2 Phibar, this is the closed-form accuracy at the realised
  # determinants.
  accuracy_A <- scenario$genetic_variance_A /
    (scenario$genetic_variance_A + tau2_A)
  accuracy_B <- scenario$genetic_variance_B /
    (scenario$genetic_variance_B + tau2_B)

  measured_A <- measure_accuracy(scenario$frequencies_A, scenario$beta,
                                 scenario$inbreeding_A, tau2_A, scenario$mean_A)
  measured_B <- measure_accuracy(scenario$frequencies_B, scenario$beta,
                                 scenario$inbreeding_B, tau2_B, scenario$mean_B)

  c(T_theoretical = accuracy_B / accuracy_A,
    T_empirical = measured_B / measured_A,
    h2_A = accuracy_A,
    h2_B = accuracy_B)
}

results <- mclapply(seeds, run_scenario, mc.cores = n_cores, mc.preschedule = TRUE)

complete <- vapply(results, function(z) {
  is.numeric(z) && length(z) == 4L && all(is.finite(z))
}, logical(1))
all_scenarios <- as.data.frame(do.call(rbind, results[complete]))

# The five heritability constraints, plus the outlier guard.
keep <- all_scenarios$h2_A >= 0.30 &
  all_scenarios$h2_B > 0.30 &
  abs(all_scenarios$h2_A - all_scenarios$h2_B) <= ACCURACY_GAP &
  all_scenarios$h2_A <= ACCURACY_RANGE[2] &
  all_scenarios$h2_B <= ACCURACY_RANGE[2] &
  is.finite(all_scenarios$T_empirical) &
  abs(all_scenarios$T_empirical) <= MAX_ABS_TRANSFERABILITY
scenarios <- all_scenarios[keep, ]

r_squared <- cor(scenarios$T_theoretical, scenarios$T_empirical)^2

cat(sprintf(paste("scenarios simulated = %d, kept = %d, r^2 = %.4f,",
                  "theoretical transferability in [%.2f, %.2f]\n"),
            nrow(all_scenarios), nrow(scenarios), r_squared,
            min(scenarios$T_theoretical), max(scenarios$T_theoretical)))

output_dir <- file.path(PRS_DATA, "theory_validation")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(scenarios = scenarios,
             r_squared = r_squared,
             reference_genetic_variance = reference_genetic_variance,
             n_scenarios_simulated = nrow(all_scenarios),
             settings = list(n = N_INDIVIDUALS, m = N_VARIANTS,
                             redraws = N_REDRAWS,
                             subpopulations = N_SUBPOPULATIONS,
                             scenarios = N_SCENARIOS,
                             differentiation_range = DIFFERENTIATION_RANGE,
                             inbreeding_range = INBREEDING_RANGE,
                             accuracy_gap = ACCURACY_GAP)),
        file.path(output_dir, "distinct_populations.rds"))

cat(sprintf("written in %.0f seconds\n",
            as.numeric(difftime(Sys.time(), started, units = "secs"))))
