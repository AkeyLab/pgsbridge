# R/theory.R
#
# The closed-form results of the paper, written as R functions. Nothing here
# simulates anything; every function is an exact evaluation of a formula from
# the Materials and methods section.
#
# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------
# Individual j in population S has genotype X_ij at variant i, counted as the
# number of copies of the reference allele, and phenotype
#
#     Y_j = alpha + sum_i beta_i X_ij + epsilon_j,
#
# where beta is the same in every population, and epsilon_j has mean zero and
# variance tau_S^2 given the genotypes. Two population-genetic quantities set
# the genotype moments: the reference-allele frequency p_i^S, and the kinship
# coefficient Phi_jk^S, the probability that an allele drawn at random from
# individual j and one drawn at random from individual k are identical by
# descent. These give
#
#     E[X_ij]        = 2 p_i^S,
#     Var[X_ij]      = 4 p_i^S (1 - p_i^S) Phi_jj^S,
#     Cov[X_ij, X_ik] = 4 p_i^S (1 - p_i^S) Phi_jk^S.
#
# ---------------------------------------------------------------------------
# The three determinants
# ---------------------------------------------------------------------------
# Only three population-specific quantities enter the results:
#
#   sigma_S^2   the additive genetic variance under Hardy-Weinberg equilibrium,
#               sigma_S^2 = 2 sum_i beta_i^2 p_i^S (1 - p_i^S);
#   Phibar_S    the average self-kinship, Phibar_S = (1/n) sum_j Phi_jj^S,
#               which equals 1/2 when the population is in Hardy-Weinberg
#               equilibrium and the sampled individuals are unrelated;
#   tau_S^2     the non-genetic variance.
#
# The genetic variance of the trait is 2 sigma_S^2 Phibar_S, and the total
# trait variance is 2 sigma_S^2 Phibar_S + tau_S^2.

# ---------------------------------------------------------------------------

#' Additive genetic variance under Hardy-Weinberg equilibrium
#'
#' Evaluates sigma^2 = 2 sum_i beta_i^2 p_i (1 - p_i).
#'
#' @param beta Numeric vector of effect sizes, one per causal variant.
#' @param p Numeric vector of reference-allele frequencies, same length as
#'   `beta`.
#' @return A single number, the additive genetic variance sigma^2.
additive_genetic_variance <- function(beta, p) {
  stopifnot(length(beta) == length(p))
  2 * sum(beta^2 * p * (1 - p))
}

#' Population mean of the trait
#'
#' Evaluates mu = alpha + 2 sum_i beta_i p_i.
#'
#' @param alpha Intercept of the trait model.
#' @param beta Numeric vector of effect sizes.
#' @param p Numeric vector of reference-allele frequencies.
#' @return A single number, the population mean mu.
population_mean <- function(alpha, beta, p) {
  stopifnot(length(beta) == length(p))
  alpha + 2 * sum(beta * p)
}

#' Best achievable prediction accuracy in one population
#'
#' The optimal predictor of the trait under squared-error loss is the
#' conditional mean given the genotypes, Yhat_j = alpha + sum_i beta_i X_ij.
#' Because it is built from the true causal variants and their true effects, it
#' carries no estimation error, and the accuracy it attains is an upper bound on
#' the accuracy of any score estimated from finite data.
#'
#' Its expected squared loss is tau^2, the non-genetic variance, because the
#' non-genetic deviation averages to zero and cannot be predicted from genotype.
#' The baseline is the constant predictor mu, whose loss is the whole trait
#' variance. Accuracy is one minus the ratio of the two losses,
#'
#'     PA = 1 - tau^2 / (2 sigma^2 Phibar + tau^2)
#'        = 2 sigma^2 Phibar / (2 sigma^2 Phibar + tau^2),
#'
#' which is the share of the trait's variance that is genetic. Under
#' Hardy-Weinberg equilibrium Phibar = 1/2 and this equals the narrow-sense
#' heritability.
#'
#' @param sigma2 Additive genetic variance sigma^2.
#' @param phi_bar Average self-kinship Phibar.
#' @param tau2 Non-genetic variance tau^2.
#' @return The best achievable prediction accuracy, a number in [0, 1).
prediction_accuracy <- function(sigma2, phi_bar, tau2) {
  genetic_variance <- 2 * sigma2 * phi_bar
  genetic_variance / (genetic_variance + tau2)
}

#' Best-case transferability from a discovery to a target population
#'
#' Transferability is the ratio of the two populations' best achievable
#' accuracies, that is, the fraction of the discovery population's ceiling that
#' is retained in the target population:
#'
#'     T_{A->B} = PA_B / PA_A
#'              =   2 sigma_B^2 Phibar_B / (2 sigma_B^2 Phibar_B + tau_B^2)
#'                * (2 sigma_A^2 Phibar_A + tau_A^2) / (2 sigma_A^2 Phibar_A).
#'
#' It equals one when the two ceilings coincide, falls below one when accuracy
#' is lost, and exceeds one when the target population is the more predictable
#' of the two. When both populations are in Hardy-Weinberg equilibrium,
#' Phibar_A = Phibar_B = 1/2 and the ratio reduces to h_B^2 / h_A^2.
#'
#' @param sigma2_A,phi_bar_A,tau2_A Determinants of the discovery population A.
#' @param sigma2_B,phi_bar_B,tau2_B Determinants of the target population B.
#' @return The transferability ratio.
transferability <- function(sigma2_A, phi_bar_A, tau2_A,
                            sigma2_B, phi_bar_B, tau2_B) {
  prediction_accuracy(sigma2_B, phi_bar_B, tau2_B) /
    prediction_accuracy(sigma2_A, phi_bar_A, tau2_A)
}

#' Non-genetic variance that achieves a target heritability
#'
#' Inverts PA = 2 sigma^2 Phibar / (2 sigma^2 Phibar + tau^2) for tau^2. Used
#' when a simulation fixes the heritability of one population and lets the other
#' populations take whatever heritability the same effect sizes produce.
#'
#' @param target_accuracy The heritability (or accuracy) to achieve, in (0, 1).
#' @param sigma2 Additive genetic variance sigma^2.
#' @param phi_bar Average self-kinship Phibar. Defaults to the Hardy-Weinberg
#'   value 1/2.
#' @return The non-genetic variance tau^2 giving that accuracy.
noise_variance_for_accuracy <- function(target_accuracy, sigma2, phi_bar = 0.5) {
  stopifnot(target_accuracy > 0, target_accuracy < 1)
  2 * sigma2 * phi_bar * (1 - target_accuracy) / target_accuracy
}

#' Common effect size that achieves a target heritability
#'
#' With a single effect size beta shared by all m causal variants,
#' sigma^2 = 2 beta^2 sum_i p_i (1 - p_i). Setting the heritability of a
#' Hardy-Weinberg population to h2 with non-genetic variance tau^2 requires
#' 2 sigma^2 * (1/2) = sigma^2 = h2 / (1 - h2) * tau^2, so
#'
#'     beta = sqrt( h2 / (1 - h2) * tau^2 / (2 sum_i p_i (1 - p_i)) ).
#'
#' @param h2 Target heritability in (0, 1).
#' @param tau2 Non-genetic variance tau^2.
#' @param p Reference-allele frequencies of the causal variants.
#' @return The single effect size beta shared by all causal variants.
shared_effect_size <- function(h2, tau2, p) {
  sigma2_target <- h2 / (1 - h2) * tau2
  sqrt(sigma2_target / (2 * sum(p * (1 - p))))
}

#' Common and rare effect sizes for a given rare-to-common ratio
#'
#' Rare causal variants are assigned an effect size s times that of common
#' causal variants. Writing beta_c for the common effect size, the additive
#' genetic variance is
#'
#'     sigma^2 = 2 beta_c^2 sum_c p_c (1 - p_c) + 2 s^2 beta_c^2 sum_r p_r (1 - p_r),
#'
#' so matching a target heritability h2 at non-genetic variance tau^2 gives
#'
#'     beta_c = sqrt( h2/(1-h2) * tau^2
#'                    / ( 2 sum_c p_c(1-p_c) + 2 s^2 sum_r p_r(1-p_r) ) ).
#'
#' @param h2 Target heritability in (0, 1).
#' @param tau2 Non-genetic variance tau^2.
#' @param p_common Reference-allele frequencies of the common causal variants.
#' @param p_rare Reference-allele frequencies of the rare causal variants.
#' @param s Rare-to-common effect-size ratio.
#' @return A list with elements `common` and `rare`, the two effect sizes.
split_effect_sizes <- function(h2, tau2, p_common, p_rare, s) {
  sigma2_target <- h2 / (1 - h2) * tau2
  denominator <- 2 * sum(p_common * (1 - p_common)) +
    s^2 * 2 * sum(p_rare * (1 - p_rare))
  beta_common <- sqrt(sigma2_target / denominator)
  list(common = beta_common, rare = s * beta_common)
}
