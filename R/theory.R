#' Additive genetic variance under Hardy-Weinberg equilibrium
#'
#' Evaluates `sigma^2 = 2 sum_i beta_i^2 p_i (1 - p_i)`, the additive genetic
#' variance a set of causal variants contributes when the population is in
#' Hardy-Weinberg equilibrium.  This is the quantity through which the
#' allele-frequency spectrum enters every other result in this package: the same
#' effect sizes acting on different allele frequencies give different additive
#' genetic variances, which is why two populations can differ in how predictable
#' a trait is even when their genetic architecture is identical.
#'
#' @param beta The length-`m` vector of effect sizes, one per causal variant.
#' @param p The length-`m` vector of reference-allele frequencies.
#'
#' @return The additive genetic variance, a scalar.
#'
#' @examples
#' # ten causal variants with equal effects
#' m_loci <- 10
#' beta <- rep(0.1, m_loci)
#' # allele frequencies drawn uniformly
#' p <- runif(m_loci, 0.05, 0.95)
#' additive_genetic_variance(beta, p)
#'
#' @export
additive_genetic_variance <- function(beta, p) {
    if (missing(beta))
        stop('`beta` is required!')
    if (missing(p))
        stop('`p` is required!')
    if (length(beta) != length(p))
        stop('`beta` and `p` must have the same length!')

    2 * sum(beta^2 * p * (1 - p))
}

#' Population mean of the trait
#'
#' Evaluates `mu = alpha + 2 sum_i beta_i p_i`, the expected trait value in a
#' population.  Simulations that measure prediction accuracy against a
#' constant-predictor baseline can use this exact mean rather than the sample
#' mean of one draw, which makes the estimated quantity match the closed form
#' rather than approximate it.
#'
#' @param alpha The intercept of the trait model.
#' @param beta The length-`m` vector of effect sizes.
#' @param p The length-`m` vector of reference-allele frequencies.
#'
#' @return The population mean, a scalar.
#'
#' @examples
#' beta <- rep(0.1, 10)
#' p <- runif(10, 0.05, 0.95)
#' population_mean(alpha = 0, beta = beta, p = p)
#'
#' @export
population_mean <- function(alpha, beta, p) {
    if (missing(alpha))
        stop('`alpha` is required!')
    if (missing(beta))
        stop('`beta` is required!')
    if (missing(p))
        stop('`p` is required!')
    if (length(beta) != length(p))
        stop('`beta` and `p` must have the same length!')

    alpha + 2 * sum(beta * p)
}

#' Best accuracy any polygenic score can reach in a population
#'
#' The optimal predictor of a trait under squared-error loss is its conditional
#' mean given the genotypes, `Yhat_j = alpha + sum_i beta_i X_ij`.  Because it is
#' built from the true causal variants and their true effects, it carries no
#' estimation error, so the accuracy it attains is an upper bound on the accuracy
#' of any score estimated from finite data.
#'
#' Its expected squared loss is the non-genetic variance `tau^2`, because the
#' non-genetic deviation averages to zero and cannot be predicted from genotype.
#' The baseline is the constant predictor, whose loss is the whole trait
#' variance.  Accuracy is one minus the ratio of the two losses,
#'
#' `PA = 2 sigma^2 Phibar / (2 sigma^2 Phibar + tau^2)`,
#'
#' the share of the trait's variance that is genetic.  Under Hardy-Weinberg
#' equilibrium the average self-kinship is `1/2` and this equals the
#' narrow-sense heritability.
#'
#' @param sigma2 The additive genetic variance, as returned by
#'   [additive_genetic_variance()].
#' @param phi_bar The average self-kinship.  Use `0.5` for a population in
#'   Hardy-Weinberg equilibrium with unrelated individuals.
#' @param tau2 The non-genetic variance.
#'
#' @return The best achievable prediction accuracy, in `[0, 1]`.
#'
#' @examples
#' # a population in Hardy-Weinberg equilibrium whose trait is half genetic
#' prediction_accuracy(sigma2 = 1, phi_bar = 0.5, tau2 = 1)
#'
#' # within-population structure raises the genetic variance and the ceiling
#' prediction_accuracy(sigma2 = 1, phi_bar = 0.6, tau2 = 1)
#'
#' @export
prediction_accuracy <- function(sigma2, phi_bar, tau2) {
    if (missing(sigma2))
        stop('`sigma2` is required!')
    if (missing(phi_bar))
        stop('`phi_bar` is required!')
    if (missing(tau2))
        stop('`tau2` is required!')

    genetic_variance <- 2 * sigma2 * phi_bar
    genetic_variance / (genetic_variance + tau2)
}

#' Intrinsic transferability between two populations
#'
#' Transferability is the ratio of the two populations' accuracy ceilings, that
#' is the fraction of the discovery population's ceiling retained in the target
#' population, `T = PA_B / PA_A`.  It equals one when the two ceilings coincide,
#' falls below one when accuracy is lost, and exceeds one when the target
#' population is the more predictable of the two.
#'
#' When both populations are in Hardy-Weinberg equilibrium the average
#' self-kinship is `1/2` in each and the ratio reduces to the ratio of
#' heritabilities, `h_B^2 / h_A^2`.
#'
#' @param sigma2_A,phi_bar_A,tau2_A The three determinants of the discovery
#'   population `A`.
#' @param sigma2_B,phi_bar_B,tau2_B The three determinants of the target
#'   population `B`.
#'
#' @return The transferability ratio.
#'
#' @examples
#' # two populations differing only in additive genetic variance
#' transferability(
#'     sigma2_A = 1.0, phi_bar_A = 0.5, tau2_A = 1,
#'     sigma2_B = 1.2, phi_bar_B = 0.5, tau2_B = 1
#' )
#'
#' # identical populations transfer perfectly
#' transferability(
#'     sigma2_A = 1, phi_bar_A = 0.5, tau2_A = 1,
#'     sigma2_B = 1, phi_bar_B = 0.5, tau2_B = 1
#' )
#'
#' @export
transferability <- function(sigma2_A, phi_bar_A, tau2_A,
                            sigma2_B, phi_bar_B, tau2_B) {
    for (argument in c('sigma2_A', 'phi_bar_A', 'tau2_A',
                       'sigma2_B', 'phi_bar_B', 'tau2_B'))
        if (eval(call('missing', as.name(argument))))
            stop('`', argument, '` is required!')

    prediction_accuracy(sigma2_B, phi_bar_B, tau2_B) /
        prediction_accuracy(sigma2_A, phi_bar_A, tau2_A)
}

#' Non-genetic variance achieving a target accuracy
#'
#' Inverts [prediction_accuracy()] for `tau^2`.  Use this when a simulation fixes
#' the heritability of one population and lets the others take whatever
#' heritability the same effect sizes produce in them.
#'
#' @param target_accuracy The accuracy, or equivalently the heritability under
#'   Hardy-Weinberg equilibrium, to achieve.  Must lie in `(0, 1)`.
#' @param sigma2 The additive genetic variance.
#' @param phi_bar The average self-kinship.  Defaults to the Hardy-Weinberg
#'   value `0.5`.
#'
#' @return The non-genetic variance giving that accuracy.
#'
#' @examples
#' # what non-genetic variance makes this trait 30 percent heritable?
#' tau2 <- noise_variance_for_accuracy(0.3, sigma2 = 1)
#' # confirm it round-trips
#' prediction_accuracy(sigma2 = 1, phi_bar = 0.5, tau2 = tau2)
#'
#' @export
noise_variance_for_accuracy <- function(target_accuracy, sigma2, phi_bar = 0.5) {
    if (missing(target_accuracy))
        stop('`target_accuracy` is required!')
    if (missing(sigma2))
        stop('`sigma2` is required!')
    if (any(target_accuracy <= 0) || any(target_accuracy >= 1))
        stop('`target_accuracy` must be strictly between 0 and 1!')

    2 * sigma2 * phi_bar * (1 - target_accuracy) / target_accuracy
}

#' Effect size shared by all causal variants achieving a target heritability
#'
#' With a single effect size `beta` shared by all `m` causal variants the
#' additive genetic variance is `2 beta^2 sum_i p_i (1 - p_i)`, so reaching
#' heritability `h2` at non-genetic variance `tau^2` in a Hardy-Weinberg
#' population requires the effect size below.
#'
#' @details
#' `beta = sqrt( h2 / (1 - h2) * tau^2 / (2 sum_i p_i (1 - p_i)) )`.
#'
#' @param h2 The target heritability, in `(0, 1)`.
#' @param tau2 The non-genetic variance.
#' @param p The reference-allele frequencies of the causal variants.
#'
#' @return The single effect size shared by every causal variant.
#'
#' @examples
#' p <- runif(100, 0.05, 0.95)
#' beta <- shared_effect_size(h2 = 0.5, tau2 = 1, p = p)
#' # the realised heritability matches the target
#' prediction_accuracy(
#'     sigma2 = additive_genetic_variance(rep(beta, length(p)), p),
#'     phi_bar = 0.5, tau2 = 1
#' )
#'
#' @export
shared_effect_size <- function(h2, tau2, p) {
    if (missing(h2))
        stop('`h2` is required!')
    if (missing(tau2))
        stop('`tau2` is required!')
    if (missing(p))
        stop('`p` is required!')
    if (any(h2 <= 0) || any(h2 >= 1))
        stop('`h2` must be strictly between 0 and 1!')

    sigma2_target <- h2 / (1 - h2) * tau2
    sqrt(sigma2_target / (2 * sum(p * (1 - p))))
}

#' Common and rare effect sizes for a given rare-to-common ratio
#'
#' Rare causal variants are given an effect `s` times that of common ones.
#' Writing `beta_c` for the common effect size, the additive genetic variance is
#' `2 beta_c^2 sum_c p_c (1 - p_c) + 2 s^2 beta_c^2 sum_r p_r (1 - p_r)`, so
#' matching a target heritability `h2` at non-genetic variance `tau^2` in a
#' Hardy-Weinberg population gives the common effect size below.
#'
#' @details
#' `beta_c = sqrt( h2/(1-h2) * tau^2 / ( 2 sum_c p_c(1-p_c) + 2 s^2 sum_r p_r(1-p_r) ) )`.
#'
#' @param h2 The target heritability, in `(0, 1)`.
#' @param tau2 The non-genetic variance.
#' @param p_common The reference-allele frequencies of the common causal
#'   variants.
#' @param p_rare The reference-allele frequencies of the rare causal variants.
#' @param s The rare-to-common effect-size ratio.
#'
#' @return A list with elements `common` and `rare`, the two effect sizes.
#'
#' @examples
#' p_common <- runif(500, 0.05, 0.95)
#' p_rare <- runif(500, 0.001, 0.01)
#' # rare variants act five times as strongly as common ones
#' split_effect_sizes(h2 = 0.5, tau2 = 1, p_common, p_rare, s = 5)
#'
#' @export
split_effect_sizes <- function(h2, tau2, p_common, p_rare, s) {
    if (missing(h2))
        stop('`h2` is required!')
    if (missing(tau2))
        stop('`tau2` is required!')
    if (missing(p_common))
        stop('`p_common` is required!')
    if (missing(p_rare))
        stop('`p_rare` is required!')
    if (missing(s))
        stop('`s` is required!')
    if (any(h2 <= 0) || any(h2 >= 1))
        stop('`h2` must be strictly between 0 and 1!')

    sigma2_target <- h2 / (1 - h2) * tau2
    denominator <- 2 * sum(p_common * (1 - p_common)) +
        s^2 * 2 * sum(p_rare * (1 - p_rare))
    beta_common <- sqrt(sigma2_target / denominator)

    list(common = beta_common, rare = s * beta_common)
}
