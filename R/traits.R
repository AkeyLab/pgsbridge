#' Draw causal variants spread evenly across a variant pool
#'
#' The pool of `n_pool` variants is cut into `m_causal` consecutive blocks of
#' width `floor(n_pool / m_causal)`, and one variant is drawn uniformly at
#' random from each block. Causal variants are therefore separated by
#' `floor(n_pool / m_causal)` variants on average. For the pools used here that
#' runs from a few hundred variants, for 1,000 causal variants drawn from the
#' common pool, to tens of thousands, for 20 causal variants drawn from the same
#' pool, so that any two of them sit far apart in the genome and carry
#' negligible linkage disequilibrium. Individual separations vary:
#' two consecutive draws can land at the end of one block and the start of the
#' next, so the spacing is a mean rather than a guaranteed minimum.
#'
#' @param n_pool Number of variants in the pool.
#' @param m_causal Number of causal variants to draw.
#' @return An integer vector of `m_causal` column indices into the pool.
#' @export
sample_causal_variants <- function(n_pool, m_causal) {
  block_width <- floor(n_pool / m_causal)
  stopifnot(block_width >= 1)
  as.integer(sample.int(block_width, m_causal, replace = TRUE) +
               (seq_len(m_causal) - 1L) * block_width)
}

#' Reference-allele frequencies of a genotype matrix
#'
#' @param X Genotype matrix, individuals in rows and variants in columns, coded
#'   as the number of copies of the reference allele.
#' @return A numeric vector of allele frequencies, one per column.
#' @export
allele_frequencies <- function(X) {
  colSums(X) / (2 * nrow(X))
}

#' Simulate a trait from genotypes and effect sizes
#'
#' Draws Y_j = alpha + sum_i beta_i X_ij + epsilon_j with
#' epsilon_j ~ Normal(0, tau2), independently across individuals.
#'
#' @param X Genotype matrix, individuals in rows and causal variants in columns.
#' @param beta Numeric vector of effect sizes, one per column of `X`.
#' @param alpha Intercept of the trait model.
#' @param tau2 Non-genetic variance.
#' @return A list with `y`, the simulated trait, and `g`, the genetic value
#'   sum_i beta_i X_ij, which is also the optimal predictor up to the intercept.
#' @export
simulate_trait <- function(X, beta, alpha, tau2) {
  g <- as.vector(X %*% beta)
  list(y = alpha + g + rnorm(nrow(X), mean = 0, sd = sqrt(tau2)),
       g = g)
}

#' Empirical accuracy of the optimal predictor in one population
#'
#' The optimal predictor is Yhat_j = alpha + g_j, built from the true causal
#' variants and their true effects, so the trait decomposes as
#' Y_j = Yhat_j + epsilon_j with epsilon independent of Yhat. Accuracy is
#' estimated by the squared Pearson correlation between the trait and the
#' predictor. Under that decomposition,
#'
#'     cor(Y, Yhat)^2 = Cov(Y, Yhat)^2 / (Var(Y) Var(Yhat))
#'                    = Var(Yhat)^2 / ((Var(Yhat) + tau^2) Var(Yhat))
#'                    = Var(Yhat) / (Var(Yhat) + tau^2),
#'
#' and since Var(Yhat) = 2 sigma^2 Phibar this is exactly the population
#' quantity that `prediction_accuracy()` evaluates in closed form. The
#' equivalent form 1 - mean((Y - Yhat)^2) / Var(Y) estimates the same quantity,
#' because mean((Y - Yhat)^2) estimates tau^2 and Var(Y) estimates
#' 2 sigma^2 Phibar + tau^2. The two agree in the limit and differ only by the
#' finite-sample behaviour of the estimators they are built from.
#'
#' The individuals within each simulated population are unrelated, so the
#' average off-diagonal kinship is zero and the sample-mean baseline implicit in
#' both forms is unbiased for the closed form.
#'
#' @param y Simulated trait values.
#' @param y_hat Optimal-predictor values for the same individuals.
#' @return The estimated prediction accuracy.
#' @export
empirical_accuracy <- function(y, y_hat) {
  as.numeric(cor(y, y_hat))^2
}

#' Empirical transferability between two populations
#'
#' The ratio of the optimal predictor's accuracy measured in the target
#' population to its accuracy measured in the discovery population. The same
#' effect sizes are used in both, so the only thing that differs is the
#' populations' genotypes and hence their allele-frequency spectra.
#'
#' @param y_discovery,y_hat_discovery Trait and optimal predictor in the
#'   discovery population.
#' @param y_target,y_hat_target Trait and optimal predictor in the target
#'   population.
#' @return The estimated transferability ratio.
#' @export
empirical_transferability <- function(y_discovery, y_hat_discovery,
                                      y_target, y_hat_target) {
  empirical_accuracy(y_target, y_hat_target) /
    empirical_accuracy(y_discovery, y_hat_discovery)
}

#' Heritability explained by two disjoint groups of causal variants
#'
#' Splits the genetic value into the part carried by common causal variants and
#' the part carried by rare causal variants, and expresses each as a share of
#' the total trait variance. Used for the heritability decomposition.
#'
#' @param X_common,X_rare Genotype matrices for the common and the rare causal
#'   variants, on the same individuals.
#' @param beta_common,beta_rare The two effect sizes.
#' @param y The simulated trait for those individuals.
#' @return A list with `common` and `rare`, the two heritability shares.
#' @export
heritability_split <- function(X_common, X_rare, beta_common, beta_rare, y) {
  g_common <- as.vector(X_common %*% rep(beta_common, ncol(X_common)))
  g_rare <- as.vector(X_rare %*% rep(beta_rare, ncol(X_rare)))
  list(common = var(g_common) / var(y),
       rare = var(g_rare) / var(y))
}
