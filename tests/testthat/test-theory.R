test_that("additive_genetic_variance matches an independent calculation", {
    # worked by hand: 2 * (0.01*0.09 + 0.04*0.25 + 0.09*0.09) = 2 * 0.0190 = 0.0380
    expect_equal(additive_genetic_variance(c(0.1, 0.2, 0.3), c(0.1, 0.5, 0.9)),
                 0.038)
    # it is the variance of a sum of independent binomials, so it must equal
    # the sum of the per-variant variances computed separately
    beta <- c(0.3, 0.7); p <- c(0.2, 0.6)
    expect_equal(additive_genetic_variance(beta, p),
                 additive_genetic_variance(beta[1], p[1]) +
                     additive_genetic_variance(beta[2], p[2]))
    # a monomorphic variant contributes nothing
    expect_equal(additive_genetic_variance(1, 0), 0)
    expect_equal(additive_genetic_variance(1, 1), 0)
    # variance is maximised at a frequency of one half
    expect_gt(additive_genetic_variance(1, 0.5), additive_genetic_variance(1, 0.2))
})

test_that("additive_genetic_variance rejects mismatched inputs", {
    expect_error(additive_genetic_variance(c(1, 2), 0.5), "same length")
    expect_error(additive_genetic_variance(), "required")
})

test_that("population_mean matches the mean of simulated genotypes", {
    set.seed(11)
    beta <- c(0.1, 0.2, -0.15)
    p <- c(0.3, 0.7, 0.5)
    n <- 200000
    X <- sapply(p, function(pi) rbinom(n, 2, pi))
    # the analytic mean should match the mean of a large simulated sample
    expect_equal(population_mean(1, beta, p), mean(1 + as.vector(X %*% beta)),
                 tolerance = 1e-3)
    # worked by hand: 1 + 2*(0.03 + 0.14 - 0.075) = 1.19
    expect_equal(population_mean(1, beta, p), 1.19)
})

test_that("prediction_accuracy is the genetic share of the trait variance", {
    # equal genetic and non-genetic variance gives a ceiling of one half
    expect_equal(prediction_accuracy(sigma2 = 1, phi_bar = 0.5, tau2 = 1), 0.5)
    # no non-genetic variance gives a ceiling of one
    expect_equal(prediction_accuracy(sigma2 = 1, phi_bar = 0.5, tau2 = 0), 1)
    # no genetic variance gives a ceiling of zero
    expect_equal(prediction_accuracy(sigma2 = 0, phi_bar = 0.5, tau2 = 1), 0)
    # the ceiling always lies in [0, 1)
    accuracy <- prediction_accuracy(runif(20, 0.1, 2), runif(20, 0.5, 1),
                                    runif(20, 0.1, 3))
    expect_true(all(accuracy >= 0 & accuracy < 1))
})

test_that("prediction_accuracy rises with genetic variance and falls with noise", {
    expect_gt(prediction_accuracy(2, 0.5, 1), prediction_accuracy(1, 0.5, 1))
    expect_gt(prediction_accuracy(1, 0.6, 1), prediction_accuracy(1, 0.5, 1))
    expect_lt(prediction_accuracy(1, 0.5, 2), prediction_accuracy(1, 0.5, 1))
})

test_that("transferability is one between identical populations", {
    expect_equal(transferability(1, 0.5, 1, 1, 0.5, 1), 1)
})

test_that("transferability is the reciprocal in the opposite direction", {
    forward <- transferability(1.0, 0.5, 1, 1.5, 0.6, 2)
    backward <- transferability(1.5, 0.6, 2, 1.0, 0.5, 1)
    expect_equal(forward * backward, 1)
})

test_that("transferability reduces to the heritability ratio under Hardy-Weinberg", {
    h2_A <- prediction_accuracy(1.0, 0.5, 1)
    h2_B <- prediction_accuracy(1.5, 0.5, 1)
    expect_equal(transferability(1.0, 0.5, 1, 1.5, 0.5, 1), h2_B / h2_A)
})

test_that("noise_variance_for_accuracy inverts prediction_accuracy", {
    for (target in c(0.1, 0.4, 0.8)) {
        tau2 <- noise_variance_for_accuracy(target, sigma2 = 1.3, phi_bar = 0.5)
        expect_equal(prediction_accuracy(1.3, 0.5, tau2), target)
    }
    expect_error(noise_variance_for_accuracy(0, 1), "between 0 and 1")
    expect_error(noise_variance_for_accuracy(1, 1), "between 0 and 1")
})

test_that("shared_effect_size achieves the heritability it is asked for", {
    set.seed(1)
    p <- runif(200, 0.05, 0.95)
    for (h2 in c(0.1, 0.4, 0.8)) {
        beta <- shared_effect_size(h2, tau2 = 1, p = p)
        sigma2 <- additive_genetic_variance(rep(beta, length(p)), p)
        expect_equal(prediction_accuracy(sigma2, 0.5, 1), h2)
    }
    expect_error(shared_effect_size(1, 1, 0.5), "between 0 and 1")
})

test_that("split_effect_sizes achieves the heritability and the ratio", {
    set.seed(2)
    p_common <- runif(300, 0.05, 0.95)
    p_rare <- runif(300, 0.001, 0.01)
    for (s in c(1, 2, 5, 10)) {
        effects <- split_effect_sizes(h2 = 0.4, tau2 = 1, p_common, p_rare, s)
        # the requested ratio between the two effect sizes
        expect_equal(effects$rare / effects$common, s)
        # and the requested heritability overall
        sigma2 <- additive_genetic_variance(
            c(rep(effects$common, length(p_common)), rep(effects$rare, length(p_rare))),
            c(p_common, p_rare))
        expect_equal(prediction_accuracy(sigma2, 0.5, 1), 0.4)
    }
})

test_that("split_effect_sizes with a ratio of one matches shared_effect_size", {
    set.seed(3)
    p_common <- runif(100, 0.05, 0.95)
    p_rare <- runif(100, 0.001, 0.01)
    effects <- split_effect_sizes(0.5, 1, p_common, p_rare, s = 1)
    expect_equal(effects$common,
                 shared_effect_size(0.5, 1, c(p_common, p_rare)))
})
