test_that("sample_causal_variants draws the right number, in range", {
    set.seed(1)
    index <- sample_causal_variants(n_pool = 10000, m_causal = 100)
    expect_length(index, 100)
    expect_true(all(index >= 1 & index <= 10000))
    expect_false(anyDuplicated(index) > 0)
})

test_that("sample_causal_variants spreads variants across the pool", {
    set.seed(2)
    index <- sample_causal_variants(n_pool = 10000, m_causal = 100)
    # one variant per consecutive block of floor(10000/100) = 100
    expect_true(all(diff(sort(index)) >= 1))
    expect_equal(mean(diff(sort(index))), 100, tolerance = 0.25)
    # and they span the pool rather than clustering
    expect_gt(max(index) - min(index), 9000)
})

test_that("allele_frequencies counts the reference allele", {
    # every individual homozygous for the reference allele
    expect_equal(allele_frequencies(matrix(2, nrow = 10, ncol = 3)), rep(1, 3))
    # every individual homozygous for the other allele
    expect_equal(allele_frequencies(matrix(0, nrow = 10, ncol = 3)), rep(0, 3))
    # all heterozygous
    expect_equal(allele_frequencies(matrix(1, nrow = 10, ncol = 3)), rep(0.5, 3))
})

test_that("simulate_trait returns the genetic value and a noisy trait", {
    set.seed(3)
    X <- matrix(rbinom(2000 * 5, 2, 0.4), 2000, 5)
    beta <- rep(0.2, 5)
    trait <- simulate_trait(X, beta, alpha = 1.5, tau2 = 0.25)
    # the genetic value is exactly X %*% beta, with no intercept and no noise
    expect_equal(trait$g, as.vector(X %*% beta))
    # the residual has the requested variance and mean zero
    residual <- trait$y - 1.5 - trait$g
    expect_equal(mean(residual), 0, tolerance = 0.05)
    expect_equal(var(residual), 0.25, tolerance = 0.05)
})

test_that("empirical_accuracy recovers the heritability it was given", {
    set.seed(4)
    m_loci <- 500
    p <- runif(m_loci, 0.1, 0.9)
    X <- matrix(rbinom(20000 * m_loci, 2, rep(p, each = 20000)), 20000, m_loci)
    for (h2 in c(0.2, 0.6)) {
        beta <- rep(shared_effect_size(h2, tau2 = 1, p = allele_frequencies(X)),
                    m_loci)
        trait <- simulate_trait(X, beta, alpha = 0, tau2 = 1)
        expect_equal(empirical_accuracy(trait$y, trait$g), h2, tolerance = 0.03)
    }
})

test_that("empirical_transferability is the ratio of the two accuracies", {
    set.seed(5)
    y_a <- rnorm(500); g_a <- y_a + rnorm(500, sd = 0.5)
    y_b <- rnorm(500); g_b <- y_b + rnorm(500, sd = 0.8)
    expect_equal(empirical_transferability(y_a, g_a, y_b, g_b),
                 empirical_accuracy(y_b, g_b) / empirical_accuracy(y_a, g_a))
})
