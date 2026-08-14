test_that("group_summary computes the mean and interval it documents", {
    df <- data.frame(value = c(1, 2, 3, 4, 5), group = "a")
    out <- group_summary(df, value = "value", groups = "group")
    expect_equal(out$value, 3)
    expect_equal(out$n, 5)
    expect_equal(out$sd, sd(1:5))
    expect_equal(out$se, sd(1:5) / sqrt(5))
    expect_equal(out$ci, out$se * qt(0.975, 4))
})

test_that("group_summary trims by the interquartile rule", {
    # 100 is far outside [q1 - 1.5 iqr, q3 + 1.5 iqr] and must be dropped
    df <- data.frame(value = c(1, 2, 3, 4, 100), group = "a")
    out <- group_summary(df, value = "value", groups = "group")
    expect_equal(out$n, 4)
    expect_equal(out$value, 2.5)
})

test_that("group_summary trims within each group separately", {
    df <- data.frame(value = c(1, 2, 3, 4, 100, 10, 11, 12, 13, 900),
                     group = rep(c("a", "b"), each = 5))
    out <- group_summary(df, value = "value", groups = "group")
    expect_equal(nrow(out), 2)
    expect_equal(out$n, c(4, 4))
    expect_equal(out$value, c(2.5, 11.5))
})

test_that("label_pair_direction maps every ordered pair", {
    ordered <- c("ASN-AFR", "AFR-ASN", "EUR-AFR", "AFR-EUR", "ASN-EUR", "EUR-ASN")
    out <- label_pair_direction(ordered)
    expect_equal(as.character(out$pop_pair),
                 c("ASN-AFR", "ASN-AFR", "EUR-AFR", "EUR-AFR", "ASN-EUR", "ASN-EUR"))
    # each unordered pair appears once in each direction
    expect_equal(as.character(out$direction),
                 rep(c("Discovery-Target", "Target-Discovery"), 3))
    expect_error(label_pair_direction("XXX-YYY"))
})
