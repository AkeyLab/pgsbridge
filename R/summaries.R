#' Group means and confidence intervals, after trimming outlying replicates
#'
#' Transferability is a ratio of two estimated accuracies. When the denominator
#' happens to be small in a replicate, the ratio can take an extreme value, and
#' a handful of such replicates would dominate an ordinary mean. Within each
#' group the replicates are therefore trimmed once by the interquartile rule:
#' with q1 and q3 the first and third quartiles and iqr = q3 - q1, a replicate is
#' kept if it lies in `[q1 - 1.5 iqr, q3 + 1.5 iqr]`. The mean, the standard
#' error, and the confidence interval are computed from the kept replicates.
#'
#' The confidence interval half-width is qt(1 - (1 - level)/2, N - 1) times the
#' standard error, where N is the number of kept replicates in the group.
#'
#' @param data A data frame of replicate-level results.
#' @param value Name of the column holding the value to summarise.
#' @param groups Character vector of column names defining the groups.
#' @param level Confidence level. Defaults to 0.95.
#' @return A data frame with one row per group, holding the grouping columns and
#'   the columns `n` (kept replicates), the mean under its original name, `sd`,
#'   `se`, and `ci` (the confidence interval half-width).
#' @export
group_summary <- function(data, value, groups, level = 0.95) {
  stopifnot(value %in% names(data), all(groups %in% names(data)))

  trimmed <- data %>%
    group_by(across(all_of(groups))) %>%
    filter({
      v <- .data[[value]]
      q <- quantile(v, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      iqr <- q[2] - q[1]
      !is.na(v) & v >= q[1] - 1.5 * iqr & v <= q[2] + 1.5 * iqr
    }) %>%
    ungroup()

  trimmed %>%
    group_by(across(all_of(groups))) %>%
    summarise(n = dplyr::n(),
              .mean = mean(.data[[value]]),
              sd = sd(.data[[value]]),
              .groups = "drop") %>%
    mutate(se = sd / sqrt(n),
           ci = se * qt(level / 2 + 0.5, pmax(n - 1, 1))) %>%
    rename(!!value := .mean) %>%
    as.data.frame()
}

#' Label an ordered population pair with its unordered pair and its direction
#'
#' The simulations record an ordered pair such as "ASN-AFR", meaning a discovery
#' population of East Asian and a target population of African. Figures 2, S2
#' and S3 draw the two directions of each unordered pair as the two ends of a
#' line, so each row needs both the unordered pair and which direction it is.
#'
#' The direction is named for the order in which the two populations appear in
#' the unordered pair label. The unordered labels are "ASN-AFR", "EUR-AFR" and
#' "ASN-EUR", so for example the ordered pair "ASN-AFR" is the
#' discovery-then-target direction, and "AFR-ASN" is the reverse.
#'
#' @param ordered_pair Character vector of ordered pair labels.
#' @return A data frame with columns `pop_pair` and `direction`.
#' @export
label_pair_direction <- function(ordered_pair) {
  unordered <- c("ASN-AFR" = "ASN-AFR", "AFR-ASN" = "ASN-AFR",
                 "EUR-AFR" = "EUR-AFR", "AFR-EUR" = "EUR-AFR",
                 "ASN-EUR" = "ASN-EUR", "EUR-ASN" = "ASN-EUR")
  forward <- c("ASN-AFR", "EUR-AFR", "ASN-EUR")
  ordered_pair <- as.character(ordered_pair)
  stopifnot(all(ordered_pair %in% names(unordered)))
  data.frame(
    pop_pair = factor(unname(unordered[ordered_pair])),
    direction = factor(ifelse(ordered_pair %in% forward,
                              "Discovery-Target", "Target-Discovery")),
    stringsAsFactors = FALSE
  )
}

#' Read the three heritability levels of a transferability experiment
#'
#' Figures 2, S2 and S3 all show the same quantity at European heritability
#' 0.1, 0.4 and 0.8. This reads the three files, stacks them, labels the pair
#' and direction, and returns the group summary the figure plots.
#'
#' @param paths A character vector of three file paths, for heritability 0.1,
#'   0.4 and 0.8 in that order.
#' @return A data frame ready for `plot_direction_pairs()`.
#' @export
read_transferability_by_heritability <- function(paths) {
  stopifnot(length(paths) == 3)
  levels_h2 <- c("0.1", "0.4", "0.8")

  stacked <- do.call(rbind, Map(function(path, h2) {
    df <- readRDS(path)
    df$h2 <- h2
    df
  }, paths, levels_h2))

  stacked <- cbind(stacked, label_pair_direction(stacked$pop_comb))
  stacked$h2 <- factor(stacked$h2, levels = levels_h2)

  group_summary(stacked, value = "emp_tau",
                groups = c("direction", "h2", "pop_pair"))
}

#' Read a rare-variant sweep for one ordered population pair
#'
#' Figures 3 and S4 to S6 vary the fraction of causal variants that are rare and
#' the rare-to-common effect-size ratio. Results are stored one file per
#' combination of the two, each file holding every ordered population pair.
#'
#' @param data_dir Directory holding the files.
#' @param h2_tag File-name tag for the European heritability: "h201", "h204" or
#'   "h208".
#' @param ordered_pair The ordered pair to extract, for example "AFR-EUR".
#' @param ratios Rare-to-common effect-size ratios to read.
#' @param fractions Rare-variant fractions to read.
#' @return A data frame with one row per (fraction, ratio), holding the mean
#'   transferability and its confidence interval half-width.
#' @export
read_rare_variant_sweep <- function(data_dir, h2_tag, ordered_pair,
                                    ratios = c(1, 2, 5, 10),
                                    fractions = seq(0.1, 0.9, by = 0.1)) {
  fraction_tags <- sprintf("0%d", round(fractions * 10))

  stacked <- do.call(rbind, lapply(ratios, function(s) {
    do.call(rbind, Map(function(fraction, tag) {
      path <- file.path(data_dir, sprintf("transdf_mloci1000_%s_rare%s_s%d.RDS",
                                          h2_tag, tag, s))
      df <- readRDS(path)
      df$fraction <- fraction
      df$ratio <- s
      df
    }, fractions, fraction_tags))
  }))

  summary <- group_summary(stacked, value = "emp_tau",
                           groups = c("pop_comb", "fraction", "ratio"))
  summary <- summary[summary$pop_comb == ordered_pair, ]
  summary$ratio <- factor(summary$ratio, levels = ratios)
  summary
}
