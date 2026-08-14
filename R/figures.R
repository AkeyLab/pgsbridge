#' Figure 1: theoretical determinants of intrinsic transferability
#'
#' Panel **a** plots empirical against theoretical transferability across
#' simulated pairs of populations that differ in all three determinants, with
#' the Pearson `r^2` annotated.  Panels **b** to **d** trace the closed-form
#' response of transferability to the additive genetic variance, the average
#' self-kinship and the non-genetic variance in turn.  Those three panels carry
#' no simulation: each curve is [transferability()] evaluated on a grid, so they
#' show the exact theoretical response with no Monte Carlo error.
#'
#' In each of panels **b** to **d** the dashed curve varies the determinant in
#' the discovery population and the solid curve varies it in the target
#' population.  A change in the discovery population enters the ratio through
#' its denominator and a change in the target through its numerator, so the two
#' curves run in opposite directions and cross at the symmetric baseline, where
#' the populations coincide and transferability is one.
#'
#' @param data_dir Directory holding the pre-computed results, that is the
#'   `data` directory of the repository.
#'
#' @return A `patchwork` object.
#'
#' @examples
#' \dontrun{
#' figure_01("data")
#' }
#'
#' @export
figure_01 <- function(data_dir) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    validation <- readRDS(file.path(data_dir, "theory_validation",
                                    "distinct_populations.rds"))
    scenarios <- validation$scenarios
    axis_range <- range(c(scenarios$T_theoretical, scenarios$T_empirical))

    panel_a <- ggplot2::ggplot(scenarios,
                               ggplot2::aes(.data$T_theoretical, .data$T_empirical)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = "grey30") +
        ggplot2::geom_point(colour = col_point, alpha = 0.55, size = 1.5,
                            stroke = 0) +
        ggplot2::annotate("text", x = axis_range[1], y = axis_range[2],
                          hjust = 0, vjust = 1,
                          label = sprintf("r^2 == %.3f", validation$r_squared),
                          parse = TRUE, size = 4) +
        ggplot2::labs(
            x = expression(paste("Predicted Intrinsic Transferability  ", T[theo])),
            y = expression(paste("Simulated Intrinsic Transferability  ", T[emp]))) +
        ggplot2::coord_fixed(xlim = axis_range, ylim = axis_range, clip = "off") +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(axis.title = ggplot2::element_text(size = 11),
                       plot.margin = ggplot2::margin(6, 8, 4, 6))

    panel_b <- sweep_panel(seq(0.1, 1.0, length.out = 300), "sigma2",
                           expression(paste("Additive Genetic Variance  ", sigma^2)),
                           expression(sigma[A]^2), expression(sigma[B]^2))
    panel_c <- sweep_panel(seq(0.5, 1.0, length.out = 300), "phi",
                           expression(paste("Mean Self-Kinship  ", bar(Phi))),
                           expression(bar(Phi)[A]), expression(bar(Phi)[B]))
    panel_d <- sweep_panel(seq(0.3, 3.0, length.out = 300), "tau2",
                           expression(paste("Non-Genetic Variance  ", tau^2)),
                           expression(tau[A]^2), expression(tau[B]^2))

    (panel_a | panel_b) / (panel_c | panel_d) +
        patchwork::plot_annotation(tag_levels = "a") &
        ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))
}

#' One panel of the Figure 1 parameter sweep
#'
#' Sweeps a single determinant away from a symmetric baseline at which the two
#' populations coincide, in the discovery population and then in the target
#' population, and draws both curves.
#'
#' @param grid Values of the determinant to sweep over.
#' @param determinant One of `"sigma2"`, `"phi"` or `"tau2"`.
#' @param x_label Axis label for the swept determinant.
#' @param label_A,label_B Legend labels for the two curves.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
sweep_panel <- function(grid, determinant, x_label, label_A, label_B) {
    baseline <- list(sigma2 = 0.3, phi = 0.6, tau2 = 1)

    evaluate <- function(value, population) {
        args <- list(sigma2_A = baseline$sigma2, phi_bar_A = baseline$phi,
                     tau2_A = baseline$tau2,
                     sigma2_B = baseline$sigma2, phi_bar_B = baseline$phi,
                     tau2_B = baseline$tau2)
        prefix <- switch(determinant,
                         sigma2 = "sigma2_", phi = "phi_bar_", tau2 = "tau2_",
                         stop('`determinant` must be "sigma2", "phi" or "tau2"!'))
        args[[paste0(prefix, population)]] <- value
        do.call(transferability, args)
    }

    curves <- rbind(
        data.frame(x = grid, population = "A",
                   value = vapply(grid, evaluate, numeric(1), population = "A")),
        data.frame(x = grid, population = "B",
                   value = vapply(grid, evaluate, numeric(1), population = "B")))

    ggplot2::ggplot(curves, ggplot2::aes(.data$x, .data$value,
                                         colour = .data$population,
                                         linetype = .data$population)) +
        ggplot2::geom_hline(yintercept = 1, linetype = 3, colour = "grey35",
                            linewidth = 0.9) +
        ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::scale_colour_manual(values = c(A = col_discovery, B = col_target),
                                     breaks = c("A", "B"),
                                     labels = c(label_A, label_B), name = NULL) +
        ggplot2::scale_linetype_manual(values = c(A = "longdash", B = "solid"),
                                       breaks = c("A", "B"),
                                       labels = c(label_A, label_B), name = NULL) +
        ggplot2::labs(x = x_label,
                      y = expression(paste("Transferability  ", T[A %->% B]))) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(axis.title = ggplot2::element_text(size = 11),
                       legend.position = c(0.04, 0.96),
                       legend.justification = c(0, 1),
                       legend.key = ggplot2::element_blank(),
                       legend.key.width = ggplot2::unit(1.6, "lines"),
                       legend.key.height = ggplot2::unit(1.0, "lines"),
                       legend.text = ggplot2::element_text(size = 10),
                       plot.margin = ggplot2::margin(6, 8, 4, 6))
}

#' Figure 2: transferability for three population pairs in both directions
#'
#' Each unordered population pair is one coloured line joining its two
#' directions, so the slope shows how much transferability depends on which
#' population supplies the score.  Panels are the three levels of total
#' heritability in the European sample.  Traits carry 1,000 common causal
#' variants.
#'
#' @param data_dir Directory holding the pre-computed results.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' figure_02("data")
#' }
#'
#' @export
figure_02 <- function(data_dir) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    plot_direction_pairs(read_transferability_by_heritability(
        file.path(data_dir, "transferability_common_variants",
                  sprintf("transferability_m1000_h2%s.RDS",
                          c("0.1", "0.4", "0.8")))))
}

#' Figure 3: transferability against the rare-variant fraction
#'
#' Transferability between the African and European samples as the fraction of
#' causal variants that are rare increases, at four rare-to-common effect-size
#' ratios, with total heritability 0.1 in the European sample.  The two panels
#' are the two directions of transfer.
#'
#' @param data_dir Directory holding the pre-computed results.
#'
#' @return A `patchwork` object.
#'
#' @examples
#' \dontrun{
#' figure_03("data")
#' }
#'
#' @export
figure_03 <- function(data_dir) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    arrange_rare_variant_pair(
        file.path(data_dir, "transferability_rare_variants", "h2_0.1"),
        h2_tag = "h201", first = "AFR", second = "EUR", y_max = 2.25)
}

#' Figure 4: heritability split between common and rare causal variants
#'
#' Because transferability in the coalescent setting reduces to a ratio of
#' heritabilities, this shows the quantity driving Figures 2 and 3 directly.
#' The layout is a grid: columns are the total heritability in the European
#' sample and rows are the rare-to-common effect-size ratio.  Within a block the
#' nine rare-variant fractions share one vertical scale; between blocks the
#' scale is free, because the two heritability levels differ by nearly an order
#' of magnitude.
#'
#' @param data_dir Directory holding the pre-computed results.
#'
#' @return A `patchwork` object.
#'
#' @examples
#' \dontrun{
#' figure_04("data")
#' }
#'
#' @export
figure_04 <- function(data_dir) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    heritabilities <- c("0.1", "0.8")
    ratios <- c(2, 5, 10)

    decomposition <- do.call(rbind, unlist(recursive = FALSE, lapply(
        heritabilities, function(h2) lapply(ratios, function(ratio) {
            tag <- sub("0\\.", "0", h2)
            df <- readRDS(file.path(data_dir, "heritability_decomposition",
                                    sprintf("h2decomp_h2%s_s%d.RDS", tag, ratio)))
            df$h2_level <- h2
            df$ratio <- ratio
            df
        }))))

    decomposition$pop <- factor(decomposition$pop, levels = c("ASN", "EUR", "AFR"))
    decomposition$type <- factor(decomposition$type, levels = c("rare", "common"))
    decomposition$ratio <- factor(decomposition$ratio, levels = ratios)
    decomposition$h2_level <- factor(decomposition$h2_level, levels = heritabilities)
    decomposition$ratio_title <- factor("Effect Size Ratio")
    decomposition$h2_title <- factor("Heritability in EUR")
    decomposition$fraction_title <- factor("Fraction of Rare Variants")

    # One block per heritability and effect-size ratio.  Within a block the nine
    # rare-variant fractions share a vertical scale, so the bars can be compared
    # across fractions; between blocks the scale is free, because the two
    # heritability levels differ by nearly an order of magnitude and one scale
    # would flatten the 0.1 column.  Only the top row carries the heritability
    # strip, and the effect-size ratio runs down the left.
    block <- function(heritability, ratio, top_row) {
        subset <- decomposition[decomposition$h2_level == heritability &
                                    decomposition$ratio == ratio, ]
        facets <- if (top_row)
            ggh4x::facet_nested(ratio_title + ratio ~ h2_title + h2_level +
                                    fraction_title + fraction, switch = "y")
        else
            ggh4x::facet_nested(ratio_title + ratio ~ fraction_title + fraction,
                                switch = "y")

        ggplot2::ggplot(subset, ggplot2::aes(.data$pop, .data$h2,
                                             fill = .data$type)) +
            ggplot2::geom_col(position = "stack", width = 0.8) +
            facets +
            ggplot2::scale_fill_manual("Variant Type", values = decomp_pal,
                                       breaks = c("common", "rare")) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0, 0.05))) +
            ggplot2::labs(x = "Pop", y = "Heritability") +
            ggplot2::theme_minimal(base_size = 8) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                    hjust = 1, size = 5.5),
                strip.text = ggplot2::element_text(size = 7),
                strip.background = ggplot2::element_rect(fill = "grey85",
                                                         colour = NA),
                strip.placement = "outside",
                panel.spacing = ggplot2::unit(1.5, "pt"))
    }

    rows <- lapply(seq_along(ratios), function(i) {
        left <- block(heritabilities[1], ratios[i], i == 1)
        right <- block(heritabilities[2], ratios[i], i == 1)
        (left | right) +
            patchwork::plot_layout(guides = "collect") &
            ggplot2::theme(legend.position = "right")
    })

    patchwork::wrap_plots(rows, ncol = 1)
}
