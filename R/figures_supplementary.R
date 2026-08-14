#' S1 Fig: transferability is robust to the number of causal variants
#'
#' The quantity of [figure_02()], which uses 1,000 causal variants, repeated at
#' 20, 50 and 100 causal variants.  The pattern is unchanged, so the number of
#' causal variants has only a modest effect on transferability.
#'
#' @param data_dir Directory holding the pre-computed results.
#'
#' @return A `patchwork` object.
#'
#' @examples
#' \dontrun{
#' figure_S1("data")
#' }
#'
#' @export
figure_S1 <- function(data_dir) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    panels <- lapply(c(20, 50, 100), function(m) {
        plot_direction_pairs(read_transferability_by_heritability(
            file.path(data_dir, "transferability_common_variants",
                      sprintf("transferability_m%d_h2%s.RDS", m,
                              c("0.1", "0.4", "0.8")))))
    })

    # Each panel keeps its own legend and the panels are tagged a, b, c, which is
    # how the figure appears in the manuscript.
    patchwork::wrap_plots(panels, ncol = 1) +
        patchwork::plot_annotation(tag_levels = "a") &
        ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))
}

#' S2 Fig: transferability is robust to the effect-size distribution
#'
#' The upper panel repeats [figure_02()], in which every causal variant is given
#' the same effect size.  The lower panel repeats the experiment with effect
#' sizes drawn from a normal distribution.  The pattern is unchanged, so the
#' result does not depend on the effects being equal.
#'
#' @param data_dir Directory holding the pre-computed results.
#'
#' @return A `patchwork` object.
#'
#' @examples
#' \dontrun{
#' figure_S2("data")
#' }
#'
#' @export
figure_S2 <- function(data_dir) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    heritabilities <- c("0.1", "0.4", "0.8")
    fixed <- plot_direction_pairs(read_transferability_by_heritability(
        file.path(data_dir, "transferability_common_variants",
                  sprintf("transferability_m1000_h2%s.RDS", heritabilities))))
    normal <- plot_direction_pairs(read_transferability_by_heritability(
        file.path(data_dir, "transferability_normal_effects",
                  sprintf("transferability_m1000_h2%s.RDS", heritabilities))))

    # No tags and one legend per panel: the caption names the top and the bottom
    # panel rather than (a) and (b), which is how the figure appears in the
    # manuscript.
    fixed / normal
}

#' S3 Fig: average minor allele frequency in the three populations
#'
#' The ordering of the average minor allele frequency across populations is what
#' drives every transferability result in the paper: the population whose shared
#' causal variants segregate at higher frequency carries more additive genetic
#' variance, and so a higher ceiling.  This figure checks that the coalescent
#' simulation reproduces the ordering seen in real data, by drawing the same
#' quantity for the simulation and for the 1000 Genomes Project.
#'
#' Both series use the pooled definition of a common variant: minor allele
#' frequency at least `maf_common` in the three populations taken together.
#' This is the definition the trait simulations draw from, and it is the one
#' used throughout the paper.  It matters which is used, because requiring a
#' variant to be common in every population separately selects against variation
#' private to any one of them and reverses the ordering.
#'
#' The 1000 Genomes series is read from alternate-allele counts over the
#' African, European and East Asian superpopulations, computed from the phase 3
#' genotypes by `scripts/empirical/1000g_frequencies.sh`.
#'
#' @param data_dir Directory holding the pre-computed results.
#' @param maf_common The minor allele frequency defining a common variant.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' figure_S3("data")
#' }
#'
#' @export
figure_S3 <- function(data_dir, maf_common = 0.01) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    populations <- c("AFR", "EUR", "ASN")
    minor <- function(p) pmin(p, 1 - p)

    # The simulated common pool is already the pooled-common set: every
    # population contributes 1,000 diploids, so the pooled frequency is the mean
    # of the three and the pool was filtered on it.
    simulated <- readRDS(file.path(data_dir, "allele_frequencies",
                                   "common_af.RDS"))
    simulated_mean <- c(
        AFR = mean(minor(simulated$pafr_common)),
        EUR = mean(minor(simulated$peur_common)),
        ASN = mean(minor(simulated$pasi_common)))

    # 1000 Genomes: alternate-allele counts per superpopulation, over biallelic
    # single-nucleotide variants.  The superpopulations differ in size, so the
    # pooled frequency is a count ratio rather than a mean of the three.
    empirical <- readRDS(file.path(data_dir, "allele_frequencies",
                                   "thousand_genomes_common_af.RDS"))
    counts <- empirical$counts
    allele_number <- empirical$allele_number
    pooled <- rowSums(counts) / sum(allele_number)
    keep <- minor(pooled) >= maf_common
    empirical_mean <- vapply(c(AFR = "AFR", EUR = "EUR", ASN = "EAS"),
                             function(column) {
                                 mean(minor(counts[keep, column] /
                                            allele_number[[column]]))
                             }, numeric(1))

    summary_df <- rbind(
        data.frame(population = populations, mean_maf = simulated_mean[populations],
                   series = "Simulation", stringsAsFactors = FALSE),
        data.frame(population = populations, mean_maf = empirical_mean[populations],
                   series = "1000 Genome Project", stringsAsFactors = FALSE))
    summary_df$population <- factor(summary_df$population, levels = populations)
    summary_df$series <- factor(summary_df$series,
                                levels = c("Simulation", "1000 Genome Project"))

    ggplot2::ggplot(summary_df,
                    ggplot2::aes(x = .data$population, y = .data$mean_maf,
                                 colour = .data$series, group = .data$series)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2.5, shape = 21, fill = "white",
                            stroke = 1.1) +
        ggplot2::scale_colour_manual(
            values = unname(okabe_ito[c("vermillion", "blue")]), name = NULL) +
        ggplot2::labs(x = "Pop", y = "Minor Allele Frequency") +
        ggplot2::expand_limits(y = 0) +
        ggplot2::theme_minimal(base_size = PGS_BASE_SIZE) +
        ggplot2::theme(legend.position = c(0.98, 0.04),
                       legend.justification = c(1, 0))
}


#' S4 to S6 Figs: rare-variant sweeps for the remaining population pairs
#'
#' Figure 3 shows the African and European pair at heritability 0.1.  These show
#' the other pairs, and the other heritabilities.  `figure_S4()` gives the
#' African with East Asian and European with East Asian pairs at heritability
#' 0.1; `figure_S5()` and `figure_S6()` give all three pairs at heritability 0.4
#' and 0.8 respectively.
#'
#' @param data_dir Directory holding the pre-computed results.
#'
#' @return A `patchwork` object.
#'
#' @examples
#' \dontrun{
#' figure_S4("data")
#' }
#'
#' @name figures_rare_variants
NULL

#' @rdname figures_rare_variants
#' @export
figure_S4 <- function(data_dir) {
    rare_variant_figure(data_dir, "0.1", "h201",
                        list(c("AFR", "ASN"), c("EUR", "ASN")))
}

#' @rdname figures_rare_variants
#' @export
figure_S5 <- function(data_dir) {
    rare_variant_figure(data_dir, "0.4", "h204",
                        list(c("AFR", "EUR"), c("AFR", "ASN"), c("EUR", "ASN")))
}

#' @rdname figures_rare_variants
#' @export
figure_S6 <- function(data_dir) {
    rare_variant_figure(data_dir, "0.8", "h208",
                        list(c("AFR", "EUR"), c("AFR", "ASN"), c("EUR", "ASN")))
}

#' Assemble a rare-variant sweep figure from a list of population pairs
#'
#' Population pairs are laid out two to a row, each keeping its own legend, so a
#' figure with three pairs fills the first row and puts the third below.
#' `wrap_elements` keeps `patchwork` from tagging every leaf panel of the nested
#' arrangements, so each population pair receives exactly one tag.
#'
#' @param data_dir Directory holding the pre-computed results.
#' @param heritability The heritability level, as it appears in the directory
#'   name.
#' @param h2_tag The heritability tag as it appears in the file names.
#' @param pairs A list of two-element character vectors naming the population
#'   pairs.
#' @param y_max Upper limit of the vertical axis, shared by every panel.
#'
#' @return A `patchwork` object.
#'
#' @keywords internal
rare_variant_figure <- function(data_dir, heritability, h2_tag, pairs,
                                y_max = 2.5) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    directory <- file.path(data_dir, "transferability_rare_variants",
                           paste0("h2_", heritability))
    panels <- lapply(pairs, function(pair) {
        arrange_rare_variant_pair(directory, h2_tag, pair[1], pair[2], y_max)
    })

    patchwork::wrap_plots(lapply(panels, patchwork::wrap_elements), ncol = 2) +
        patchwork::plot_annotation(tag_levels = "a") &
        ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))
}

#' S7 Fig: allele frequency in African against non-African samples
#'
#' Four panels in one row: the African with European comparison for common then
#' rare variants, followed by the African with East Asian comparison for the
#' same two classes.  The point of the figure is the mass of rare variants lying
#' on the horizontal axis: variants that segregate in the African sample and have
#' fallen to zero frequency in the non-African one.  Rare variants are lost this
#' way far more often than common ones, which is what makes the heritability gap
#' widen as more of the causal signal is placed on rare variants.
#'
#' This function draws with base graphics onto the device that is already open,
#' so it returns nothing.  Open a device first.
#'
#' @param data_dir Directory holding the pre-computed results.
#' @param rare_axis_limit Upper limit of both axes in the rare-variant panels,
#'   which occupy only a small corner of the unit square.
#'
#' @return Nothing; called for the plot it draws.
#'
#' @examples
#' \dontrun{
#' pdf("S7_Fig.pdf", height = 8, width = 8)
#' figure_S7("data")
#' dev.off()
#' }
#'
#' @export
figure_S7 <- function(data_dir, rare_axis_limit = 0.02) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    common <- readRDS(file.path(data_dir, "allele_frequencies", "common_af.RDS"))
    rare <- readRDS(file.path(data_dir, "allele_frequencies", "rare_af.RDS"))
    density_ramp <- grDevices::colorRampPalette(
        c("white", unname(okabe_ito["blue"])))

    old <- graphics::par(mfrow = c(1, 4), mar = c(4.2, 4.2, 2.4, 1.2))
    on.exit(graphics::par(old), add = TRUE)

    graphics::smoothScatter(common$pafr_common, common$peur_common,
                            xlab = "AFR Common", ylab = "EUR Common",
                            colramp = density_ramp)
    graphics::mtext("a", side = 3, adj = 0, line = 0.8, font = 2, cex = 1.2)
    graphics::smoothScatter(rare$pafr_rare, rare$peur_rare,
                            xlim = c(0, rare_axis_limit),
                            ylim = c(0, rare_axis_limit),
                            xlab = "AFR Rare", ylab = "EUR Rare",
                            colramp = density_ramp)
    graphics::smoothScatter(common$pafr_common, common$pasi_common,
                            xlab = "AFR Common", ylab = "ASN Common",
                            colramp = density_ramp)
    graphics::mtext("b", side = 3, adj = 0, line = 0.8, font = 2, cex = 1.2)
    graphics::smoothScatter(rare$pafr_rare, rare$pasi_rare,
                            xlim = c(0, rare_axis_limit),
                            ylim = c(0, rare_axis_limit),
                            xlab = "AFR Rare", ylab = "ASN Rare",
                            colramp = density_ramp)
    invisible(NULL)
}

#' S8 Fig: African variants that are monomorphic elsewhere
#'
#' The percentage of variants segregating in the African sample that have become
#' monomorphic in a non-African one, for common and for rare variants.  This puts
#' a number on the pattern [figure_S7()] shows graphically.
#'
#' The stored frequencies count the reference allele, which is held fixed across
#' populations, so a variant monomorphic in the target has frequency zero or one
#' there.  Both contribute exactly zero additive genetic variance, so both count
#' as lost.
#'
#' @param data_dir Directory holding the pre-computed results.
#'
#' @return A `patchwork` object.
#'
#' @examples
#' \dontrun{
#' figure_S8("data")
#' }
#'
#' @export
figure_S8 <- function(data_dir) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    common <- readRDS(file.path(data_dir, "allele_frequencies", "common_af.RDS"))
    rare <- readRDS(file.path(data_dir, "allele_frequencies", "rare_af.RDS"))
    type_pal <- unname(decomp_pal[c("common", "rare")])
    names(type_pal) <- c("Common", "Rare")

    # A variant counts as lost when it is present in the African sample and has
    # drifted to zero frequency in the other one.  Variants that have instead
    # gone to fixation in the other sample are not counted, which is what the
    # published figure does; note that they too contribute no genetic variance
    # there, so this is the narrower of the two reasonable rules.
    percent_lost <- function(african, other) {
        round(100 * sum(african > 0 & other == 0) / length(african), 2)
    }

    loss_panel <- function(percent_common, percent_rare, target, tag) {
        df <- data.frame(Type = factor(c("Common", "Rare"),
                                       levels = c("Common", "Rare")),
                         percent = c(percent_common, percent_rare))
        ggplot2::ggplot(df, ggplot2::aes(x = .data$Type, y = .data$percent,
                                         fill = .data$Type)) +
            ggplot2::geom_col(width = 0.6) +
            ggplot2::geom_text(ggplot2::aes(label = .data$percent),
                               vjust = -0.4, size = 3.5) +
            ggplot2::scale_fill_manual(values = type_pal) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0, 0.12))) +
            ggplot2::labs(
                x = "Type",
                y = sprintf("Allele frequency %% in AFR becomes 0 in %s", target),
                title = sprintf("From AFR to %s", target), tag = tag) +
            ggplot2::theme_minimal(base_size = PGS_BASE_SIZE) +
            ggplot2::theme(legend.position = "none",
                           plot.title = ggplot2::element_text(hjust = 0.5),
                           plot.tag = ggplot2::element_text(face = "bold",
                                                            size = 14))
    }

    loss_panel(percent_lost(common$pafr_common, common$peur_common),
               percent_lost(rare$pafr_rare, rare$peur_rare), "EUR", "a") |
        loss_panel(percent_lost(common$pafr_common, common$pasi_common),
                   percent_lost(rare$pafr_rare, rare$pasi_rare), "ASN", "b")
}
