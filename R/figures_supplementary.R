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

#' S7 Fig: rare variants are lost outside Africa more often than common ones
#'
#' The rare-variant results rest on a single asymmetry: a causal variant that
#' segregates in the African sample is far more likely to have drifted to zero
#' frequency outside Africa than a common one, and a variant at zero frequency
#' contributes no additive genetic variance at all in that population. This
#' figure states that asymmetry directly, for both variant pools.
#'
#' Earlier versions of the manuscript carried two figures here, a smoothed
#' scatter of the joint allele-frequency distribution and a bar chart of the
#' variants lost from Africa outward. Neither stated the claim the text makes,
#' so both were replaced by this one.
#'
#' @param data_dir Directory holding the pre-computed results.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' figure_S7("data")
#' }
#'
#' @export
figure_S7 <- function(data_dir) {
    if (missing(data_dir))
        stop('`data_dir` is required!')

    common <- readRDS(file.path(data_dir, "allele_frequencies", "common_af.RDS"))
    rare <- readRDS(file.path(data_dir, "allele_frequencies", "rare_af.RDS"))

    pools <- list(
        Common = list(AFR = common$pafr_common, EUR = common$peur_common,
                      ASN = common$pasi_common),
        Rare = list(AFR = rare$pafr_rare, EUR = rare$peur_rare,
                    ASN = rare$pasi_rare))
    populations <- c("AFR", "EUR", "ASN")

    absent <- do.call(rbind, lapply(names(pools), function(pool) {
        do.call(rbind, lapply(populations, function(population) {
            data.frame(pool = pool, population = population,
                       percent = 100 * mean(pools[[pool]][[population]] == 0))
        }))
    }))
    absent$population <- factor(absent$population, levels = populations)
    absent$pool <- factor(absent$pool, levels = c("Common", "Rare"))

    ggplot2::ggplot(absent, ggplot2::aes(.data$population, .data$percent,
                                         fill = .data$population)) +
        ggplot2::geom_col(width = 0.65) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", .data$percent)),
                           vjust = -0.45, size = 3.4) +
        ggplot2::facet_wrap(~ pool, scales = "free_y") +
        ggplot2::scale_fill_manual(
            values = unname(okabe_ito[c("vermillion", "blue", "green")])) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.16))) +
        ggplot2::labs(x = "Population",
                      y = "Percentage of variants at zero frequency") +
        ggplot2::theme_minimal(base_size = PGS_BASE_SIZE) +
        ggplot2::theme(legend.position = "none",
                       strip.background = ggplot2::element_rect(fill = "grey92",
                                                                colour = NA),
                       strip.text = ggplot2::element_text(face = "bold"))
}

#' S8 Fig: intrinsic transferability under five causal-variant ascertainment schemes
#'
#' Figure 3 decides which causal variants are common and which are rare on the
#' three populations pooled. That is a modelling choice, and a variant that is
#' rare in the pooled sample need not be rare in any single population. This
#' figure repeats the Figure 3 sweep under five rules for making that decision.
#'
#' Writing `m_S` for the minor allele count in population `S`, out of the 2,000
#' alleles it contributes, and `m_global` for the count out of all 6,000, the
#' threshold of 0.01 is 20 within a population and 60 globally. The rules are
#'
#' \describe{
#'   \item{Model 1, all three}{common when `m_S >= 20` in every population, rare
#'     when `m_S < 20` in every population. A variant may be absent in one or two
#'     populations and still count as rare. Variants common in one population and
#'     rare in another belong to neither pool.}
#'   \item{Model 2, European}{the rule applied to the European sample alone,
#'     among the variants that segregate there.}
#'   \item{Model 3, African}{the same, applied to the African sample.}
#'   \item{Model 4a, pooled}{common when `m_global >= 60`. This is the rule used
#'     for Figure 3, and it is also the global allele-frequency rule, since
#'     pooling the samples and computing one frequency over 6,000 alleles is the
#'     same calculation.}
#'   \item{Model 4b, random}{no ascertainment at all. Causal variants are drawn
#'     uniformly from the whole panel, and the Model 4a split is applied
#'     afterwards only to decide which of them receive the larger effect size.}
#' }
#'
#' Ascertaining on European frequencies roughly doubles the asymmetry between the
#' two directions of transfer, requiring the same class in all three populations
#' removes it, random sampling reproduces the pooled result, and ascertaining on
#' African frequencies reverses its sign. In every scheme the value is the one
#' the closed form gives from the additive genetic variance the ascertained
#' variants carry in each population, so the ascertainment sets which population
#' the mechanism favours rather than whether it operates.
#'
#' The first four schemes choose what fraction of the causal variants is rare, so
#' that fraction is their horizontal axis. Model 4b does not choose it: drawing
#' without reference to frequency realises a rare fraction near 0.88 whatever
#' else is done, so plotting it on the same axis would put its points on top of
#' one another. It is therefore drawn underneath, on its own axis, against the
#' number of causal variants. The vertical axis is shared throughout.
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

    scheme_labels <- c(model1_allthree = "Model 1 all three",
                       model2_european = "Model 2 European",
                       model3_african  = "Model 3 African",
                       model4a_pooled  = "Model 4a pooled",
                       model4b_random  = "Model 4b random")

    results <- utils::read.csv(file.path(data_dir, "ascertainment",
                                         "transferability_by_scheme.csv"))
    results <- results[results$pop_comb %in% c("AFR-EUR", "EUR-AFR"), ]

    results$ratio <- factor(results$ratio, levels = c(1, 2, 5, 10))
    results$direction <- factor(
        ifelse(results$pop_comb == "AFR-EUR", "AFR to EUR", "EUR to AFR"),
        levels = c("AFR to EUR", "EUR to AFR"))
    results$scheme <- factor(unname(scheme_labels[results$model]),
                             levels = unname(scheme_labels))

    y_max <- 3.6
    hidden <- results$emp_tau + results$ci > y_max | results$emp_tau - results$ci < 0
    if (any(hidden))
        stop("the vertical limit would hide ", sum(hidden), " point(s)")

    random <- results$model == "model4b_random"
    ascertained <- results[!random, ]
    unascertained <- results[random, ]

    common_layers <- function(plot) {
        plot +
            ggplot2::geom_hline(yintercept = 1, linetype = "dotted",
                                colour = "black") +
            ggplot2::scale_colour_manual(values = ratio_pal) +
            ggplot2::coord_cartesian(ylim = c(0, y_max)) +
            ggplot2::labs(colour = "Rare/Common Effect Size Ratio",
                          y = "Intrinsic Transferability") +
            ggplot2::theme_minimal(base_size = PGS_BASE_SIZE) +
            ggplot2::theme(legend.position = "bottom",
                           legend.margin = ggplot2::margin(0, 0, 0, 0),
                           panel.spacing = ggplot2::unit(0.5, "lines"),
                           strip.text.y = ggplot2::element_text(
                               size = PGS_BASE_SIZE - 2))
    }

    # The four schemes that choose a rare fraction, on that shared axis.
    upper <- common_layers(
        ggplot2::ggplot(ascertained,
                        ggplot2::aes(.data$fraction, .data$emp_tau,
                                     colour = .data$ratio))) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data$emp_tau - .data$ci,
                         ymax = .data$emp_tau + .data$ci), width = 0.06) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::geom_point(size = 1.0) +
        ggplot2::facet_grid(scheme ~ direction) +
        ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                                    limits = c(0, 1)) +
        ggplot2::labs(x = "Fraction of Rare Variants")

    # Model 4b, on its own axis, because it realises a rare fraction rather than
    # choosing one.
    realised <- sprintf("%.3f", mean(unascertained$rare_share))
    lower <- common_layers(
        ggplot2::ggplot(unascertained,
                        ggplot2::aes(.data$n_causal, .data$emp_tau,
                                     colour = .data$ratio))) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data$emp_tau - .data$ci,
                         ymax = .data$emp_tau + .data$ci), width = 0.05) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::geom_point(size = 1.0) +
        ggplot2::facet_grid(scheme ~ direction) +
        ggplot2::scale_x_log10(
            breaks = sort(unique(unascertained$n_causal)),
            labels = format(sort(unique(unascertained$n_causal)),
                            big.mark = ",", trim = TRUE)) +
        ggplot2::labs(x = paste0("Number of Causal Variants (realised rare ",
                                 "fraction ", realised, ")"))

    patchwork::wrap_plots(upper, lower, ncol = 1, heights = c(4, 1.35)) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
}
