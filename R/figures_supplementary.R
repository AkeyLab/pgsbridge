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
