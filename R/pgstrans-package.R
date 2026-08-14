#' pgstrans: intrinsic transferability of polygenic scores between populations
#'
#' A polygenic score built in one population usually predicts less accurately in
#' another.  This package separates the part of that loss which is fixed by
#' population history from the part better data and methods could still recover,
#' by computing the best accuracy any score could reach.
#'
#' The best accuracy attainable in a population `S` is the share of its trait
#' variance that is genetic,
#'
#' `PA_S = 2 sigma_S^2 Phibar_S / (2 sigma_S^2 Phibar_S + tau_S^2)`,
#'
#' set by three quantities: the additive genetic variance `sigma_S^2` fixed by
#' the allele-frequency spectrum, the average self-kinship `Phibar_S` summarising
#' within-population structure, and the non-genetic variance `tau_S^2`.  Because
#' the optimal predictor is given rather than estimated, this is a ceiling that no
#' score built from data can exceed.  Intrinsic transferability from a discovery
#' population `A` to a target population `B` is the ratio of their ceilings,
#' `PA_B / PA_A`.
#'
#' The functions fall into four groups.  [prediction_accuracy()] and
#' [transferability()] evaluate the closed form, with
#' [additive_genetic_variance()], [shared_effect_size()] and
#' [split_effect_sizes()] supplying its inputs.  [simulate_trait()] and
#' [empirical_transferability()] measure the same quantities from simulated
#' traits.  [load_variant_pools()] and [draw_causal_genotypes()] read genotypes
#' from PLINK files without loading pools of millions of variants into memory.
#' [group_summary()], [plot_direction_pairs()] and [plot_rare_variant_sweep()]
#' turn replicate output into the figures.
#'
#' @keywords internal
#'
#' @importFrom stats cor median qt quantile rnorm runif sd var
#' @importFrom utils head
#' @importFrom rlang .data :=
#' @import ggplot2
#' @import ggh4x
#' @import patchwork
#' @import dplyr
"_PACKAGE"

# Column names used inside dplyr and ggplot2 calls, which R CMD check would
# otherwise report as undefined global variables.
utils::globalVariables(c(
    "se", ".mean", "n", "sd", "ci", "emp_tau", "pop_comb", "pop_pair",
    "direction", "h2", "fraction", "ratio", "x", "on_left", "label",
    "label_y", "heritability_title"
))
