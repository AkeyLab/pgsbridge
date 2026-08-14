#' What can and cannot be reproduced from this repository
#'
#' Most of the manuscript can be regenerated from what is distributed here.
#' Four things cannot, and this page records exactly which, and why, so that a
#' reader is never left guessing whether a mismatch is their mistake.
#'
#' @section What reproduces:
#'
#' All twelve figures redraw from the replicate-level results in `data/` by
#' calling [figure_01()] through [figure_S8()], with no genotypes and no
#' cluster.  The simulation scripts regenerate those results from the coalescent
#' genotypes: [figure_02()]'s inputs have been checked directly against the
#' distributed files and agree within Monte Carlo error at every heritability.
#'
#' @section S3 Fig cannot be reproduced:
#'
#' The published S3 Fig draws two series, one from the coalescent simulation and
#' one from the 1000 Genomes Project.  [figure_S3()] draws only the first.
#'
#' The reason is that no 1000 Genomes allele frequencies are distributed with
#' this repository, and no script here computes them.  The published series
#' therefore cannot be regenerated from what is here; it would have to be
#' recomputed from the 1000 Genomes Project data itself.
#'
#' The simulated series is also not an exact match.  Read off the published
#' figure it runs at roughly 0.125, 0.089 and 0.086 for the African, European
#' and East Asian samples.  From the distributed frequencies the pooled
#' definition of a common variant gives 0.1244, 0.1101 and 0.1034, and the
#' per-population definition gives 0.2152, 0.2398 and 0.2334.  The African value
#' matches the pooled definition closely and the other two match neither, so a
#' third variant definition lies behind the published curve.  Which one is not
#' recorded anywhere available.
#'
#' Note that the two definitions order the populations differently.  Under the
#' pooled definition the African sample has the highest average minor allele
#' frequency; requiring a variant to be common in every population selects
#' against variation private to any one of them and reverses that.  The trait
#' simulations draw from the pooled panel.
#'
#' @section The published Figure 2 image predates the distributed results:
#'
#' [figure_02()] reproduces the *data* behind Figure 2, but not the *image* in
#' the manuscript.  At a heritability of 0.1 the published labels read 0.996,
#' 0.979 and 0.970 for the three pairs, whereas the distributed replicates give
#' 1.292, 1.231 and 1.128.  At a heritability of 0.8 the two agree to three
#' decimal places, and at 0.4 they agree in one direction only.
#'
#' The reason is that the published image was drawn from an earlier run and was
#' never redrawn after the replicates were regenerated.  The code and the data
#' here are consistent with each other; it is the manuscript image that is out
#' of date.
#'
#' @section Figure 1a is sensitive to the package versions:
#'
#' `scripts/simulation/01_theory_validation.R` seeds every scenario
#' individually, so it does not depend on the number of cores, but rerunning it
#' does not reproduce the distributed
#' `data/theory_validation/distinct_populations.rds` exactly.  A rerun keeps 385
#' of 400 scenarios and gives a Pearson `r^2` of 0.99771, against 382 and
#' 0.99747 in the distributed file.  None of the regenerated scenarios matches a
#' distributed one to within `1e-6`, so this is a different draw rather than a
#' rounding difference.
#'
#' The reason is that no package versions are pinned, and the random streams of
#' `bnpsd` or of R itself have changed since the distributed file was produced.
#' The consequence is small but visible: panel a would annotate `r^2 = 0.998`
#' where the manuscript states 0.997.
#'
#' @section The coalescent genotypes cannot be regenerated exactly:
#'
#' `scripts/coalescent/01_simulate_genotypes.py` accepts a `--seed`, but the run
#' that produced the genotypes behind Figures 2 to 4 and S1 to S8 was made
#' before that option existed and its seed was not recorded.  Any rerun
#' therefore produces a statistically equivalent but numerically different set
#' of genotypes, and every downstream result inherits that.  Results here are
#' reproducible in distribution, never bit for bit.
#'
#' @section The rare-variant pool is only approximately recovered:
#'
#' `scripts/simulation/04_traits_rare_and_common.R` and
#' `scripts/simulation/05_heritability_decomposition.R` define the rare pool as
#' the complement of the pooled common panel, which is the definition the
#' Methods give and which yields 2,817,430 variants, matching the fixation-index
#' analysis.  The surviving intermediate files from the published run instead
#' hold 907,919 rare variants, so a further filter was applied whose exact form
#' is not recorded.
#'
#' Rerunning the reconstruction at the Figure 3 setting reproduces the published
#' asymmetry between African and non-African samples, and four of the six
#' ordered pairs agree within Monte Carlo error, with the remaining two about
#' 2.3 standard errors away.  The qualitative result is therefore reproduced and
#' the precise values are not.
#'
#' @name reproducibility
#' @keywords internal
NULL
