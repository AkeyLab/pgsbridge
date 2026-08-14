POPULATION_STEMS <- c(AFR = "afr", EUR = "eur", ASN = "asn")

#' Open the genotype panels and locate the two variant pools
#'
#' @param genotype_dir Directory holding the PLINK files written by
#'   02_prepare_genotypes.sh.
#' @return A list with
#'   `genotypes`, a named list of three `BEDMatrix` handles on the full panels;
#'   `common_index` and `rare_index`, column indices into those panels;
#'   `n_individuals`, the number of individuals per population.
#' @export
load_variant_pools <- function(genotype_dir) {
    for (pkg in c('BEDMatrix', 'genio'))
        if (!requireNamespace(pkg, quietly = TRUE))
            stop('package `', pkg, '` is required to read genotypes!')

  genotypes <- list()
  positions <- list()

  for (population in names(POPULATION_STEMS)) {
    stem <- POPULATION_STEMS[[population]]
    genotypes[[population]] <- BEDMatrix::BEDMatrix(
      file.path(genotype_dir, paste0(stem, "_all_variants")))
    positions[[population]] <- genio::read_bim(
      file.path(genotype_dir, paste0(stem, "_all_variants.bim")))$pos
  }

  # The three full panels are cut from one pooled file, so they must carry the
  # same variants in the same order. Check rather than assume.
  for (population in names(POPULATION_STEMS)[-1]) {
    if (!identical(positions[[population]], positions[[1]])) {
      stop("the three populations do not share one variant order, so a column ",
           "index would not refer to the same variant in each")
    }
  }

  common_positions <- genio::read_bim(
    file.path(genotype_dir,
              paste0(POPULATION_STEMS[[1]], "_common.bim")))$pos

  common_index <- match(common_positions, positions[[1]])
  if (anyNA(common_index)) {
    stop("a common variant has no matching position in the full panel")
  }
  rare_index <- setdiff(seq_along(positions[[1]]), common_index)

  list(genotypes = genotypes,
       common_index = sort(common_index),
       rare_index = rare_index,
       n_individuals = nrow(genotypes[[1]]))
}

#' Draw causal variants from one pool and read only those columns
#'
#' @param pools The value returned by `load_variant_pools()`.
#' @param pool Either "common" or "rare".
#' @param m_causal How many causal variants to draw.
#' @return A list with `index`, the drawn column indices into the full panel,
#'   and `genotypes`, a named list of three matrices holding those columns.
#' @export
draw_causal_genotypes <- function(pools, pool, m_causal) {
  available <- switch(pool,
                      common = pools$common_index,
                      rare = pools$rare_index,
                      stop("pool must be \"common\" or \"rare\""))
  if (m_causal == 0) {
    empty <- lapply(pools$genotypes,
                    function(X) matrix(numeric(0), nrow = pools$n_individuals))
    return(list(index = integer(0), genotypes = empty))
  }
  index <- available[sample_causal_variants(length(available), m_causal)]
  list(index = index,
       genotypes = lapply(pools$genotypes, function(X) X[, index, drop = FALSE]))
}

#' Allele frequencies over a whole pool, read in chunks
#'
#' Computes the reference-allele frequency of every variant in a pool without
#' ever holding the pool in memory.
#'
#' @param pools The value returned by `load_variant_pools()`.
#' @param population Which population to compute frequencies for.
#' @param pool Either "common" or "rare".
#' @param chunk_size How many variants to read at a time.
#' @return A numeric vector of allele frequencies, in pool order.
#' @export
pool_allele_frequencies <- function(pools, population, pool,
                                    chunk_size = 20000L) {
  index <- switch(pool,
                  common = pools$common_index,
                  rare = pools$rare_index,
                  stop("pool must be \"common\" or \"rare\""))
  X <- pools$genotypes[[population]]
  frequencies <- numeric(length(index))
  for (start in seq(1L, length(index), by = chunk_size)) {
    stop_at <- min(start + chunk_size - 1L, length(index))
    block <- X[, index[start:stop_at], drop = FALSE]
    frequencies[start:stop_at] <- colSums(block) / (2 * nrow(block))
  }
  frequencies
}
