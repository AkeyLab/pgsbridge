# simulation/08_allele_frequencies.R
#
# Computes the per-population allele frequency of every simulated variant, split
# into the common and the rare pool. Produces the data behind S1, S7 and S8 Figs.
#
# These three figures are all about where a variant's frequency sits in one
# population relative to another, and in particular about variants that
# segregate in the African sample and have fallen to zero frequency elsewhere.
# All of them need nothing more than the frequencies themselves, so they are
# computed once here and read by each figure script.
#
# Output
# ------
#   common_af.RDS   a list of three vectors, one per population, holding the
#                   reference-allele frequency of every common variant
#   rare_af.RDS     the same for the rare pool
#   simulated_common_variants.RDS   the common-pool frequencies again, named
#                   af_eur_001, af_afr_001 and af_asi_001, which is what S1 Fig
#                   reads
#
# In all of these files the vectors are aligned: position k refers to the same
# variant in all three populations.
#
# Usage:  Rscript simulation/08_allele_frequencies.R

source("config.R")
source(file.path(PRS_ROOT, "R", "traits.R"))
source(file.path(PRS_ROOT, "R", "variant_pools.R"))

pools <- load_variant_pools(PRS_GENOTYPES)

cat(sprintf("%d common and %d rare variants\n",
            length(pools$common_index), length(pools$rare_index)))

output_dir <- file.path(PRS_DATA, "allele_frequencies")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# The rare pool runs to millions of variants, so the frequencies are
# accumulated in chunks rather than by loading the genotypes.
frequencies <- list()
for (pool in c("common", "rare")) {
  frequencies[[pool]] <- lapply(c(AFR = "AFR", EUR = "EUR", ASN = "ASN"),
                                function(population) {
    cat(sprintf("  computing %s frequencies for %s\n", pool, population))
    flush.console()
    pool_allele_frequencies(pools, population, pool)
  })
}

saveRDS(list(pafr_common = frequencies$common$AFR,
             peur_common = frequencies$common$EUR,
             pasi_common = frequencies$common$ASN),
        file.path(output_dir, "common_af.RDS"))

saveRDS(list(pafr_rare = frequencies$rare$AFR,
             peur_rare = frequencies$rare$EUR,
             pasi_rare = frequencies$rare$ASN),
        file.path(output_dir, "rare_af.RDS"))

# The simulated frequencies used by S1 Fig, over the common pool.
saveRDS(list(af_eur_001 = frequencies$common$EUR,
             af_afr_001 = frequencies$common$AFR,
             af_asi_001 = frequencies$common$ASN),
        file.path(output_dir, "simulated_common_variants.RDS"))

# A short summary, printed so that the numbers quoted in the text can be checked
# without opening the files.
minor <- function(p) pmin(p, 1 - p)

cat("\nMean minor allele frequency over the common pool\n")
for (population in c("AFR", "EUR", "ASN")) {
  cat(sprintf("  %s  %.4f\n", population,
              mean(minor(frequencies$common[[population]]))))
}

cat("\nPercentage of variants segregating in AFR and absent elsewhere\n")
for (target in c("EUR", "ASN")) {
  for (pool in c("common", "rare")) {
    african <- frequencies[[pool]]$AFR
    other <- frequencies[[pool]][[target]]
    segregating <- african > 0 & african < 1
    monomorphic <- other == 0 | other == 1
    cat(sprintf("  %-6s AFR to %s: %.2f%%\n", pool, target,
                100 * sum(segregating & monomorphic) / length(african)))
  }
}

cat("\nwritten to", output_dir, "\n")
