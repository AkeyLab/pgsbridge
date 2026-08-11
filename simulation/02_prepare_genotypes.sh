#!/bin/bash
#
# Turn the VCF written by 01_simulate_genotypes.py into the PLINK files that the
# trait simulations and the fixation-index calculation read.
#
# What this produces
# ------------------
#   all_populations                 all 3,000 individuals, all variants
#   all_populations_common          all 3,000 individuals, common variants only
#   {afr,eur,asn}_all_variants      one population, all variants
#   {afr,eur,asn}_common            one population, common variants only
#
# Definitions
# -----------
# Common variants are those with minor allele frequency at least 0.01 in the
# three populations pooled. Rare variants are the remaining segregating
# variants, that is those whose pooled minor allele frequency is below 0.01.
# The rare pool is formed later, inside the trait-simulation scripts, as the
# complement of the common pool on the pooled sample, with no per-population
# condition, so a rare variant may be monomorphic in one or more populations.
#
# Every per-population file is extracted with --keep-allele-order, so that a
# variant's reference allele is the same in all three populations. Without this
# PLINK would reset the coded allele to whichever is minor within each
# population, and allele frequencies would not be comparable across them.
#
# The simulator does not assign variant identifiers, so every variant carries
# the identifier ".". Downstream matching between the full and the common panel
# is done by base-pair position, not by identifier.
#
# Requirements: PLINK 1.90 and PLINK 2.00, both on the PATH as `plink` and
# `plink2`.
#
# Usage:  bash 02_prepare_genotypes.sh <input.vcf> <output_directory>

set -euo pipefail

VCF=${1:?usage: 02_prepare_genotypes.sh <input.vcf> <output_directory>}
OUT=${2:?usage: 02_prepare_genotypes.sh <input.vcf> <output_directory>}

# Resolve the VCF before changing directory, so a relative path still works.
VCF=$(readlink -f "$VCF")

MAF_COMMON=0.01          # minor allele frequency defining a common variant
N_PER_POPULATION=1000    # diploid individuals per population

mkdir -p "$OUT"
cd "$OUT"

# --- 1. VCF to PLINK binary -------------------------------------------------
# The VCF holds all three populations in one file, so positions are aligned
# across them by construction.
plink --vcf "$VCF" --make-bed --out all_populations

# --- 2. Sample lists, one per population ------------------------------------
# 01_simulate_genotypes.py draws African samples first, then European, then
# East Asian, and msprime names them msp_0, msp_1, and so on in that order. The
# lists below therefore split the individuals into consecutive blocks in that
# same order. The first two columns of the .fam file are the family and
# individual identifiers, which is what --keep expects.
head -n "$N_PER_POPULATION" all_populations.fam | awk '{print $1, $2}' > afr.samples
sed -n "$((N_PER_POPULATION + 1)),$((2 * N_PER_POPULATION))p" \
    all_populations.fam | awk '{print $1, $2}' > eur.samples
sed -n "$((2 * N_PER_POPULATION + 1)),$((3 * N_PER_POPULATION))p" \
    all_populations.fam | awk '{print $1, $2}' > asn.samples

# --- 3. The common-variant panel --------------------------------------------
# Filtering on the pooled sample, so a variant is common if it is common across
# the three populations taken together.
plink2 --bfile all_populations \
       --maf "$MAF_COMMON" \
       --make-bed \
       --out all_populations_common

# --- 4. Split each panel by population --------------------------------------
for POP in afr eur asn; do
    plink --bfile all_populations \
          --keep "${POP}.samples" \
          --keep-allele-order \
          --make-bed \
          --out "${POP}_all_variants"

    plink --bfile all_populations_common \
          --keep "${POP}.samples" \
          --keep-allele-order \
          --make-bed \
          --out "${POP}_common"
done

echo
echo "Genotype preparation complete. Variant counts:"
for STEM in all_populations all_populations_common afr_all_variants afr_common; do
    printf '  %-24s %s variants\n' "$STEM" "$(wc -l < "${STEM}.bim")"
done
