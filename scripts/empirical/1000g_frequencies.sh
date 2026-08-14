#!/bin/bash
#
# Alternate-allele counts per superpopulation for 1000 Genomes phase 3, over the
# variants that are common in the pooled African, European and East Asian
# sample.  Writes data/allele_frequencies/thousand_genomes_common_af.RDS, which
# is what the 1000 Genomes series of S3 Fig is drawn from.
#
# The genotypes themselves are not distributed with this repository: they are
# the public phase 3 release,
#
#   ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#   integrated_call_samples_v3.20130502.ALL.panel
#
# available from the 1000 Genomes Project.  Point ONE_KG_DIR at a directory
# holding them.  The shipped RDS means the figure redraws without any of this.
#
# Counts rather than frequencies: counts are the raw quantity, and they compress
# to a third of the size.  Allele numbers are fixed because phase 3 has no
# missing genotypes: 2 x 661 African, 2 x 503 European, 2 x 504 East Asian.
#
# Usage:  ONE_KG_DIR=/path/to/1kg-phase3 bash scripts/empirical/1000g_frequencies.sh
#
# Requires bcftools (tested with 1.19).  Runtime is a few hours on one core and
# well under an hour if the 22 chromosomes are farmed out to a cluster.

set -euo pipefail

ONE_KG_DIR="${ONE_KG_DIR:?set ONE_KG_DIR to the directory holding the phase 3 VCFs}"
OUT_DIR="${OUT_DIR:-$(pwd)/1000g_work}"
MAF_COMMON="${MAF_COMMON:-0.01}"
PANEL="$ONE_KG_DIR/integrated_call_samples_v3.20130502.ALL.panel"

mkdir -p "$OUT_DIR/per_chromosome"
cd "$OUT_DIR"

# --- 1. The three superpopulations -------------------------------------------
awk 'NR > 1 && ($3 == "AFR" || $3 == "EUR" || $3 == "EAS") { print $1 "\t" $3 }' \
    "$PANEL" > groups.tsv
awk '{ print $1 }' groups.tsv > samples.txt

# --- 2. Per-chromosome counts ------------------------------------------------
# Biallelic single-nucleotide variants only, to match what the coalescent
# simulation produces.
for CHROM in $(seq 1 22); do
    VCF="$ONE_KG_DIR/ALL.chr${CHROM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    bcftools view -S samples.txt --force-samples -m2 -M2 -v snps -Ou "$VCF" \
      | bcftools +fill-tags -Ou -- -t AC -S groups.tsv \
      | bcftools query -f '%INFO/AC_AFR\t%INFO/AC_EUR\t%INFO/AC_EAS\n' \
      | gzip -1 > "per_chromosome/chr${CHROM}.tsv.gz"
done

# --- 3. Keep the variants common in the pooled sample ------------------------
# The pooled definition is the one used everywhere else in the paper.  Requiring
# a variant to be common in every population separately selects against
# variation private to any one of them and reverses the ordering of the three
# populations, so the choice matters.
for CHROM in $(seq 1 22); do zcat "per_chromosome/chr${CHROM}.tsv.gz"; done \
  | awk -v AN=3336 -v THRESHOLD="$MAF_COMMON" 'BEGIN { OFS = "\t" }
      { allele_count = $1 + $2 + $3
        frequency = allele_count / AN
        minor = (frequency < 1 - frequency) ? frequency : 1 - frequency
        if (minor >= THRESHOLD) print $1, $2, $3 }' \
  | gzip -1 > common_counts.tsv.gz

# --- 4. The file the figure reads --------------------------------------------
Rscript -e '
  allele_number <- c(AFR = 1322L, EUR = 1006L, EAS = 1008L)
  flat <- scan(gzfile("common_counts.tsv.gz"), what = integer(), quiet = TRUE)
  counts <- matrix(flat, ncol = 3, byrow = TRUE,
                   dimnames = list(NULL, names(allele_number)))
  saveRDS(list(counts = counts, allele_number = allele_number),
          "thousand_genomes_common_af.RDS", compress = "xz")
  cat(sprintf("%s common variants written\n", format(nrow(counts), big.mark = ",")))
'

echo "Now copy $OUT_DIR/thousand_genomes_common_af.RDS into data/allele_frequencies/"
