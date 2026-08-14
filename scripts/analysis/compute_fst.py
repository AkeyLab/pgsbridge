#!/usr/bin/env python3
"""Pairwise fixation index between the three simulated populations.

The fixation index F_ST measures how much of the total genetic variation at a
variant sits between two populations rather than within them. It is reported in
the paper to show that the coalescent simulation produces between-population
differentiation comparable to that observed among human continental groups.

Estimator
---------
Hudson's estimator, in the form recommended by Bhatia, Patterson, Sankararaman
and Price (2013, Genome Research 23:1514). At a variant with reference-allele
frequencies p1 and p2 estimated from a1 and a2 sampled alleles,

    numerator_l   = (p1 - p2)^2 - p1(1 - p1)/(a1 - 1) - p2(1 - p2)/(a2 - 1)
    denominator_l = p1(1 - p2) + p2(1 - p1)

The two correction terms remove the upward bias that finite sampling would
otherwise introduce into the squared frequency difference. The genome-wide
value is the ratio of the sums, sum_l numerator_l / sum_l denominator_l, rather
than the average of the per-variant ratios. Averaging ratios is badly behaved
here because the per-variant denominator is near zero at variants that are
nearly monomorphic in both populations.

Allele coding
-------------
Frequencies must refer to the same allele in every population. The .bim files
are checked against the first population's coded allele, and any variant whose
coded allele differs has its frequency replaced by one minus itself. If the
genotypes were prepared with --keep-allele-order, as 02_prepare_genotypes.sh
does, no variant needs this correction and the check simply passes.

Usage
-----
    python compute_fst.py --dir <genotype directory> \
                          --stems afr_all_variants eur_all_variants asn_all_variants \
                          --labels AFR EUR ASN
"""

import argparse
import os

import numpy as np

# PLINK packs four individuals into each byte, two bits each, least significant
# bits first. The two-bit codes are 0 = homozygous for the coded allele,
# 1 = missing, 2 = heterozygous, 3 = homozygous for the other allele.
_DOSAGE = {0: 2, 1: 0, 2: 1, 3: 0}
_NON_MISSING = {0: 1, 1: 0, 2: 1, 3: 1}


def _lookup_tables(n_individuals):
    """Per-byte sums of dosage and of non-missing count, masked to real people.

    The final byte of each variant is padded when the number of individuals is
    not a multiple of four. The padding bits read as code 0, which would
    otherwise be counted as homozygous individuals, so the last byte gets its
    own tables that ignore the padded slots.
    """
    full_dosage = np.zeros(256, dtype=np.int16)
    full_present = np.zeros(256, dtype=np.int16)
    for byte in range(256):
        codes = [(byte >> (2 * k)) & 3 for k in range(4)]
        full_dosage[byte] = sum(_DOSAGE[c] for c in codes)
        full_present[byte] = sum(_NON_MISSING[c] for c in codes)

    n_in_last_byte = n_individuals % 4
    if n_in_last_byte == 0:
        return full_dosage, full_present, full_dosage, full_present

    last_dosage = np.zeros(256, dtype=np.int16)
    last_present = np.zeros(256, dtype=np.int16)
    for byte in range(256):
        codes = [(byte >> (2 * k)) & 3 for k in range(n_in_last_byte)]
        last_dosage[byte] = sum(_DOSAGE[c] for c in codes)
        last_present[byte] = sum(_NON_MISSING[c] for c in codes)
    return full_dosage, full_present, last_dosage, last_present


def count_individuals(stem):
    with open(stem + ".fam") as handle:
        return sum(1 for _ in handle)


def read_coded_alleles(stem):
    """Return the coded allele of each variant, column five of the .bim file."""
    with open(stem + ".bim") as handle:
        return np.array([line.split()[4] for line in handle])


def read_allele_frequencies(stem, n_individuals):
    """Return the coded-allele frequency and the sampled-allele count per variant."""
    bytes_per_variant = (n_individuals + 3) // 4
    raw = np.fromfile(stem + ".bed", dtype=np.uint8)
    if not (raw[0] == 0x6C and raw[1] == 0x1B and raw[2] == 0x01):
        raise ValueError(f"{stem}.bed is not a variant-major PLINK file")

    body = raw[3:]
    n_variants, remainder = divmod(body.size, bytes_per_variant)
    if remainder:
        raise ValueError(f"{stem}.bed has a length inconsistent with its .fam")
    body = body.reshape(n_variants, bytes_per_variant)

    full_dosage, full_present, last_dosage, last_present = _lookup_tables(n_individuals)

    dosage_sum = np.zeros(n_variants, dtype=np.int64)
    present_count = np.zeros(n_variants, dtype=np.int64)
    chunk = 300_000
    for start in range(0, n_variants, chunk):
        block = body[start:start + chunk]
        dosage_sum[start:start + chunk] = (
            full_dosage[block[:, :-1]].sum(axis=1) + last_dosage[block[:, -1]])
        present_count[start:start + chunk] = (
            full_present[block[:, :-1]].sum(axis=1) + last_present[block[:, -1]])

    frequency = dosage_sum / (2.0 * present_count)
    return frequency, 2 * present_count


def hudson_fst(p1, a1, p2, a2, keep):
    """Genome-wide Hudson F_ST over the variants selected by the mask `keep`."""
    numerator = (p1 - p2) ** 2 - p1 * (1 - p1) / (a1 - 1) - p2 * (1 - p2) / (a2 - 1)
    denominator = p1 * (1 - p2) + p2 * (1 - p1)
    return numerator[keep].sum() / denominator[keep].sum()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", required=True,
                        help="directory holding the PLINK files")
    parser.add_argument("--stems", nargs=3, required=True,
                        help="the three PLINK file stems, without extension")
    parser.add_argument("--labels", nargs=3, default=["AFR", "EUR", "ASN"],
                        help="population labels used in the output")
    parser.add_argument("--maf-common", type=float, default=0.01,
                        help="minor allele frequency defining a common variant")
    args = parser.parse_args()

    stems = [os.path.join(args.dir, s) for s in args.stems]
    labels = args.labels

    reference_alleles = read_coded_alleles(stems[0])
    frequencies, allele_counts = {}, {}
    n_variants = None

    for label, stem in zip(labels, stems):
        n_individuals = count_individuals(stem)
        frequency, n_alleles = read_allele_frequencies(stem, n_individuals)
        if n_variants is None:
            n_variants = frequency.size
        elif frequency.size != n_variants:
            raise ValueError("the three populations have different variant counts")

        needs_flip = read_coded_alleles(stem) != reference_alleles
        frequencies[label] = np.where(needs_flip, 1.0 - frequency, frequency)
        allele_counts[label] = n_alleles
        print(f"{label}: {n_individuals:,} individuals, {frequency.size:,} variants, "
              f"{int(needs_flip.sum()):,} realigned to the reference coding",
              flush=True)

    # Classify variants by their minor allele frequency in the pooled sample,
    # weighting the three populations equally.
    pooled = sum(frequencies[label] for label in labels) / len(labels)
    pooled_maf = np.minimum(pooled, 1 - pooled)

    variant_sets = [
        ("all segregating variants", np.ones(n_variants, dtype=bool)),
        (f"common, pooled minor allele frequency at least {args.maf_common}",
         pooled_maf >= args.maf_common),
        (f"rare, pooled minor allele frequency in (0, {args.maf_common})",
         (pooled_maf > 0) & (pooled_maf < args.maf_common)),
    ]

    pairs = [(labels[0], labels[1]), (labels[0], labels[2]), (labels[1], labels[2])]

    for description, mask in variant_sets:
        print(f"\nHudson F_ST over {description}  ({int(mask.sum()):,} variants)")
        for a, b in pairs:
            # Drop variants that are monomorphic for the same allele in both
            # populations: they contribute zero to both sums and carry no
            # information about differentiation.
            both_monomorphic = ((frequencies[a] == frequencies[b]) &
                                ((frequencies[a] == 0) | (frequencies[a] == 1)))
            keep = mask & ~both_monomorphic
            value = hudson_fst(frequencies[a], allele_counts[a],
                               frequencies[b], allele_counts[b], keep)
            print(f"  {a} versus {b}: {value:.4f}")


if __name__ == "__main__":
    main()
