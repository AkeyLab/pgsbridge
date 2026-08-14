# Coalescent genotypes

Three steps produce the genotypes that every trait simulation reads, and check
that they carry realistic between-population structure.

```
sbatch scripts/coalescent/01_simulate_genotypes.slurm       writes sim.vcf
bash   scripts/coalescent/02_prepare_genotypes.sh sim.vcf genotypes/
python scripts/coalescent/03_compute_fst.py --dir genotypes \
       --stems afr_all_variants eur_all_variants asn_all_variants \
       --labels AFR EUR ASN
```

## 1. The demographic model

`01_simulate_genotypes.py` simulates 150 Mb of sequence under the out-of-Africa
history of Fu et al. (2014) and Chen et al. (2020), drawing 1,000 diploid
genomes from each of three populations. All three are written to one VCF, so
positions are aligned across them by construction. Sample order is African,
then European, then East Asian.

Population sizes, in diploid individuals:

| Population | Size |
| --- | --- |
| ancestral hominin | 7,310 |
| Neanderthal lineage | 1,000 |
| African, after its expansion | 15,388 |
| out-of-Africa founding population | 2,758 |
| European, after the split | 1,620 |
| East Asian, after the split | 821 |
| African, present day | 424,000 |
| European, present day | 512,000 |
| East Asian, present day | 1,370,990 |

Event times, in years before the present, converted to generations at 25 years
per generation:

| Event | Years |
| --- | --- |
| human and Neanderthal split | 700,000 |
| African expansion | 316,000 |
| Neanderthal introgression pulse | 55,000 |
| out of Africa | 51,000 |
| European and East Asian split | 28,000 |
| onset of recent growth | 5,115 |

Migration is per generation: 20e-5 between Africa and the out-of-Africa
population, then 1.7e-5 African with European, 0.58e-5 African with East Asian
and 5.9e-5 European with East Asian once the three exist. The per-base
recombination rate is 1e-8 and the mutation rate 1.25e-8.

Note what the ordering of the last two event times implies. The introgression
pulse is dated older than the out-of-Africa merge, so going backwards in time
the branch it targets has already merged into the African population before the
pulse fires. No Neanderthal ancestry reaches the sampled genomes, and no
Neanderthal samples are drawn.

The script requires `msprime` 0.7 and its 0.7-era interface; it does not run
under `msprime` 1.x without changes.

## 2. Splitting the variants

`02_prepare_genotypes.sh` converts the VCF to PLINK format and writes, for each
population, a panel of every segregating variant and a panel of the common ones.

Common variants are those with minor allele frequency at least 0.01 in the three
populations **pooled**, not in each separately. The rare pool is formed later,
inside the trait simulations, as the complement of the common panel on the same
pooled sample, so a rare variant may be monomorphic in one or more populations.

Every per-population file is extracted with `--keep-allele-order`, so a variant's
reference allele is the same in all three. Without it PLINK would reset the coded
allele to whichever is minor within each population and the frequencies would not
be comparable across them.

## 3. Checking the differentiation

`03_compute_fst.py` computes Hudson's fixation index pairwise between the three
populations, as the ratio of sums over variants rather than the mean of
per-variant ratios. Over common variants it gives roughly 0.14 for African with
European, 0.18 for African with East Asian and 0.12 for European with East
Asian, which is the range observed among human continental groups.

## Output

The trait simulations read these files through `load_variant_pools()`, which
expects them under the directory named by `PRS_GENOTYPES` in
`scripts/config.R`:

```
{afr,eur,asn}_all_variants.{bed,bim,fam}
{afr,eur,asn}_common.{bed,bim,fam}
```
