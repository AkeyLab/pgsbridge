# Population-genetic and environmental determinants of polygenic score transferability

Code and simulation results for the manuscript.

Kaiqian Zhang, John D. Storey, Joshua M. Akey.
Lewis-Sigler Institute for Integrative Genomics, Princeton University.

## What the paper does

A polygenic score built in one population usually predicts less accurately in
another. This work asks how much of that loss is fixed by population history and
how much better data and methods could still recover. It does so by analysing
the best case: the optimal predictor, in which the causal variants and their
effects are known exactly, are the same in every population, and there is no
linkage disequilibrium among them. The accuracy that predictor attains is a
ceiling no estimated score can exceed.

The closed form for that ceiling in a population `S` is

```
PA_S = 2 sigma_S^2 Phibar_S / (2 sigma_S^2 Phibar_S + tau_S^2),
```

the share of the trait's variance that is genetic, set by three quantities: the
additive genetic variance `sigma_S^2` fixed by the allele-frequency spectrum,
the average self-kinship `Phibar_S` summarising within-population structure, and
the non-genetic variance `tau_S^2`. Transferability from a discovery population
`A` to a target population `B` is the ratio `PA_B / PA_A`.

The code has two parts. The first validates the closed form against forward
simulation and traces its response to each determinant. The second applies it to
African, European and East Asian populations simulated under a coalescent model
of out-of-Africa demographic history.

## Repository layout

```
config.R              every file path used anywhere, in one place
R/                    shared functions, sourced by the scripts
  theory.R            the closed-form results as R functions
  traits.R            simulating traits and measuring accuracy
  variant_pools.R     loading the common and rare variant pools
  summaries.R         group means and confidence intervals
  plots.R             the plot shapes that several figures share
  theme.R             the palette and font size shared by every figure
simulation/           produces the results; numbered in running order
analysis/             the fixation-index calculation
figures/              one script per figure in the paper
data/                 pre-computed results, so figures run without simulating
output/               where figures and new results are written
```

## Reproducing the figures

The `data/` directory holds the results of every simulation, so the figures can
be redrawn without running anything expensive. From the repository root:

```
Rscript figures/figure_01.R
Rscript figures/figure_02.R
Rscript figures/figure_03.R
Rscript figures/figure_04.R
Rscript figures/figure_S1.R
...
Rscript figures/figure_S8.R
```

Each script writes a PDF and a PNG into `output/`. Nothing else needs to be set
up, and no genotype data is required.

## Reproducing the simulations

Running the simulations from scratch takes a cluster and several days. The
scripts in `simulation/` are numbered in the order they must run, and each has a
matching `.slurm` file for submission.

```
sbatch simulation/01_simulate_genotypes.slurm       coalescent genotypes, written as a VCF
bash   simulation/02_prepare_genotypes.sh sim.vcf genotypes/   VCF to PLINK, split by population
sbatch simulation/03_theory_validation.slurm        Figure 1a
sbatch simulation/04_traits_common_variants.slurm   Figure 2, S2 Fig
sbatch simulation/05_traits_normal_effects.slurm    S3 Fig
sbatch simulation/06_traits_rare_and_common.slurm   Figure 3, S4 to S6 Figs
sbatch simulation/07_heritability_decomposition.slurm   Figure 4
sbatch simulation/08_allele_frequencies.slurm       S1, S7, S8 Figs
```

Step 3 is self-contained: it simulates its own genotypes with the `bnpsd`
package and needs no input. Steps 4 to 8 read the PLINK files written by step 2,
whose location is set by `PRS_GENOTYPES` in `config.R`.

The fixation index quoted in the text is computed separately:

```
python analysis/compute_fst.py --dir genotypes \
    --stems afr_all_variants eur_all_variants asn_all_variants \
    --labels AFR EUR ASN
```

## Which script makes which figure

| Figure | Figure script | Simulation that produced its input |
| --- | --- | --- |
| Fig 1 | `figures/figure_01.R` | `simulation/03_theory_validation.R` for panel a; panels b to d are evaluated from the closed form with no simulation |
| Fig 2 | `figures/figure_02.R` | `simulation/04_traits_common_variants.R` |
| Fig 3 | `figures/figure_03.R` | `simulation/06_traits_rare_and_common.R` |
| Fig 4 | `figures/figure_04.R` | `simulation/07_heritability_decomposition.R` |
| S1 Fig | `figures/figure_S1.R` | `simulation/08_allele_frequencies.R` |
| S2 Fig | `figures/figure_S2.R` | `simulation/04_traits_common_variants.R` |
| S3 Fig | `figures/figure_S3.R` | `simulation/04_traits_common_variants.R` and `simulation/05_traits_normal_effects.R` |
| S4 Fig | `figures/figure_S4.R` | `simulation/06_traits_rare_and_common.R` |
| S5 Fig | `figures/figure_S5.R` | `simulation/06_traits_rare_and_common.R` |
| S6 Fig | `figures/figure_S6.R` | `simulation/06_traits_rare_and_common.R` |
| S7 Fig | `figures/figure_S7.R` | `simulation/08_allele_frequencies.R` |
| S8 Fig | `figures/figure_S8.R` | `simulation/08_allele_frequencies.R` |

## What is in `data/`

| Directory | Contents |
| --- | --- |
| `theory_validation/` | the simulated population pairs behind Figure 1a, with their theoretical and empirical transferability |
| `transferability_common_variants/` | transferability replicates for 20, 50, 100 and 1,000 common causal variants, at three heritabilities |
| `transferability_normal_effects/` | the same experiment at 1,000 common causal variants only, with normally distributed effect sizes, at three heritabilities |
| `transferability_rare_variants/` | transferability replicates across nine rare-variant fractions and four effect-size ratios, one directory per heritability |
| `heritability_decomposition/` | the common and rare heritability shares behind Figure 4 |
| `allele_frequencies/` | per-population allele frequencies for the common and rare variant pools |

Every transferability file is a data frame with a column `emp_tau`, one
transferability value per replicate, and a column `pop_comb` naming the ordered
population pair. The label `ASN-AFR` means a discovery population of East Asian
and a target population of African.

## Software

R 4.5.1, with `bnpsd`, `BEDMatrix`, `genio`, `ggplot2`, `ggh4x`,
`patchwork`, `dplyr` and `scales`.

Python 3 with `numpy` for the fixation index. The coalescent simulator needs
`msprime` 0.7 and does not run unchanged under `msprime` 1.x.

PLINK 1.90 and PLINK 2.00 for preparing the genotypes.

If R packages are installed outside the default library path, set `PRS_RLIB` to
that path before running anything.

## Data availability

The theory simulation is self-contained. The coalescent genotypes are produced
by `simulation/01_simulate_genotypes.py`; they are large and are not
redistributed here. All simulation results needed to redraw the figures are
included in `data/`.

## Licence

MIT. See `LICENSE`.
