# pgsbridge

The `pgsbridge` ('polygenic score bridge') R package computes the best
accuracy a polygenic score can attain in a population, and how much of that
accuracy carries over to another population. More specifically, `pgsbridge`
evaluates the closed form for that ceiling, simulates traits on real or
simulated genotypes to measure it directly, and draws the figures of the
accompanying manuscript.

A polygenic score built in one population usually predicts less accurately in
another. Some of that loss is fixed by population history and some could still
be recovered by better data and methods. Separating the two needs a reference
point, and the natural one is the best case: the optimal predictor, in which the
causal variants and their effects are known exactly, are identical across
populations, and carry no linkage disequilibrium. The accuracy it attains is the
share of a population's trait variance that is genetic,

```
PA_S = 2 sigma_S^2 Phibar_S / (2 sigma_S^2 Phibar_S + tau_S^2)
```

set by the additive genetic variance `sigma_S^2` fixed by the allele-frequency
spectrum, the average self-kinship `Phibar_S` summarising within-population
structure, and the non-genetic variance `tau_S^2`. Intrinsic transferability
from a discovery population `A` to a target population `B` is the ratio of their
ceilings, `PA_B / PA_A`.

## Installation

Install the development version from GitHub:

```R
install.packages("remotes")
remotes::install_github("KaiqianZhang/prs-transferability")
```

The genotype readers are optional and only the simulations need them:

```R
install.packages(c("BEDMatrix", "genio"))
```

## Example

The ceiling in a single population, and how it responds to each determinant:

```R
library(pgsbridge)

# a population in Hardy-Weinberg equilibrium whose trait is half genetic
prediction_accuracy(sigma2 = 1, phi_bar = 0.5, tau2 = 1)
#> [1] 0.5

# within-population structure raises the genetic variance, and the ceiling
prediction_accuracy(sigma2 = 1, phi_bar = 0.6, tau2 = 1)
#> [1] 0.5454545

# environmental variance lowers it
prediction_accuracy(sigma2 = 1, phi_bar = 0.5, tau2 = 2)
#> [1] 0.3333333
```

Transferability compares two such ceilings. It exceeds one when the target
population is the more predictable of the two:

```R
# the target has the larger additive genetic variance
transferability(
    sigma2_A = 1.0, phi_bar_A = 0.5, tau2_A = 1,
    sigma2_B = 1.2, phi_bar_B = 0.5, tau2_B = 1
)
#> [1] 1.090909

# identical populations transfer perfectly
transferability(1, 0.5, 1, 1, 0.5, 1)
#> [1] 1
```

## Example with simulated traits

The same quantity can be measured rather than evaluated. Choose effect sizes
that give a target heritability, simulate a trait, and compare the optimal
predictor's measured accuracy against the closed form:

```R
library(pgsbridge)
set.seed(1)

# 1,000 causal variants and 5,000 individuals in Hardy-Weinberg equilibrium
m_loci <- 1000
n_ind <- 5000
p <- runif(m_loci, 0.05, 0.95)
X <- matrix(rbinom(n_ind * m_loci, 2, rep(p, each = n_ind)), n_ind, m_loci)

# one effect size shared by every causal variant, giving 40 percent heritability
beta <- rep(shared_effect_size(h2 = 0.4, tau2 = 1, p = p), m_loci)

# simulate the trait and score the optimal predictor
trait <- simulate_trait(X, beta, alpha = 0, tau2 = 1)
empirical_accuracy(trait$y, trait$g)
#> [1] 0.4034

# which matches the closed form at the realised allele frequencies
prediction_accuracy(
    sigma2 = additive_genetic_variance(beta, allele_frequencies(X)),
    phi_bar = 0.5, tau2 = 1
)
#> [1] 0.4013
```

## Reproducing the manuscript

The repository carries the analysis scripts and the pre-computed simulation
results alongside the package. Clone it, then from the repository root:

```
Rscript scripts/figures/figure_01.R      ...      Rscript scripts/figures/figure_S8.R
```

Each script writes a PDF and a PNG into `output/`. No genotypes and no cluster
are needed; the replicate-level results the figures summarise are in `data/`.

Regenerating those results from scratch needs a cluster and several days. The
scripts in `scripts/simulation/` are numbered in the order they run, each with a
matching `.slurm` wrapper, and `vignette("reproducing-the-figures")` walks
through both paths.

## Citations

Zhang K, Storey JD, Akey JM. "Population-genetic and environmental determinants
of polygenic score transferability."

Ochoa A, Storey JD. "Estimating FST and kinship for arbitrary population
structures." PLoS Genet 17(1): e1009241 (2021).
[doi:10.1371/journal.pgen.1009241](https://doi.org/10.1371/journal.pgen.1009241)
