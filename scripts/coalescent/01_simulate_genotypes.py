#!/usr/bin/env python
"""Simulate African, European and East Asian genotypes under a coalescent model
of human demographic history.

The demographic model is the out-of-Africa history of Fu et al. (2014) and Chen
et al. (2020): an ancestral human population that expands in Africa, a founding
bottleneck as one lineage leaves Africa, a later split of that lineage into
European and East Asian branches, recent exponential growth in all three, and
ongoing migration between them.

A Neanderthal lineage is carried in the model, together with a pulse of
introgression into the out-of-Africa branch. No Neanderthal samples are drawn,
and the pulse transfers no ancestry to the sampled genomes: it is dated older
than the out-of-Africa merge, so the branch it targets has already merged into
the African population by the time it fires. The three sampled populations are
therefore unadmixed with respect to the Neanderthal lineage.

The output is a VCF holding all three populations in one file, with positions
aligned across them, which is what the downstream allele-frequency comparisons
require. Sample ordering in the VCF is African, then European, then East Asian.

Run 02_prepare_genotypes.sh next to convert this VCF into the PLINK files that
the trait simulations read.

Requires msprime 0.7 and its 0.7-era interface (PopulationConfiguration,
MassMigration, simulate). It does not run under msprime 1.x without changes.

Example
-------
    python 01_simulate_genotypes.py --out sim.vcf --seed 20220926
"""

import argparse
import math

import msprime

# Sample sizes, in haploid genomes. Each population contributes 2,000 haploid
# genomes, that is 1,000 diploid individuals.
N_HAPLOIDS_PER_POPULATION = 2000

# Length of sequence to simulate, in base pairs, and the per-base per-generation
# recombination and mutation rates.
SEQUENCE_LENGTH = 15e7
RECOMBINATION_RATE = 1e-8
MUTATION_RATE = 1.25e-8


def build_simulation(n_neanderthal, n_african, n_european, n_east_asian, seed):
    """Return an msprime tree sequence for the out-of-Africa model.

    Population indices follow the order in which the configurations are listed:
    0 Neanderthal, 1 African, 2 European, 3 East Asian. Going backwards in time
    the East Asian population merges into the European one to form the
    out-of-Africa population B, which in turn merges into the African
    population, which finally merges with the Neanderthal lineage.
    """
    # --- effective population sizes, in diploid individuals ----------------
    N_ancestral = 7310          # ancestral hominin population
    N_neanderthal = 1000        # Altai and Vindija lineage
    N_african_early = 15388     # African population after its expansion
    N_bottleneck = 2758         # out-of-Africa founding population B
    N_european_early = 1620     # European population after the B split
    N_east_asian_early = 821    # East Asian population after the B split
    N_african_now = 424000      # present-day sizes, after recent growth
    N_european_now = 512000
    N_east_asian_now = 1370990

    # --- event times, converted from years to generations -------------------
    generation_time = 25
    T_human_neanderthal_split = 700e3 / generation_time
    T_african_expansion = 316e3 / generation_time
    T_out_of_africa = 51e3 / generation_time
    T_introgression_pulse = 55e3 / generation_time
    T_european_east_asian_split = 28e3 / generation_time
    T_recent_growth = 5.115e3 / generation_time

    # --- growth rates -------------------------------------------------------
    # Europeans and East Asians grow at a modest rate between their split and
    # the onset of recent growth, then at a faster rate to their present sizes.
    r_european_early = 0.0027
    r_east_asian_early = 0.0031
    N_european_at_growth = N_european_early * math.exp(
        r_european_early * (T_european_east_asian_split - T_recent_growth))
    N_east_asian_at_growth = N_east_asian_early * math.exp(
        r_east_asian_early * (T_european_east_asian_split - T_recent_growth))

    r_african = math.log(N_african_now / N_african_early) / T_recent_growth
    r_european = math.log(N_european_now / N_european_at_growth) / T_recent_growth
    r_east_asian = math.log(N_east_asian_now / N_east_asian_at_growth) / T_recent_growth

    # --- migration rates, per generation ------------------------------------
    m_african_bottleneck = 20e-5    # between Africa and population B
    m_african_european = 1.7e-5
    m_african_east_asian = 0.58e-5
    m_european_east_asian = 5.9e-5
    m_neanderthal = 0               # no continuous Neanderthal migration

    population_configurations = [
        msprime.PopulationConfiguration(initial_size=N_neanderthal),
        msprime.PopulationConfiguration(initial_size=N_african_now,
                                        growth_rate=r_african),
        msprime.PopulationConfiguration(initial_size=N_european_now,
                                        growth_rate=r_european),
        msprime.PopulationConfiguration(initial_size=N_east_asian_now,
                                        growth_rate=r_east_asian),
    ]

    # Neanderthal samples, if any are requested, are drawn from an ancient time
    # point rather than the present.
    T_neanderthal_sampling = 50300 / generation_time
    samples = (
        [msprime.Sample(population=0, time=T_neanderthal_sampling)] * n_neanderthal
        + [msprime.Sample(population=1, time=0)] * n_african
        + [msprime.Sample(population=2, time=0)] * n_european
        + [msprime.Sample(population=3, time=0)] * n_east_asian
    )

    migration_matrix = [
        [0,                    m_neanderthal,        m_neanderthal,         m_neanderthal],
        [m_neanderthal,        0,                    m_african_european,    m_african_east_asian],
        [m_neanderthal,        m_african_european,   0,                     m_european_east_asian],
        [m_neanderthal,        m_african_east_asian, m_european_east_asian, 0],
    ]

    # Events are given in order of increasing time before the present.
    demographic_events = [
        # End of the recent growth phase: all three populations revert to their
        # earlier sizes and growth rates.
        msprime.PopulationParametersChange(
            time=T_recent_growth, initial_size=N_african_early,
            growth_rate=0, population_id=1),
        msprime.PopulationParametersChange(
            time=T_recent_growth, initial_size=N_european_at_growth,
            growth_rate=r_european_early, population_id=2),
        msprime.PopulationParametersChange(
            time=T_recent_growth, initial_size=N_east_asian_at_growth,
            growth_rate=r_east_asian_early, population_id=3),

        # The European and East Asian split. Going backwards, East Asian merges
        # into European to form the out-of-Africa population B, migration is
        # reset so that only Africa and B exchange migrants, and B takes the
        # bottleneck size.
        msprime.PopulationParametersChange(
            time=T_european_east_asian_split, initial_size=N_european_early,
            growth_rate=0, population_id=2),
        msprime.PopulationParametersChange(
            time=T_european_east_asian_split, initial_size=N_east_asian_early,
            growth_rate=0, population_id=3),
        msprime.MassMigration(
            time=T_european_east_asian_split, source=3, destination=2,
            proportion=1.0),
        msprime.MigrationRateChange(time=T_european_east_asian_split, rate=0),
        msprime.MigrationRateChange(
            time=T_european_east_asian_split, rate=m_african_bottleneck,
            matrix_index=(1, 2)),
        msprime.MigrationRateChange(
            time=T_european_east_asian_split, rate=m_african_bottleneck,
            matrix_index=(2, 1)),
        msprime.PopulationParametersChange(
            time=T_european_east_asian_split, initial_size=N_bottleneck,
            growth_rate=0, population_id=2),

        # A brief pulse of Neanderthal introgression into population B.
        msprime.MigrationRateChange(
            time=T_introgression_pulse, rate=1e-3, matrix_index=(2, 0)),
        msprime.MigrationRateChange(
            time=T_introgression_pulse + 30, rate=0, matrix_index=(2, 0)),

        # The out-of-Africa event: population B merges back into Africa.
        msprime.MassMigration(
            time=T_out_of_africa, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_out_of_africa, rate=0),

        # The African expansion, and finally the split from the Neanderthal
        # lineage.
        msprime.PopulationParametersChange(
            time=T_african_expansion, initial_size=N_ancestral, population_id=1),
        msprime.MassMigration(
            time=T_human_neanderthal_split, source=0, destination=1,
            proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_human_neanderthal_split, initial_size=N_ancestral,
            population_id=1),
    ]

    # msprime requires the event list in non-decreasing time order. The
    # introgression pulse is dated 55,000 years before the present and the
    # out-of-Africa merge 51,000, so the pulse is the older of the two and sorts
    # after it. Note what this ordering means for the model: going backwards in
    # time, population 2 has already merged into population 1 by the time the
    # pulse fires, so no lineage remains in population 2 to receive Neanderthal
    # ancestry, and the sampled genomes carry none. The sort is stable, so
    # events sharing a time keep the order written above.
    demographic_events.sort(key=lambda event: event.time)

    return msprime.simulate(
        Ne=N_ancestral,
        length=SEQUENCE_LENGTH,
        recombination_rate=RECOMBINATION_RATE,
        mutation_rate=MUTATION_RATE,
        samples=samples,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        random_seed=seed,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="sim.vcf",
                        help="path of the VCF to write (default: sim.vcf)")
    parser.add_argument("--seed", type=int, default=None,
                        help="random seed; set it to make the run reproducible")
    args = parser.parse_args()

    simulation = build_simulation(
        n_neanderthal=0,
        n_african=N_HAPLOIDS_PER_POPULATION,
        n_european=N_HAPLOIDS_PER_POPULATION,
        n_east_asian=N_HAPLOIDS_PER_POPULATION,
        seed=args.seed,
    )
    with open(args.out, "w") as vcf_file:
        simulation.write_vcf(vcf_file, 2)


if __name__ == "__main__":
    main()
