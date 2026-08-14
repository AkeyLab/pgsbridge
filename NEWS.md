# pgstrans 1.0.0

* First release, accompanying the manuscript.
* Closed-form functions for the accuracy ceiling and for intrinsic
  transferability: `prediction_accuracy()`, `transferability()`,
  `additive_genetic_variance()`, `population_mean()`,
  `noise_variance_for_accuracy()`, `shared_effect_size()`,
  `split_effect_sizes()`.
* Trait simulation and measurement on given genotypes: `simulate_trait()`,
  `empirical_accuracy()`, `empirical_transferability()`,
  `sample_causal_variants()`, `allele_frequencies()`, `heritability_split()`.
* Genotype access that never holds a whole variant pool in memory:
  `load_variant_pools()`, `draw_causal_genotypes()`,
  `pool_allele_frequencies()`.
* Summaries and figure builders: `group_summary()`, `label_pair_direction()`,
  `read_transferability_by_heritability()`, `read_rare_variant_sweep()`,
  `plot_direction_pairs()`, `plot_rare_variant_sweep()`,
  `arrange_rare_variant_pair()`.
