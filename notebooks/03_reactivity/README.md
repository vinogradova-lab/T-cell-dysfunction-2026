# Reactivity analysis

## Data processing: 

Five biological replicates were prepared:

- EVRU1_32

- EV1-46A_10_plex

- EV15C

- HK1-117

- EV2_68A

Each experiment was processed individually using following filters:

- Peptides with more than 2 internal missed cleavage sites were removed.

- Half-tryptic peptides were removed.

- Peptides with low average of reporter ion intensities (<5,000 for all technical replicates) were removed ("floating control")

- Peptides with high variation between all technical replicate channels (CV >0.5) were removed ("floating control")

- Peptides were required to be quantified in at least two unique biological replicates for analysis


## Reactivity change analysis

Three filters were applied to detect reactivity changes:

**bio_replicate_expression_variation_filter**

- To account for potential donor variations in protein expression level, proteins were required to have at least one peptide R ratio within 1.5-fold of the protein expression level measured in TMT-exp experiments (if available) and were excluded from the analysis if all peptide R ratios were greater than 2.0 or less than 0.5.


**expression_filter** and **median_or_ratio filter**

- For proteins with two or more quantified peptides, a cysteine was considered for potential change in reactivity if its peptide R value differed more than two-fold from the protein expression level measured by TMT-exp data, with an additional
requirement that the maximum peptide R ratio differed more than 2-fold from the minimum peptide R ratio. 

- Proteins with one peptide were not considered.