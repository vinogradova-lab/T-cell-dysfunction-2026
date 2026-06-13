# Low input proteomics

### Quality control filters

- Peptides with more than 2 internal missed cleavage sites were removed.

- Half-tryptic peptides were removed.

- Contaminant peptides were removed

- Peptides with low average of reporter ion intensities (<5,000 for all technical replicates) were removed ("floating control")

- Peptides with high variation between all technical replicate channels (CV >0.5) were removed ("floating control")

- Proteins were required to be quantified in at least two unique biological replicates for analysis

- A curated list of proteins (PDCD1, LAG3, TIGIT, HAVCR2 , CD39, NR4A1, NR4A3, GZMB, GSR, CD62L, CTLA4, PRKCH, PRKCQ, and SAMHD1) were included in the final table if they were quantified in a single experiment. For the experiment containing human TILs, the requirement that proteins be quantified in two replicates was not applied. For experiments containing in vivo samples, hemoglobins were excluded.


### Median normalization

Normalization factors were calculated by dividing the global median (the median of all channel medians) by each individual channel's median. Each channel was then multiplied by its corresponding normalization factor.