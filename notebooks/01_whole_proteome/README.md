# Whole proteome analysis

## Data processing

Whole proteome data consists of six biological replicates.

- 041521_MOBN_20_WP_2024

- 051121_BN_WP_rep3_2024

- 070521_EV1-46B_10plex_2024

- 091521_EV1-55B_WP_2024

- 20201024_ev_vi_expt7_WP_2024

- 20240408_EV2_68B_2024 

Each experiment was processed individually using following filters:

- Peptides with more than 2 internal missed cleavage sites were removed.

- Half-tryptic peptides were removed.

- Contaminant peptides were removed

- Peptides with low average of reporter ion intensities (<5,000 for all technical replicates) were removed ("floating control")

- Peptides with high variation between all technical replicate channels (CV >0.5) were removed ("floating control")

- Proteins were required to have at least two unique quantified peptides to pass into the final list

- Proteins were required to be quantified in at least two unique biological replicates for analysis

- Proteins detected in two replicates with high variability (29 proteins where the protein was not changing (FC < 1.5) in one replicate and the ratio of the two replicates was greater than two) were also removed ("two replicate variability filter")

## Principal component analysis

Prior to PCA, protein ratio values were log2 transformed and only proteins quantified in all experiments were used

PCA was performed with scikit-learn version 1.4.2 for Python version 3.1.1.

## Volcano plots

p-values were calculated with T-test for the means of two independent samples and volcano plots show uncorrected −log10 p-values on y axis and fold change on x axis. Significant peptides (< 0.05, -log10(0.05) = 1.3) which show a > 1.5 fold change are highlighted.

**Redox-related**

Multiple reference lists of genes of interest were created by searching data from the Gene Ontology Consortium. 

The basic gene ontology, [go-basic.obo](https://geneontology.org/docs/download-ontology/) downloaded on June 20, 2023, was queried using [GOATOOLS v1.2.3](https://github.com/tanghaibao/goatools).

Titles, definitions, and synonyms of GO terms were searched for substrings. The reference gene lists are comprised of genes mapping to GO-terms or children terms that matched the query.

**01_go_terms_matching_query.csv** contains the gene ontology terms and associated genes that matched the query.

**02_genes_matching_query.csv** contains a list of unique genes matching to the query.


## Whole proteome vs RNA correlation analysis

R represents Pearson’s R. Mean was used to combine replicates before calculating fold change. A fold change cutoff of 1.5 was used to color proteins by quadrant. RNA differential expression analysis was performed by Jahan Rahman in DESeq2. Fold change values were not filtered by significance.
