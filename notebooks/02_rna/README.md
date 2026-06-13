# RNA analysis

## Data preparation

Count matrix prepared by Jahan Rahman using

- TrimGalore

- FastQC

- STAR

- bedtools

## GSEA

Bulk RNA-Seq dataset was filtered to genes with a sum of least 1000 reads across 15 samples (14,391 genes).

Genes are ranked based on log2(FC)

### Gene sets

**Chronic Activation Score** The Chronic Activation set is 43 genes selected from the chronic activation list described in [Tooley et al. 2024](https://doi.org/10.1016/j.xcrm.2024.101640)

**Integrated Stress Response Signature** ISR signature genes were pulled from Table S1 of [Chandel 2023](https://doi.org/10.1038/s41586-023-06423-8) This gene list comes from a second study. 99 were identified from Bulk RNA-Seq using Clustering by Inferred Co-expression analysis. 35 additional genes were added through manual curation. Mouse genes were mapped to human orthologs using Ensembl v113. 2 genes were manually mapped using NCBI evidence.

## Heatmaps

A z-score was calculated for each observation by z = (x-μ)/σ where x  is the normalized count, μ is the mean normalized count of the gene across all samples, and σ  is the standard deviation of the normalized count of the gene across all samples. Hierarchical clustering was performed using ComplexHeatmap [Gu et al., 2022](https://doi.org/10.1002/imt2.43) library version 2.10.0 for R version 4.1.1. The distance matrix was calculated using the “euclidean” method and clustering was performed using the “complete” method.