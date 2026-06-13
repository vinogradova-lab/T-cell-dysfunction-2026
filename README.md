# T cell dysfunction

Code for analysis and visualization of molecular profiling data from the 2026 manuscript on T cell dysfunction.

## Setup

### Python (conda)

To install Python dependencies:

```sh
conda env create -f environment.yml
conda activate t-cell-dysfunction
```

### R (renv)

R dependencies are pinned in `renv.lock` (R 4.1.1). Install them with:

```r
renv::restore()
```

## Layout

The `notebooks/` folder is organized by data type:

- `01_whole_proteome`
- `02_rna` — bulk RNA-Seq
- `03_reactivity`
- `04_metabolomics`
- `05_low_input_proteomics`

Shared helper functions live in `bin/` (Python and R).

Refer to `README.md` in each subfolder of `notebooks/` for additional analysis details.

## Authors

- [Henry Sanford](mailto:hsanford@rockefeller.edu)

- [Nathalie Ropek](mailto:nropek@protonmail.com)