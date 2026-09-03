# Detergent soluble / insoluble proteomics

Unenriched TMT-16plex proteomics of detergent-soluble and detergent-insoluble fractions (run
`HK3-124`, acquired 2026-08-10), profiling protein aggregation across three T cell states:
D2, D8A, D8C x soluble, insoluble x 2 replicates = 12 samples, 5330 proteins, plus four empty
`No_*` TMT channels that are discarded.

**This file covers two runs.** Everything up to the horizontal rule describes HK3-124, the
original run. The section after it describes the HK3-129 + HK3-130 revision run, which has
two donors per group instead of two channel duplicates, and is analysed with a two-donor
filter and no fraction normalization.

## Data processing

Input is `03_combfiles_forpca_channelratio_or_rawsignal_wp.csv` (TMT channel ratios, linear
scale) from the standard unenriched proteomics pipeline, with its default filters including a
minimum of two peptides per protein - not the percent-of-control file, which carries two `inf`
rows (H3C1, PLG) and whose per-row rescaling cancels out of both the fold change and the t
statistic anyway. Zeros are TMT non-detections and become `NaN`, so they are omitted from the
statistics rather than treated as measured values (0-19 per channel).

## Normalization

Per the fractionation protocol (slide 5 of `MMM 2026.8.17 Hiro.pptx`), both fractions were
normalized to 0.9 mg/mL before labelling, so **equal protein mass from each fraction was
TMT-labelled** - the channel ratio reports each protein's share *within its own fraction*, not
its abundance per cell. The slide's "Final dilution" factor `F_g` is the recovered
soluble:insoluble mass ratio, from the measured concentration `C` and resuspension volume `V`
(400 uL soluble, 200 uL insoluble). Insoluble channels are divided by it, soluble channels
left untouched, so the fraction contrast becomes a per-cell ratio:

```
F_g    = (C_sol,g * V_sol) / (C_insol,g * V_insol) = (C_sol,g / C_insol,g) * 2
r'     = r_insol / F_g   (insoluble)        r' = r_sol   (soluble)
log2FC = log2(mean r'_insol / mean r'_sol) = log2(mean r_insol / mean r_sol) - log2(F_g)
```

| state | C_sol (mg/mL) | C_insol (mg/mL) | Sol/Insol | F (Final dilution) |
|---|---|---|---|---|
| D2  | 2.224303294 | 1.055160728 | 2.1080232   | 4.216046399 |
| D8A | 2.554519133 | 1.947635969 | 1.311599895 | 2.623199791 |
| D8C | 2.554519133 | 1.501398348 | 1.70142663  | 3.402853261 |

Soluble-fraction state contrasts are therefore unaffected; insoluble state contrasts shift by
`log2(F_control / F_condition)` (D8A/D2 +0.68, D8C/D2 +0.31, D8C/D8A -0.38); and fraction
contrasts centre at log2(I/S) = -4.11 (D2), -2.91 (D8A), -2.97 (D8C) - most proteins are
largely soluble, and D8A/D8C sit ~1.2 log2 units above D2, the aggregation signal. Read those
off-centre contrasts as global shifts, not widespread individual regulation. The `(Adju)`
sheets of the existing `02_combfiles_percctrl_wp_HK.xlsx` and `..._HK_EV.xlsx` do the same
division, verified on `uniprot` across all 5330 rows.

## Principal component analysis

Channel ratios were log2 transformed and only proteins quantified in all 12 channels were used
(5304 of 5330). PCA was performed with scikit-learn, `n_components=3`, no scaling beyond the
mean-centering sklearn applies. PC1 separates fraction, PC2 separates state. It is run on both
tables, and **the un-normalized PCA is the one to show**: the per-fraction constant is a pure
translation on the log2 scale, so it loads onto PC1 and inflates it (86.9% -> 94.79%) without
moving the PC2 scores that separate the states. Note that the two D2 insoluble replicates are
much less concordant than any other pair (PC1 -128.9 vs -88.7).

## Volcano plots

Nine comparisons: the fraction contrast `<state>-Insoluble` vs `<state>-Soluble` for each of
D2/D8A/D8C, plus the state contrasts D8A/D2, D8C/D2 and D8C/D8A run separately within the
soluble and the insoluble fraction. p-values were calculated with a T-test for the means of
two independent samples and volcano plots show uncorrected -log10 p-values on the y axis and
fold change on the x axis. Proteins with p < 0.05 (-log10(0.05) = 1.3) which show a > 1.5 fold
change are highlighted. BH-adjusted p is stored but not used to call significance.

**Caveat: n = 2.** The replicates are TMT channel duplicates within a single MS run, so the
t-test has 2 degrees of freedom and the tight soluble replicates collapse the pooled SD: 96%
of proteins pass p < 0.05 in the D2 fraction contrast and BH removes almost none (5102 ->
5099). **Read these volcanoes on fold change, not on significance.**

## GSEA

Pre-ranked GSEA (gseapy) against MSigDB C5 GO biological process, on each of the three
insoluble/soluble fraction contrasts, following `notebooks/01_whole_proteome`: 10000
permutations, seed 42. Proteins are ranked by log2(Insoluble / Soluble), so **positive NES
means the set is enriched among proteins gaining intensity in the insoluble fraction**; the
notebook asserts this by checking the leading-edge genes of the top positive-NES sets sit
above the overall median log2(I/S) and the top negative-NES sets below it. GSEA is rank-based,
so the global offset of these contrasts does not affect the result. Each plot shows the top 10
sets per direction at FDR < 0.01, so the insoluble-up side is always represented.

All three states enrich the same insoluble-up biology - alternative mRNA splicing via the
spliceosome, snRNA processing, histone H2A acetylation, nucleotide excision repair - while
oxidative phosphorylation and ATP synthesis stay soluble in D8A/D8C, and nucleotide/amino acid
biosynthesis in D2.

## Outputs

`solubility_analysis.ipynb` writes CSVs to `data/solubility/` (untracked) and
`solubility_visualization.Rmd` draws the figures; `solubility.py` holds the loading,
normalization, statistics and GSEA, with the statistical and annotation functions copied
verbatim from `bin/analysis_utils.py` and `low_input.py` but re-implemented locally to avoid
the goatools dependency. Volcanoes are 2.75 x 2.9 in, not the repo's 2 x 4 in: the axis label
carries two `<state>-<fraction>` names and clips at 2 in.

`supp_data/Data S4_20260817.xlsx` reports the normalized data and all nine comparisons:
`S4-1` is the normalized channel ratio matrix, `S4-2` the fraction comparisons, `S4-3`/`S4-4`
the state comparisons within each fraction, each contributing log2_FC, p_value, p_adj and
Regulation columns. Every sheet carries the UniProt FUNCTION comment, queried with
`unipressed` and cached to `uniprot_functions.csv`; contaminant and keratin rows are dropped
as in `low_input.write_percent_control_to_excel`, leaving 5328 proteins.

`volcano_protein_counts.csv` is the protein-count audit, asserted in both layers: `n_tested +
n_dropped_group_too_small + n_dropped_other == n_input` per comparison (0 to 4 proteins drop
where converting zeros to NaN left a group with fewer than two values), and the Rmd checks via
`ggplot_build()` that the points drawn equal `n_tested`.

---

# Revision run: HK3-129 + HK3-130

Two further TMT runs of the same fractionation design, acquired for the revisions and
combined into one channel-ratio table: `HK3-129` (12-plex, acquired 2026-08-22) and
`HK3-130` (16-plex, 2026-08-23). Each file carries the full D2/D8A/D8C x soluble/insoluble
design as TMT channel duplicates, so a state x fraction group has **four channels: two donors
(`d1` = HK3-129, `d2` = HK3-130) x two duplicates**. 5569 proteins before filtering.

Unlike HK3-124, whose two replicates were channel duplicates within one run, **the replicate
unit here is the donor**.

Input is
`/Users/henrysanford/dev/test_data/texh_revisions/detergent_proteomics/input/03_combined_files/1/03_combfiles_forpca_channelratio_or_rawsignal_wp.csv`,
snapshotted to `data/solubility/revision/01_input/` with its sha256 recorded in the notebook.
Everything runs in the `polars` conda env.

## Two upstream problems in the combined output

Both are worked around by reading the channel-level file directly, and neither affects the
analysis below.

1. **The combined-condition files are unusable.** `parameter_dict.json` carries a stray
   newline (`"D\n2-Insoluble"`) and capitalises `Insoluble` against the lowercase channel
   names, so the combining step matched no columns for any insoluble condition. All three
   insoluble columns of `01_combcond_...csv` and `04_combcond_...csv` are empty.
2. **The pipeline's `min4rep` filter is asymmetric.** It dropped exactly the 411 proteins
   absent from HK3-130 but kept the 317 absent from HK3-129, so its 5158-protein output is
   not a symmetric replicate filter. We apply our own instead.

## Two-donor filter

A protein is kept only if it was quantified in **both donors in every one of the six**
state x fraction groups; a donor counts as quantified when at least one of its two TMT
duplicates is non-zero. A protein seen in only one MS file has no replication at all,
whatever its channel count suggests.

**4823 of 5569 proteins pass** (746 dropped). The audit is written to
`two_donor_filter_audit.csv`, and `two_donor_filter_proteins_of_interest.csv` records the
per-group donor counts for the seven proteins the manuscript follows: MAP2K3, MAP2K4, LONP1,
HSPA9, DNAJA3, NDUFA9, VDAC1. **All seven are quantified in all four channels of all
six groups in both runs, so the filter costs none of them**; the notebook asserts this rather
than assuming it.

## Normalization

**None.** No soluble/insoluble protein concentrations were recorded for HK3-129 or HK3-130,
so there is no "Final dilution" factor `F_g` to divide out as there was for HK3-124. The
fraction contrasts therefore keep a large negative global offset (median log2(I/S) is -3.19
to -3.52) - read that as a global shift, not as widespread individual regulation. State
contrasts within one fraction are unaffected either way.

## Principal component analysis

Channel ratios were log2 transformed and only proteins quantified in all 24 channels were
used (4554 of 4823). scikit-learn, `n_components=3`, no scaling beyond sklearn's
mean-centering. PC1 = 93.88%, PC2 = 1.62%, PC3 = 1.31%. Only the un-normalized PCA exists
here, since there is nothing to normalize. The Rmd shapes points by donor (**circle = `d1`,
square = `d2`**) as well as colouring by group, so run structure is visible alongside the
fraction split. The barplot below draws every replicate point as a plain circle, unshaped by
donor.

## Volcano plots

The same nine comparisons as HK3-124: the fraction contrast `<state>-Insoluble` vs
`<state>-Soluble` for each of D2/D8A/D8C, plus the state contrasts D8A/D2, D8C/D2 and D8C/D8A
run separately within each fraction. Independent two-sample t-test on linear channel ratios,
significance at raw p < 0.05 **and** |log2FC| > log2(1.5); BH-adjusted p is stored but not
used to call significance. The seven proteins of interest are drawn on top of the cloud and
labelled.

**Caveat on the degrees of freedom.** The t-test uses all four channels per group (df = 6),
but those four are two donors x two technical TMT duplicates, so it treats four values as
four independent observations when there are really two. p-values are anti-conservative. As
with HK3-124, read these volcanoes primarily on fold change.

## Scatterplots

log2FC against log2FC with a dashed band one log2(1.5) either side of identity, styled after
`plot_concordance` in `notebooks/07_spin_reactivity`. Drawn from the **state vs state**
family - log2(Insoluble/Soluble) of one state against another (D8A/D2, D8C/D2, D8C/D8A).

`solubility.py` also computes the matched **soluble vs insoluble** contrast (the same state
comparison measured within each fraction) purely as a numerical cross-check: a matched pair is
algebraically the same difference,

```
(I_D8C - S_D8C) - (I_D2 - S_D2)  ==  (I_D8C - I_D2) - (S_D8C - S_D2)
```

so the two share a `delta` and therefore an identical off-band protein set - two rotations of
one quantity. `solubility_analysis.ipynb` asserts that identity numerically; neither this
rotation nor the raw state-vs-state log2FC panel is drawn as a figure.

### The state-vs-state family in % insoluble

The soluble and insoluble shares of a state sum to 100% by construction, so
log2(Insoluble/Soluble) is only a reparameterisation of "% insoluble" - the state-vs-state
scatter and the barplots below report one quantity on two scales, and the percent is the
readable half. `run_scatter` therefore writes `percent_insoluble_x` / `percent_insoluble_y`
onto the three state-vs-state CSVs (`percent_insoluble()`, `100 * I / (I + S)` from the group
means), and the Rmd draws that family on the barplot's scale.

**The classification is unchanged.** On percent axes `|delta| > log2(1.5)` is the odds curve

```
y / (100 - y) = 1.5 * x / (100 - x)
```

which `run_scatter` asserts the percent columns reproduce exactly, so no point can sit on the
wrong side of the drawn curve and the off-band protein set, the corner counts and every
column of `scatter_summary.csv` are identical between the two versions - only the axes
differ. `pearson_r_percent` is reported alongside `pearson_r` because r is scale-dependent
and each panel quotes the scale it draws (0.95/0.96/0.97 on percent, against 0.92/0.94/0.96
on log2). The identity line and the band are `geom_line` over a grid rather than
`geom_abline`, which is laid down in transformed space and would land wrongly under
`scale_*_log10`.

The axis is linear 0-100, the barplot's own scale. 77% of the 4823 proteins sit below 20%
insoluble (median 8-10%, tail to 98%), so the cloud is dense in the bottom-left corner; that
density is the finding, and reading it against the barplot axis is the point of the panel.

**A version coloured by cysteine reactivity change.** `percent_insoluble_reactivity/` holds the
same three panels coloured by a property measured in a *different assay*: whether the protein
carries a cysteine reactivity change in `notebooks/03_reactivity`. It is the only figure in this
notebook that reads another assay's output, and it asks whether a protein whose cysteines
changed reactivity also moved between the detergent fractions.

The highlight set is read by `load_reactivity_changes` from
`data/reactivity/reactivity_changes/output/cysteine_reactivity_changes.csv`, so **notebook 03
must have been run first** - that file is a gitignored build product, and the loader says so when
it is missing. The set is the **union across D4A/D4C/D8A/D8C**: 180 proteins, 136 of them
detected in the 4823-protein table. Unioned rather than state-matched because this run has no D4
states to match against, and a per-panel highlight set would change the coloured points from
panel to panel, which is not the comparison being drawn. The join is on the gene symbol, not on
`uniprot` - the reactivity table carries the peptide's accession and this one the protein-level
accession, and they disagree for isoform entries often enough to lose proteins silently.

Counts land in `scatter_summary.csv` as `n_reactivity_change`, `n_reactivity_source` and
`n_reactivity_named`.

Labels come from `reactivity_named`: `REACTIVITY_NAMED_PROTEINS` in `solubility.py` - MAP2K3,
MAP2K4, LONP1, HSPA9, PSMC5, DNAJA3, NDUFA9, VDAC1, eight in all. This is **not** a subset of
the coloured category - some of the eight carry a reactivity change and some do not - so these
labels are drawn in black rather than in the category colour, which would leave the rest grey
on grey. The reactivity swarm takes one categorical viridis sample, `#440154`.

**What the three panels say: a negative result.** The reactivity-change group does not shift in
a direction of its own in any panel, and is not enriched for off-band proteins either.

| panel | reactivity group | background | off-band, group vs background | Fisher |
|---|---|---|---|---|
| D2 -> D8A | +2.21 pp | +1.84 pp | 14.0% vs 11.8% | OR 1.22, p = 0.42 |
| D2 -> D8C | -0.25 pp | -0.40 pp | 7.4% vs 6.5% | OR 1.14, p = 0.72 |
| D8A -> D8C | -1.99 pp | -2.07 pp | 5.9% vs 7.9% | OR 0.73, p = 0.52 |

Median percentage-point shifts are within a quarter of a point of the background in all three,
and a Mann-Whitney on the shifts is far from significant everywhere (p = 0.40, 0.69, 0.84). Nor
is the group enriched for proteins that leave the band *at all, in either direction* - no panel
comes near significance, and the D8A -> D8C panel points the other way.

**The band width changed this conclusion, so it is worth recording.** At the log2(1.5) band this
family used previously, D2 -> D8A did show a nominal off-band enrichment (46.3% vs 36.2%,
OR 1.52, p = 0.019). It does not survive the move to a 2-fold band (OR 1.22, p = 0.42), which is
what one would expect of an effect carried by the proteins sitting nearest the band edge. The
1.5-fold number should not be quoted.

**An interactive companion.** A 3-inch panel can name eight of its 4823 points, which leaves
128 of the 136 reactivity-change proteins as anonymous dots. `run_reactivity_html` writes a plotly
version of the same three panels into the same folder - one `.html` beside each `.svg`/`.png`,
drawn from the same `joined` frame `run_scatter` returns, so the two cannot disagree about which
protein is in which category. Hovering a point gives the protein and accession, its category,
both percentages, the percentage-point shift, its per-condition reactivity changes and its
UniProt FUNCTION text; clicking a legend entry isolates a category, which is the quickest way to
look at the 136 alone. The R colours are carried across as their exact `col2rgb` hex equivalents
(`gray88` -> `#E0E0E0`, `gray55` -> `#8C8C8C`), since plotly takes only CSS names.

Three places it deliberately departs from the static panel, all for legibility on screen:

- **Trace order.** px emits one trace per category and plotly paints them in that order, so the
  grey cloud is drawn *first* and the highlights on top of it. `legendrank` puts the legend back
  in the static panel's order, so only the painting order differs.
- **Marker size.** Highlights are 7 px against the cloud's 5. The static panel draws everything
  at one size on purpose - there a highlighted protein is a member of the cloud it is being
  compared against - but this file exists to look things up in rather than to compare two
  populations fairly, and the reactivity trace is findable at 5 px only by clicking it off.
- **Labels.** The same eight names, from the same `reactivity_named` column, but plotly has no
  ggrepel. Most of the named proteins sit in the crowded bottom-left, so a fixed offset stacks
  their labels; instead they are fanned across the right half only (the cloud runs bottom-left
  to top-right, so rightwards is the open direction), ordered by y so the labels keep their
  points' vertical order and the arrows do not cross.

These are exploration output; the SVG/PNG panels remain the figures of record. They need
`plotly` - present in the `polars` env, which has no tracked spec file - and a network
connection, since plotly.js is loaded from the CDN rather than inlined (each file is already
4.6 MB of embedded annotation text; inlining the library would add ~4 MB more).

The UniProt text is read straight from the `uniprot_functions.csv` cache rather than through
`uniprot_functions()`: that helper re-queries the entire set the moment one accession is missing,
and one always is - the cache is built after `drop_contaminants` and so never holds
`contaminant_INT-STD1`. 4590 of the 4823 carry a function; the remainder have no FUNCTION comment
in UniProt, which the hover renders as an absent line rather than an empty heading.

The main caveat is that the two assays measure different things on different peptides. A
reactivity change is a cysteine-level call from TMT-ABPP; solubility is a protein-level fraction
ratio. A protein can change reactivity at one cysteine with no change in where the bulk of it
partitions, and the panel would show exactly what it shows. The union across conditions also
mixes D4 and D8 timepoints into a highlight set drawn on D2/D8 axes, so the highlight is not
timepoint-specific.

**Two routes to the same percent.** A scatter point is the share of the channel means; a
barplot bar is the mean of the per-channel shares. Both are legitimate and they are not
identical, so `check_percent_agreement` reconciles them and writes `percent_agreement.csv`:
median 0.04 pp apart, 21 of 14469 protein x state pairs above 2 pp, worst 6.8 pp (P01579,
D8A). The outliers all sit near 40-50% insoluble, where a protein is split evenly enough
between the fractions for channel scatter to move the two averages apart. The assertion
tolerance is 10 pp, set to catch a structural break rather than this.

The relative-share caveat below applies here too: equal protein *mass* was labelled from each
fraction and this run has no Final dilution factor, so "10% insoluble" is a consistent
relative position, not a literal claim that a tenth of the protein is insoluble. Compare
across states, not against 50.

## Proteins of interest

One panel, one CSV, written by `fraction_shares` from `BARPLOT_PROTEIN_SETS`:

| panel | CSV | proteins |
|---|---|---|
| manuscript panel | `fraction_shares_proteins_of_interest.csv` | MAP2K3, MAP2K4, LONP1, HSPA9, PSMC5, DNAJA3, NDUFA9, VDAC1 |

The panel is `PROTEINS_OF_INTEREST` plus PSMC5, one of the 19S regulatory ATPase subunits -
that is the only difference between the two lists. GZMB is intentionally excluded.

`fraction_shares.csv` gives each protein's soluble and insoluble share of its state's total,
one row per replicate channel, so the two bars of a state sum to 100%. Soluble duplicate `i`
is paired with insoluble duplicate `i` of the same donor - the duplicate index carries no
meaning of its own across fractions, but it is the only pairing available, and it keeps every
plotted point a genuine share whose partner sums with it to 100.

Each bar therefore carries **all four channels: two donors x two TMT duplicates**, drawn as
plain circles (`shape 21`) - donor is not shape-coded here, unlike the PCA above. The error
bar is the SE of those four points; two of every four are technical duplicates, so it is an
SE over channels rather than over donors and understates donor-to-donor variability, which
the points do not resolve visually.

The bars are drawn at width 0.6 against a 0.8 dodge, matching the dodged faceted bargraph
in `notebooks/04_metabolomics`; at the shared `BARPLOT_WIDTH` of 0.75 the soluble and
insoluble bars of a state touch. `BARPLOT_WIDTH` itself is unchanged, since
`notebooks/03_reactivity` still uses it. Bars are drawn at `alpha = 0.8`, as in notebooks
01, 02 and 05; the channel points keep full opacity, so they stay legible against their own
bar. The figures are written to `barplots/figures/both_fractions/`, one SVG and PNG per panel.

**These shares are relative, not absolute.** Equal protein *mass* was labelled from each
fraction and no Final dilution factor exists for this run, so the channel ratio reports a
protein's share within its own fraction. A 90% soluble bar does not mean 90% of the protein
is soluble; it means the protein sits that much higher in the soluble fraction's distribution
than in the insoluble fraction's. Compare bars across states, not against 50.

Because percent-of-control is a per-row rescaling of the channel ratio, the constant divides
out of numerator and denominator alike and the share is identical from either table.

## GSEA

Identical to HK3-124: `gseapy.prerank` against MSigDB `c5.go.bp`, 10000 permutations,
seed 42, on each of the three fraction contrasts, ranked by log2(Insoluble/Soluble) so a
positive NES means enrichment among proteins gaining intensity in the insoluble fraction. The
notebook asserts that direction from the leading-edge genes. 139 (D2), 175 (D8A) and 159
(D8C) sets reach FDR < 0.01; the plots show the top 10 per direction.

The result reproduces HK3-124: the same insoluble-up biology - alternative mRNA splicing via
the spliceosome, snRNA processing/metabolism, nucleotide excision repair, spliceosomal complex
assembly - against oxidative phosphorylation, nucleotide metabolism and redox homeostasis
staying soluble.

## Outputs

Written to `data/solubility/revision/` (untracked):

- `03_results/channel_ratios_two_donor_filtered.csv` - the 4823-protein analysis matrix
- `03_results/two_donor_filter_audit.csv`, `..._proteins_of_interest.csv`
- `03_results/pca_unnormalized/`, `volcano_plots/`, `scatterplots/`, `barplots/`, `gsea/`
- `03_results/barplots/percent_agreement.csv` - the barplot-vs-scatter percent reconciliation
- `03_results/barplots/fraction_shares.csv` and one `fraction_shares_<panel>.csv` per entry of
  `BARPLOT_PROTEIN_SETS`

`scatterplots/figures/` holds one folder per kind of panel - `percent_insoluble/` for the
percent version of the state-vs-state family and `percent_insoluble_reactivity/` for the
version coloured by cysteine reactivity change - so each carries one comparable set of three.
The file name is the same in each, since the folder says which panel it is.

`percent_insoluble_reactivity/` carries a `notes.txt` written from the Rmd alongside the
figures, holding the reading above so it survives `data/` being untracked, plus three `.html`
files, the interactive companions written by `run_reactivity_html` from the Python layer.

`supp_data/Data S4b_20260826.xlsx` mirrors `Data S4`: `S4b-1` the filtered channel ratio
matrix, `S4b-2` the fraction comparisons, `S4b-3`/`S4b-4` the state comparisons within each
fraction, every sheet annotated with the UniProt FUNCTION comment. Contaminant and keratin
rows are dropped, leaving 4822 proteins.

## Code

`solubility.py` is parameterised by run through `RunConfig`; `HK3_124` is the default
everywhere, so the original sections are unchanged, and `REVISION` selects this run. The
statistical functions are the ones already copied verbatim from `bin/analysis_utils.py` -
`analysis_utils` itself cannot be imported in the `polars` env because it imports `goatools`
at module level. Replicate columns are always selected by explicit name via `group_columns()`,
never by `filter(like=...)`.

In the Rmd, `plot_volcano()` and `plot_gsea()` take the HK3-124 globals as argument defaults
so their original call sites still work, and the revision section passes its own directories.
