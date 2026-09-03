"""
SOLUBILITY PROTEOMICS (HK3-124, and the HK3-129 + HK3-130 revision run)

Compares detergent-soluble vs detergent-insoluble fractions across T cell states
(D2, D8A, D8C). Two runs, selected via the `cfg` argument every function takes:

  HK3_124  - original single 16-plex run, two TMT duplicates per group. Default everywhere.
  REVISION - combined HK3-129 (12-plex) + HK3-130 (16-plex) run; replicate unit is the
             DONOR (d1/d2), each with two TMT duplicates, so four channels per group.

Statistical functions and PCA are copied verbatim from bin/analysis_utils.py (as in
notebooks/05_low_input_proteomics/low_input.py) so this module avoids the goatools import.

Figures are drawn in solubility_visualization.Rmd from the CSVs this module writes, except
`run_reactivity_html`, which writes an interactive plotly companion for exploration.
"""

from dataclasses import dataclass
from pathlib import Path
import math
import textwrap

import numpy as np
import pandas as pd
import scipy.stats as stat
from sklearn.decomposition import PCA

# --------------------------------------------------------------------------------------
# Source data
# --------------------------------------------------------------------------------------

SOURCE_CSV = Path(
    "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Exhaustion manuscript/"
    "02_Mass-spectrometry data/"
    "18_unenriched proteomics of detergent soluble and insoluble fraction/"
    "HK3-124_WP_two peptides/03_combined_files/HK3-124/"
    "03_combfiles_forpca_channelratio_or_rawsignal_wp.csv"
)

CHANNEL_SUFFIX = "_processed_census-out_20260810_HK3-124_WP"

# Channel-level file, not the combined-condition files (01_*, 04_*) - a parameter_dict.json
# typo leaves all three insoluble columns empty in those.
REVISION_SOURCE_CSV = Path(
    "/Users/henrysanford/dev/test_data/texh_revisions/detergent_proteomics/input/"
    "03_combined_files/1/03_combfiles_forpca_channelratio_or_rawsignal_wp.csv"
)

ID_COLS = ["uniprot", "protein", "description"]
STATES = ["D2", "D8A", "D8C"]
FRACTIONS = ["Soluble", "Insoluble"]
REPLICATES = ["1", "2"]

# The seven proteins the manuscript follows; labelled on the volcano and scatter panels.
# GZMB is intentionally excluded so labelled figures agree with the barplot panel.
PROTEINS_OF_INTEREST = [
    "MAP2K3",
    "MAP2K4",
    "LONP1",
    "HSPA9",
    "DNAJA3",
    "NDUFA9",
    "VDAC1",
]

# The barplot panel: the seven above plus PSMC5, a 19S regulatory ATPase subunit.
BARPLOT_PANEL = [
    "MAP2K3",
    "MAP2K4",
    "LONP1",
    "HSPA9",
    "PSMC5",
    "DNAJA3",
    "NDUFA9",
    "VDAC1",
]

# Positive controls: proteins with well-characterised solubility behaviour, drawn as their
# own panel rather than folded into the manuscript's eight. GZMB is the same granzyme that
# was on PROTEINS_OF_INTEREST before the revision (see comment above); HSP90B1 and PRF1 are
# the ER/secretory-granule anchors.
POSITIVE_CONTROLS = ["GZMB", "HSP90B1", "PRF1"]

# `fraction_shares` writes one CSV per entry, as `fraction_shares_<key>.csv`.
BARPLOT_PROTEIN_SETS = {
    "proteins_of_interest": BARPLOT_PANEL,
    "positive_controls": POSITIVE_CONTROLS,
}


@dataclass(frozen=True)
class RunConfig:
    """Everything that differs between the two runs.

    `source_template` builds the sample column name as it appears in the source CSV, before
    the per-file suffix. The two runs disagree on both the separator and the case of the
    fraction, so the template carries `{fraction}` and `{fraction_lower}`.

    `file_of` maps each replicate label to the suffix of the MS file it came from; for a
    single-file run every replicate maps to the same suffix.
    """

    name: str
    source: Path
    replicates: tuple
    file_of: dict
    file_labels: dict
    source_template: str
    donors: tuple = ()

    def donor_of(self, replicate):
        """The donor a replicate belongs to, or None for a single-donor run."""
        for donor in self.donors:
            if replicate.startswith(f"{donor}_"):
                return donor
        return None


HK3_124 = RunConfig(
    name="HK3-124",
    source=SOURCE_CSV,
    replicates=tuple(REPLICATES),
    file_of={r: CHANNEL_SUFFIX for r in REPLICATES},
    file_labels={CHANNEL_SUFFIX: "HK3-124"},
    source_template="{state}_{fraction}_{rep}",
)

_HK3_129 = "_processed_census-out_20260822_HK3-129_WP"
_HK3_130 = "_processed_census-out_20260823_HK3-130_WP"

REVISION = RunConfig(
    name="HK3-129+HK3-130",
    source=REVISION_SOURCE_CSV,
    replicates=("d1_1", "d1_2", "d2_1", "d2_2"),
    file_of={
        "d1_1": _HK3_129,
        "d1_2": _HK3_129,
        "d2_1": _HK3_130,
        "d2_2": _HK3_130,
    },
    file_labels={_HK3_129: "HK3-129", _HK3_130: "HK3-130"},
    source_template="{state}-{fraction_lower}_{rep}",
    donors=("d1", "d2"),
)

# --------------------------------------------------------------------------------------
# Normalization
# --------------------------------------------------------------------------------------
# Both fractions were labelled after normalizing to equal protein MASS, so the channel ratio
# reports each protein's share within its own fraction, not per-cell abundance. `FINAL_DILUTION`
# (recovered soluble:insoluble mass ratio, "MMM 2026.8.17 Hiro.pptx" slide 5) corrects for this:
# dividing each insoluble channel by it puts both fractions on a common per-cell basis, matching
# the (Adju) sheets of the existing Excel workbooks.
FINAL_DILUTION = {
    "D2": 4.216046399,
    "D8A": 2.623199791,
    "D8C": 3.402853261,
}

# Measured protein concentrations (mg/mL) behind the factors above, kept for provenance.
PROTEIN_CONC = {
    "D2": {"Soluble": 2.224303294, "Insoluble": 1.055160728},
    "D8A": {"Soluble": 2.554519133, "Insoluble": 1.947635969},
    "D8C": {"Soluble": 2.554519133, "Insoluble": 1.501398348},
}
FRACTION_VOLUME_UL = {"Soluble": 400.0, "Insoluble": 200.0}

# --------------------------------------------------------------------------------------
# Comparisons
# --------------------------------------------------------------------------------------
# (condition, control, contrast_type). Fold change is log2(condition / control).
COMPARISONS = (
    # Fraction contrasts: log2(Insoluble / Soluble) within each state.
    [(f"{s}-Insoluble", f"{s}-Soluble", "fraction") for s in STATES]
    # State contrasts within each fraction.
    + [
        (f"D8A-{f}", f"D2-{f}", f"state_{f.lower()}")
        for f in FRACTIONS
    ]
    + [
        (f"D8C-{f}", f"D2-{f}", f"state_{f.lower()}")
        for f in FRACTIONS
    ]
    + [
        (f"D8C-{f}", f"D8A-{f}", f"state_{f.lower()}")
        for f in FRACTIONS
    ]
)

# Volcanoes sharing an axis range. Fraction contrasts sit far from zero, state contrasts
# near zero, so they are scaled as two families in the Rmd.
FAMILY = {"fraction": "fraction", "state_soluble": "state", "state_insoluble": "state"}


def sample_columns(cfg=HK3_124):
    """The sample columns, in a stable order (12 for HK3-124, 24 for the revision run)."""
    return [
        f"{state}-{fraction}_{rep}"
        for state in STATES
        for fraction in FRACTIONS
        for rep in cfg.replicates
    ]


def group_columns(group, cfg=HK3_124):
    """The replicate columns belonging to a `<state>-<fraction>` group.

    Explicit names, never `filter(like=...)` - `"D2-Insoluble"` ends in `"Soluble"`.
    """
    cols = [f"{group}_{rep}" for rep in cfg.replicates]
    assert len(cols) == len(cfg.replicates), f"expected {len(cfg.replicates)} for {group}"
    return cols


def donor_columns(group, donor, cfg=REVISION):
    """The channels of one donor within a `<state>-<fraction>` group."""
    cols = [f"{group}_{r}" for r in cfg.replicates if cfg.donor_of(r) == donor]
    assert cols, f"no channels for donor {donor} in {group}"
    return cols


def channel_map(cfg=HK3_124):
    """Source column name -> `<state>-<fraction>_<replicate>`."""
    return {
        cfg.source_template.format(
            state=state,
            fraction=fraction,
            fraction_lower=fraction.lower(),
            rep=rep,
        )
        + cfg.file_of[rep]: f"{state}-{fraction}_{rep}"
        for state in STATES
        for fraction in FRACTIONS
        for rep in cfg.replicates
    }


# --------------------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------------------


def peptide_columns(cfg=HK3_124):
    """The per-MS-file peptide count columns, named `pep_num` or `pep_num_<file label>`."""
    labels = list(cfg.file_labels.values())
    if len(labels) == 1:
        return ["pep_num"]
    return [f"pep_num_{label}" for label in labels]


def load_channel_ratios(path=None, cfg=HK3_124):
    """Read the combined channel-ratio table and standardise its column names.

    Renamed to `<state>-<fraction>_<replicate>` - the hyphen is deliberate, since
    `get_pca_plot` derives the condition as `channel_name.split("_")[0]`. Zeros (TMT
    non-detections) become NaN so they're omitted from statistics rather than measured.
    """
    path = cfg.source if path is None else path
    df = pd.read_csv(path, index_col=0)

    # The unused No_* TMT slots are empty.
    df = df.drop(columns=[c for c in df.columns if c.startswith("No_")])

    rename = dict(channel_map(cfg))
    pep_cols = peptide_columns(cfg)
    for suffix, name in zip(cfg.file_labels, pep_cols):
        rename[f"pepNum{suffix}"] = name

    missing = set(rename) - set(df.columns)
    assert not missing, f"missing expected source columns: {sorted(missing)}"
    df = df.rename(columns=rename)

    samples = sample_columns(cfg)
    df[samples] = df[samples].replace(0, np.nan)
    return df[ID_COLS + pep_cols + samples]


def normalize_fractions(df, final_dilution=None, cfg=HK3_124):
    """Divide each insoluble channel by its state's Final dilution factor.

    Only applicable to HK3-124: no soluble/insoluble protein concentrations were recorded
    for HK3-129 or HK3-130, so the revision run is analysed un-normalized.
    """
    final_dilution = FINAL_DILUTION if final_dilution is None else final_dilution
    out = df.copy()
    for state in STATES:
        for rep in cfg.replicates:
            col = f"{state}-Insoluble_{rep}"
            out[col] = out[col] / final_dilution[state]
    return out


# --------------------------------------------------------------------------------------
# Replicate filtering
# --------------------------------------------------------------------------------------


def two_donor_filter(df, cfg=REVISION, out_dir=None, proteins_of_interest=None):
    """Keep proteins quantified in BOTH donors in EVERY state x fraction group.

    A donor counts as quantified when at least one of its two TMT duplicates is not NaN -
    the replicate unit is the donor, not the duplicate. Returns (filtered_df, audit_df);
    `audit_df` records, per protein and group, how many donors it was seen in.
    """
    proteins = (
        PROTEINS_OF_INTEREST if proteins_of_interest is None else proteins_of_interest
    )
    groups = [f"{state}-{fraction}" for state in STATES for fraction in FRACTIONS]

    donors_seen = pd.DataFrame(
        {
            group: sum(
                df[donor_columns(group, donor, cfg)].notna().any(axis=1).astype(int)
                for donor in cfg.donors
            )
            for group in groups
        },
        index=df.index,
    )
    keep = (donors_seen == len(cfg.donors)).all(axis=1)

    audit = pd.concat(
        [df[ID_COLS], donors_seen.add_suffix("_donors"), keep.rename("kept")], axis=1
    )

    # Whether the filter costs any protein of interest is reported, not assumed.
    poi_audit = audit[audit["protein"].isin(proteins)]
    dropped_pois = sorted(poi_audit.loc[~poi_audit["kept"], "protein"])
    print(
        f"two-donor filter: {int(keep.sum())} of {len(df)} proteins kept "
        f"({int((~keep).sum())} dropped)"
    )
    if dropped_pois:
        print(f"  WARNING: proteins of interest DROPPED: {dropped_pois}")
    else:
        print(f"  all {len(poi_audit)} proteins of interest retained")
    missing = sorted(set(proteins) - set(audit["protein"]))
    if missing:
        print(f"  WARNING: proteins of interest absent from the source table: {missing}")

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        audit.to_csv(out_dir / "two_donor_filter_audit.csv", index=False)
        poi_audit.to_csv(out_dir / "two_donor_filter_proteins_of_interest.csv", index=False)

    return df[keep].reset_index(drop=True), audit


# --------------------------------------------------------------------------------------
# Statistics - copied verbatim from bin/analysis_utils.py
# --------------------------------------------------------------------------------------


def get_p_value(row, cond_1, cond_2):
    ttest_result = stat.ttest_ind(row[cond_1], row[cond_2], nan_policy="omit")
    return ttest_result[1]


def get_expr(row):
    """define regulation of row based on fold change and p value"""
    p_value = row["-log10_pval"]
    p_value_cutoff = -1 * math.log10(0.05)
    fc_cutoff = math.log2(1.5)
    if (row["log2_FC"] > fc_cutoff) & (p_value > p_value_cutoff):
        return "Significant Up"
    if (row["log2_FC"] > fc_cutoff) & (p_value < p_value_cutoff):
        return "Not Significant Up"
    if (row["log2_FC"] < -fc_cutoff) & (p_value > p_value_cutoff):
        return "Significant Down"
    if (row["log2_FC"] < -fc_cutoff) & (p_value < p_value_cutoff):
        return "Not Significant Down"
    if (
        (row["log2_FC"] > -fc_cutoff)
        & (row["log2_FC"] < fc_cutoff)
        & (p_value > p_value_cutoff)
    ):
        return "Significant but <1.5 FC"
    else:
        return "Not Significant"


# --------------------------------------------------------------------------------------
# Volcano
# --------------------------------------------------------------------------------------


def run_comparison(df, condition, control, contrast_type, cfg=HK3_124):
    """Volcano statistics for one comparison.

    Returns (long_df, counts) where `counts` accounts for every protein in the input:
    n_tested + the drop reasons must sum to n_input.
    """
    cond_cols = group_columns(condition, cfg)
    ctrl_cols = group_columns(control, cfg)

    work = df.set_index(ID_COLS)

    n_input = len(work)
    n_cond = work[cond_cols].notna().sum(axis=1)
    n_ctrl = work[ctrl_cols].notna().sum(axis=1)

    fc = work[cond_cols].mean(axis=1) / work[ctrl_cols].mean(axis=1)
    log2_fc = np.log2(fc)
    p_value = work.apply(get_p_value, axis=1, args=(cond_cols, ctrl_cols))

    # ttest_ind returns exactly 0.0 when both groups have zero variance, which becomes
    # -log10(p) = inf and would survive dropna() and blow up the plot's axis limits.
    n_p_zero = int((p_value == 0).sum())
    assert n_p_zero == 0, (
        f"{condition} vs {control}: {n_p_zero} proteins have p == 0 (zero variance in "
        "both groups); these would plot at infinite y and must be handled explicitly"
    )

    volcano = pd.DataFrame({"p_value": p_value, "log2_FC": log2_fc})
    keep = np.isfinite(volcano["p_value"]) & np.isfinite(volcano["log2_FC"])

    # Account for every dropped protein.
    dropped = ~keep
    n_group_too_small = int((dropped & ((n_cond < 2) | (n_ctrl < 2))).sum())
    n_other = int(dropped.sum()) - n_group_too_small

    volcano = volcano[keep].copy()
    volcano["-log10_pval"] = -1 * np.log10(volcano["p_value"])
    # BH correction.
    volcano["-log10_pval_adj"] = -1 * np.log10(
        stat.false_discovery_control(volcano["p_value"])
    )
    volcano["Regulation"] = volcano.apply(get_expr, axis=1)
    volcano = volcano.reset_index()

    volcano["condition"] = condition
    volcano["control_condition"] = control
    volcano["contrast_type"] = contrast_type
    volcano["family"] = FAMILY[contrast_type]
    # Renamed from the repo's "-log10_pval": a leading hyphen is mangled by R's read.csv
    # into "X.log10_pval", which is a recurring source of breakage in the existing Rmds.
    volcano = volcano.rename(
        columns={
            "-log10_pval": "neg_log10_pval",
            "-log10_pval_adj": "neg_log10_pval_adj",
        }
    )

    counts = {
        "condition": condition,
        "control_condition": control,
        "contrast_type": contrast_type,
        "n_input": n_input,
        "n_tested": int(keep.sum()),
        "n_dropped_group_too_small": n_group_too_small,
        "n_dropped_other": n_other,
        "n_p_zero": n_p_zero,
        "n_significant_up": int((volcano["Regulation"] == "Significant Up").sum()),
        "n_significant_down": int((volcano["Regulation"] == "Significant Down").sum()),
        "median_log2_FC": float(volcano["log2_FC"].median()),
    }
    assert (
        counts["n_tested"]
        + counts["n_dropped_group_too_small"]
        + counts["n_dropped_other"]
        == n_input
    ), f"{condition} vs {control}: protein count does not reconcile"

    return volcano, counts


def run_comparisons(df, out_dir, comparisons=None, cfg=HK3_124):
    """Run every comparison and write the long volcano table plus the count audit."""
    comparisons = COMPARISONS if comparisons is None else comparisons
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    long_frames, count_rows = [], []
    for condition, control, contrast_type in comparisons:
        volcano, counts = run_comparison(df, condition, control, contrast_type, cfg)
        long_frames.append(volcano)
        count_rows.append(counts)

    long_df = pd.concat(long_frames, axis=0, ignore_index=True)
    counts_df = pd.DataFrame(count_rows)

    assert len(long_df) == counts_df["n_tested"].sum(), (
        "long volcano table does not match the per-comparison tested counts"
    )

    long_df.to_csv(out_dir / "long_volcano_df.csv", index=False)
    counts_df.to_csv(out_dir / "volcano_protein_counts.csv", index=False)
    return long_df, counts_df


# --------------------------------------------------------------------------------------
# Scatterplots
# --------------------------------------------------------------------------------------
# log2FC against log2FC, following `concordance_summary`/`plot_concordance` in
# notebooks/07_spin_reactivity: a point on the diagonal behaves the same in both contrasts,
# and the dashed band one BAND_CUTOFF either side of identity marks a fixed fold difference.

# Scatter-band width only, distinct from the project's significance cutoff (log2(1.5), in
# `get_expr` and `FC_CUTOFF`) - deliberately stricter, since `delta` has no p-value of its own.
BAND_CUTOFF = math.log2(2)

# (x_comparison, y_comparison, kind). Each comparison is a (condition, control) pair naming
# a row family of the long volcano table.
SCATTER_PAIRS = (
    # State vs state: does a protein sit further into the insoluble fraction at D8A/D8C
    # than it does at D2?
    [
        (
            ("D2-Insoluble", "D2-Soluble"),
            (f"{state}-Insoluble", f"{state}-Soluble"),
            "state_vs_state",
        )
        for state in ["D8A", "D8C"]
    ]
    + [
        (
            ("D8A-Insoluble", "D8A-Soluble"),
            ("D8C-Insoluble", "D8C-Soluble"),
            "state_vs_state",
        )
    ]
    # Soluble vs insoluble: the same state contrast measured in each fraction.
    + [
        (
            (f"{cond}-Soluble", f"{ctrl}-Soluble"),
            (f"{cond}-Insoluble", f"{ctrl}-Insoluble"),
            "soluble_vs_insoluble",
        )
        for cond, ctrl in [("D8A", "D2"), ("D8C", "D2"), ("D8C", "D8A")]
    ]
)


def _comparison_slice(long_df, comparison, suffix):
    condition, control = comparison
    sub = long_df[
        (long_df["condition"] == condition) & (long_df["control_condition"] == control)
    ]
    assert len(sub), f"no rows for {condition} vs {control}"
    return sub[ID_COLS + ["log2_FC", "p_value", "Regulation"]].rename(
        columns={
            "log2_FC": f"log2_FC_{suffix}",
            "p_value": f"p_value_{suffix}",
            "Regulation": f"Regulation_{suffix}",
        }
    )


def scatter_label(comparison):
    """`D8C-Insoluble/D8C-Soluble` - the axis label for one comparison."""
    condition, control = comparison
    return f"{condition}/{control}"


def scatter_short_label(comparison, kind):
    """The one word that distinguishes an axis from its partner.

    A state-vs-state pair holds the state constant within each axis and varies it between
    them, so the state names the axis; a soluble-vs-insoluble pair does the reverse. Used
    for the legend and the corner counts, where "higher on y" says nothing.
    """
    condition = comparison[0]
    state, fraction = condition.split("-")
    return state if kind == "state_vs_state" else fraction


def _log2_odds(percent):
    """log2 of the insoluble:soluble odds implied by a percent - the inverse of
    `percent_insoluble`."""
    return np.log2(percent / (100 - percent))


def percent_insoluble(log2_fc):
    """The insoluble share of the state total implied by log2(Insoluble / Soluble).

    100 * I / (I + S) taken from the group means, so the scatter reports the same quantity
    as the `fraction_shares` barplot and the two figures share one unit. Because the shares
    of a state sum to 100 by construction, log2(I/S) is only a reparameterisation of this
    number - the percent is the readable half of the pair.
    """
    odds = np.exp2(log2_fc)
    return 100 * odds / (1 + odds)


OTHER_LABEL = "Other"


# --------------------------------------------------------------------------------------
# Reactivity changes - the cross-assay highlight
# --------------------------------------------------------------------------------------
# Whether a protein whose cysteines changed reactivity also changed detergent solubility -
# the highlight set comes from notebook 03, not from anything measured here.

# Written by notebooks/03_reactivity/reactivity_analysis.ipynb. Resolved from this file, not
# the working directory, so the notebook and a bare `python -c` import it identically. Gitignored
# build artefact, not a tracked input - `load_reactivity_changes` says so when it is missing.
REACTIVITY_CHANGES_CSV = (
    Path(__file__).resolve().parents[2]
    / "data/reactivity/reactivity_changes/output/cysteine_reactivity_changes.csv"
)

REACTIVITY_LABEL = "Reactivity change"

# Text labels on the percent-insoluble-reactivity panel. Named explicitly, not derived from
# `PROTEINS_OF_INTEREST`, since the two lists are allowed to diverge.
REACTIVITY_NAMED_PROTEINS = [
    "MAP2K3",
    "MAP2K4",
    "LONP1",
    "HSPA9",
    "PSMC5",
    "DNAJA3",
    "NDUFA9",
    "VDAC1",
]


def load_reactivity_changes(path=None):
    """The proteins carrying a cysteine reactivity change, unioned across conditions.

    Returns a sorted list of gene symbols. The source table is one row per (condition,
    protein) over D4A/D4C/D8A/D8C; the union is taken because the solubility run has no D4
    states to match against and a per-state highlight would change the coloured set from panel
    to panel, which is not the comparison being drawn.

    Joined on the gene symbol, not on `uniprot`: the reactivity table carries the peptide's
    accession and the solubility table the protein-level one, and they disagree for isoform
    entries often enough that a uniprot join silently loses proteins.
    """
    path = REACTIVITY_CHANGES_CSV if path is None else Path(path)
    assert path.exists(), (
        f"{path} not found - run notebooks/03_reactivity/reactivity_analysis.ipynb first "
        "(everything under data/ is gitignored)"
    )
    df = pd.read_csv(path)
    expected = {"condition", "uniprot", "protein", "direction_of_reactivity_change"}
    missing = expected - set(df.columns)
    assert not missing, f"{path} is missing {sorted(missing)}"

    proteins = sorted(df["protein"].dropna().unique())
    counts = df.groupby("condition")["protein"].nunique().to_dict()
    breakdown = "  ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"{len(proteins)} proteins with a reactivity change   {breakdown}")
    return proteins


def reactivity_hover_detail(path=None):
    """One row per protein: the per-condition reactivity changes, as hover text.

    Returns `protein`, `reactivity_detail`, one line per row of the source table (e.g.
    `D8A Higher C91, C93`) joined with `<br>`. One line per SOURCE ROW, not per condition -
    a protein can have several rows per condition with different residues/directions
    (e.g. PHF1/D8A), so collapsing by condition would lose information.
    """
    path = REACTIVITY_CHANGES_CSV if path is None else Path(path)
    assert path.exists(), (
        f"{path} not found - run notebooks/03_reactivity/reactivity_analysis.ipynb first "
        "(everything under data/ is gitignored)"
    )
    df = pd.read_csv(path).sort_values(["protein", "condition"])

    lines = (
        df["condition"].astype(str)
        + " "
        + df["direction_of_reactivity_change"].astype(str)
        + " "
        + df["residues_reactivity_changes"].fillna("").astype(str).str.replace(",", ", ")
    ).str.strip()

    detail = (
        df.assign(line=lines)
        .groupby("protein")["line"]
        .apply(lambda s: "<br>".join(s))
        .reset_index(name="reactivity_detail")
    )

    # Nothing dropped on the way through: every source row is still a line, and every protein
    # still has a row.
    n_lines = detail["reactivity_detail"].str.count("<br>").sum() + len(detail)
    assert n_lines == len(df), f"{n_lines} hover lines from {len(df)} source rows"
    assert len(detail) == df["protein"].nunique()
    return detail


def run_scatter(long_df, x_comparison, y_comparison, kind, out_dir,
                proteins_of_interest=None, reactivity_proteins=None):
    """Join two comparisons on protein and score their agreement.

    Returns (joined_df, summary). `concordance` is "Concordant" inside the 2-fold band and
    otherwise names the axis the protein is higher on, so the R layer can colour and label
    without recomputing anything.
    """
    proteins = (
        PROTEINS_OF_INTEREST if proteins_of_interest is None else proteins_of_interest
    )
    reactivity = (
        load_reactivity_changes() if reactivity_proteins is None else reactivity_proteins
    )
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    x_label, y_label = scatter_label(x_comparison), scatter_label(y_comparison)
    joined = _comparison_slice(long_df, x_comparison, "x").merge(
        _comparison_slice(long_df, y_comparison, "y"), on=ID_COLS, how="inner"
    )

    x_short = scatter_short_label(x_comparison, kind)
    y_short = scatter_short_label(y_comparison, kind)

    joined["delta"] = joined["log2_FC_y"] - joined["log2_FC_x"]
    joined["concordance"] = np.where(
        joined["delta"] > BAND_CUTOFF,
        f"Higher in {y_short}",
        np.where(joined["delta"] < -BAND_CUTOFF, f"Higher in {x_short}", "Concordant"),
    )
    joined["protein_of_interest"] = joined["protein"].isin(proteins)
    joined["x_label"] = x_label
    joined["y_label"] = y_label
    joined["kind"] = kind

    # State-vs-state contrasts sum to 100% of the state by construction, so they're equally
    # readable as "% insoluble here vs there" - the barplot's unit. Carried on the same rows
    # as the log2 values so the R layer can't pair a point with the wrong concordance call.
    if kind == "state_vs_state":
        joined["percent_insoluble_x"] = percent_insoluble(joined["log2_FC_x"])
        joined["percent_insoluble_y"] = percent_insoluble(joined["log2_FC_y"])
        # On percent axes the band is the odds curve y/(100-y) = 2 * x/(100-x).
        # Assert it is exactly `delta`, so a point can never sit on the wrong side of the
        # dashed curve the Rmd draws.
        odds_delta = _log2_odds(joined["percent_insoluble_y"]) - _log2_odds(
            joined["percent_insoluble_x"]
        )
        assert np.allclose(odds_delta, joined["delta"]), (
            f"{x_label} vs {y_label}: percent insoluble does not reproduce delta"
        )

        # The percentage-point shift, used by the reactivity HTML's hover text - the log2 odds
        # `delta` is not itself in percentage points.
        joined["percent_point_delta"] = (
            joined["percent_insoluble_y"] - joined["percent_insoluble_x"]
        )

        # Reactivity changes, on a flag + display-label pattern.
        joined["reactivity_change"] = joined["protein"].isin(reactivity)
        joined["reactivity_label"] = np.where(
            joined["reactivity_change"], REACTIVITY_LABEL, OTHER_LABEL
        )

        # Not a subset of `reactivity_change` - some named proteins lack a reactivity change -
        # so the R layer draws these labels in black rather than the category colour.
        joined["reactivity_named"] = joined["protein"].isin(REACTIVITY_NAMED_PROTEINS)

    summary = {
        "kind": kind,
        "x_label": x_label,
        "y_label": y_label,
        "x_short": x_short,
        "y_short": y_short,
        "n_joined": len(joined),
        "pearson_r": float(joined["log2_FC_x"].corr(joined["log2_FC_y"])),
        "median_delta": float(joined["delta"].median()),
        "fc_cutoff": BAND_CUTOFF,
        "n_concordant": int((joined["concordance"] == "Concordant").sum()),
        "n_higher_on_x": int((joined["concordance"] == f"Higher in {x_short}").sum()),
        "n_higher_on_y": int((joined["concordance"] == f"Higher in {y_short}").sum()),
        "n_proteins_of_interest": int(joined["protein_of_interest"].sum()),
    }
    if kind == "state_vs_state":
        # Reported alongside `pearson_r`, not instead of it: r on the percent scale is a
        # different number from r on the log2 scale, and the figure quotes whichever scale
        # it draws.
        summary["pearson_r_percent"] = float(
            joined["percent_insoluble_x"].corr(joined["percent_insoluble_y"])
        )
        summary["n_reactivity_change"] = int(joined["reactivity_change"].sum())
        # The size of the source set, so the panel can annotate "136 of 180 detected" without
        # the R layer ever having to open the reactivity table.
        summary["n_reactivity_source"] = len(reactivity)
        summary["n_reactivity_named"] = int(joined["reactivity_named"].sum())
    assert (
        summary["n_concordant"] + summary["n_higher_on_x"] + summary["n_higher_on_y"]
        == summary["n_joined"]
    ), f"{x_label} vs {y_label}: concordance classes do not reconcile"

    fn = f"{x_label} vs {y_label}".replace("/", "_over_")
    joined.to_csv(out_dir / f"{fn}.csv", index=False)
    return joined, summary


def run_scatters(long_df, out_dir, pairs=None, proteins_of_interest=None,
                 reactivity_proteins=None):
    """Every scatterplot pair, plus the summary table the R layer reads."""
    pairs = SCATTER_PAIRS if pairs is None else pairs
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Loaded once here rather than per pair, so the reactivity table is read a single time and
    # every panel is guaranteed to be highlighting the same set.
    reactivity = (
        load_reactivity_changes() if reactivity_proteins is None else reactivity_proteins
    )

    rows = []
    for x_comparison, y_comparison, kind in pairs:
        joined, summary = run_scatter(
            long_df, x_comparison, y_comparison, kind, out_dir, proteins_of_interest,
            reactivity_proteins=reactivity,
        )
        rows.append(summary)
        print(
            f"{summary['x_label']:28s} vs {summary['y_label']:28s} "
            f"n={summary['n_joined']:5d}  r={summary['pearson_r']:5.2f}"
        )

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(out_dir / "scatter_summary.csv", index=False)
    return summary_df


# --------------------------------------------------------------------------------------
# Interactive companion to the percent-insoluble reactivity panels
# --------------------------------------------------------------------------------------
# Most points on the static panel are unlabelled; these HTML files let you hover a point for
# its protein, UniProt function, and reactivity status. Exploration aids, not manuscript
# figures - the SVG/PNG panels remain the output of record. Follows bin/analysis_utils.py:
# plotly express, `template="plotly_white"`, plotly.js from the CDN.

# The UniProt FUNCTION cache `write_percent_insoluble_excel` leaves behind - read directly,
# not via `uniprot_functions()`, which re-queries the whole set if one accession is missing.
UNIPROT_FUNCTIONS_CSV = "uniprot_functions.csv"

# Characters per line for hover text; plotly does not wrap hover labels itself.
HOVER_WRAP = 90

# The static panel's colours, as hex. plotly only takes CSS colours, and R's `gray88` /
# `gray55` / `lightgrey` are not CSS names - these are their exact `col2rgb` equivalents, so the
# two versions of the panel are the same picture.
HTML_COLOURS = {
    REACTIVITY_LABEL: "#440154",
    OTHER_LABEL: "#E0E0E0",          # gray88
}
CLOUD_OUTLINE = "#8C8C8C"            # gray55
IDENTITY_COL = "#D3D3D3"             # lightgrey, the Rmd's UNCHANGED_COL

# The legend reads in the static panel's order even though the traces are drawn back to front.
LEGEND_ORDER = {REACTIVITY_LABEL: 1, OTHER_LABEL: 2}

# How far a text label sits from its point, in pixels, and how wide the fan is either side of
# horizontal. Big enough to clear the bottom-left cluster most of the named proteins sit in.
LABEL_RADIUS = 85
LABEL_FAN_DEG = 72


def _percent_band(band_cutoff=None):
    """The fold-change band on percent axes, as two curves.

    The same construction as `band_grid` in the Rmd: on percent axes the band is the odds curve
    y/(100-y) = r * x/(100-x) for the band's fold ratio r. Stops short of 0 and 100, where the
    odds are undefined.
    """
    ratio = 2 ** (BAND_CUTOFF if band_cutoff is None else band_cutoff)
    x = np.linspace(0.05, 99.95, 500)
    odds = x / (100 - x)
    up = 100 * (ratio * odds) / (1 + ratio * odds)
    down = 100 * (odds / ratio) / (1 + odds / ratio)
    return x, up, down


def write_reactivity_html(joined, out_dir, functions=None, reactivity_detail=None):
    """One interactive percent-insoluble panel, coloured by reactivity change.

    Takes the `joined` frame `run_scatter` returns, so the HTML and the static panel are drawn
    from one table and cannot disagree about which proteins are in which category.
    """
    import plotly.express as px

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    kind = joined["kind"].iloc[0]
    assert kind == "state_vs_state", f"{kind} panels carry no percent columns"
    required = ["reactivity_label", "percent_insoluble_x", "percent_insoluble_y",
                "percent_point_delta"]
    missing = [c for c in required if c not in joined.columns]
    assert not missing, f"{missing} missing - regenerate the CSV with run_scatter"

    x_label, y_label = joined["x_label"].iloc[0], joined["y_label"].iloc[0]
    x_short, y_short = x_label.split("-")[0], y_label.split("-")[0]

    df = joined.copy()
    n_before = len(df)
    if functions is not None:
        df = df.merge(functions, on="uniprot", how="left")
    if reactivity_detail is not None:
        df = df.merge(reactivity_detail, on="protein", how="left")
    # A duplicated accession or protein symbol on either side would multiply points silently.
    assert len(df) == n_before, f"annotation join changed the row count {n_before} -> {len(df)}"

    for col in ["uniprot_function", "reactivity_detail"]:
        if col not in df.columns:
            df[col] = ""
    # 232 of the cached functions are legitimately empty strings rather than nulls, so both have
    # to collapse to the same thing.
    df["uniprot_function"] = df["uniprot_function"].fillna("").astype(str)
    df["reactivity_detail"] = df["reactivity_detail"].fillna("").astype(str)
    df["function_wrapped"] = [
        "<br>".join(textwrap.wrap(f, HOVER_WRAP)) if f else "" for f in df["uniprot_function"]
    ]

    fig = px.scatter(
        df,
        x="percent_insoluble_x",
        y="percent_insoluble_y",
        color="reactivity_label",
        color_discrete_map=HTML_COLOURS,
        # Draw order, not legend order - grey cloud first or it buries the highlighted points.
        # `legendrank` below puts the legend back the other way round.
        category_orders={
            "reactivity_label": [OTHER_LABEL, REACTIVITY_LABEL]
        },
        custom_data=["protein", "uniprot", "reactivity_label", "percent_point_delta",
                     "reactivity_detail", "function_wrapped"],
        template="plotly_white",
    )

    # `%{customdata[n]}` on an empty string collapses the line to nothing but leaves the label
    # it was introduced with, so the reactivity and function blocks carry their own headings
    # inside the value rather than as literal text here.
    df["reactivity_detail"] = np.where(
        df["reactivity_detail"] == "", "", "<br>" + df["reactivity_detail"]
    )
    df["function_wrapped"] = np.where(
        df["function_wrapped"] == "", "", "<br>" + df["function_wrapped"]
    )
    for trace in fig.data:
        sub = df[df["reactivity_label"] == trace.name]
        # px keeps row order within a category, but the hover text is only correct if it does -
        # a reordered trace would put every protein's name on someone else's point.
        assert np.allclose(sub["percent_insoluble_x"], trace.x), (
            f"{trace.name}: customdata does not line up with the drawn points"
        )
        trace.customdata = sub[
            ["protein", "uniprot", "reactivity_label", "percent_point_delta",
             "reactivity_detail", "function_wrapped"]
        ].to_numpy()
        trace.hovertemplate = (
            "<b>%{customdata[0]}</b>  (%{customdata[1]})<br>"
            "%{customdata[2]}<br>"
            f"% insoluble, {x_short}: " "%{x:.1f}<br>"
            f"% insoluble, {y_short}: " "%{y:.1f}<br>"
            "shift: %{customdata[3]:+.1f} pp"
            "%{customdata[4]}"
            "%{customdata[5]}"
            "<extra></extra>"
        )
        # The static panel draws every point at one size, because there a highlighted protein
        # is a member of the cloud being compared against. Here the highlighted trace is drawn
        # a little larger for findability, since this file exists to be looked things up in
        # rather than to make a fair visual comparison of the two populations.
        highlighted = trace.name != OTHER_LABEL
        trace.marker.update(
            size=7 if highlighted else 5,
            line=dict(width=0.5, color="black" if highlighted else CLOUD_OUTLINE),
            opacity=1.0 if highlighted else 0.45,
        )
        trace.legendrank = LEGEND_ORDER[trace.name]

    # Identity line and the fold-change band. `hoverinfo="skip"` so a curve never steals the hover
    # from a point underneath it.
    fig.add_shape(type="line", x0=0, y0=0, x1=100, y1=100,
                  line=dict(color=IDENTITY_COL, width=1), layer="below")
    band_x, band_up, band_down = _percent_band()
    for band in (band_up, band_down):
        fig.add_scatter(x=band_x, y=band, mode="lines", hoverinfo="skip",
                        showlegend=False, line=dict(color="black", width=1, dash="dash"))

    # plotly has no ggrepel - labels are fanned across the right half (the cloud's open
    # direction) with the fan angle swept by y-order, so arrows don't cross.
    named = df[df["reactivity_named"]].sort_values("percent_insoluble_y")
    assert len(named) == int(df["reactivity_named"].sum())
    for i, row in enumerate(named.itertuples()):
        frac = i / max(len(named) - 1, 1)
        angle = np.radians(-LABEL_FAN_DEG + 2 * LABEL_FAN_DEG * frac)
        fig.add_annotation(
            x=row.percent_insoluble_x,
            y=row.percent_insoluble_y,
            text=row.protein,
            showarrow=True,
            arrowhead=0,
            arrowwidth=0.7,
            arrowcolor="black",
            ax=LABEL_RADIUS * np.cos(angle),
            ay=-LABEL_RADIUS * np.sin(angle),
            font=dict(size=11, color="black"),
            bgcolor="rgba(255,255,255,0.75)",
            borderpad=1,
        )

    n_changed = int((df["reactivity_label"] == REACTIVITY_LABEL).sum())
    fig.update_layout(
        title=dict(
            text=f"% insoluble, {x_label} vs {y_label}"
                 f"<br><sub>{n_changed} reactivity-change proteins detected;"
                 f" {len(named)} labelled; dashed band = {2 ** BAND_CUTOFF:.0f}-fold;"
                 " click a legend entry to isolate a category</sub>",
            font=dict(size=13),
        ),
        xaxis_title=f"% insoluble, {x_short}",
        yaxis_title=f"% insoluble, {y_short}",
        legend_title_text="",
        width=850,
        height=800,
        hoverlabel=dict(align="left"),
    )
    fig.update_xaxes(range=[0, 100], constrain="domain")
    fig.update_yaxes(range=[0, 100], scaleanchor="x", scaleratio=1)

    drawn = sum(len(t.x) for t in fig.data if t.mode != "lines")
    assert drawn == len(df), f"{drawn} points drawn of {len(df)}"

    fn = f"{x_label} vs {y_label}".replace("/", "_over_")
    path = out_dir / f"{fn}.html"
    fig.write_html(path, include_plotlyjs="cdn")
    return path


def run_reactivity_html(long_df, scatter_dir, functions_csv=None, reactivity_proteins=None):
    """The three interactive percent-insoluble reactivity panels.

    Written into the same `figures/percent_insoluble_reactivity/` folder as the static panels,
    under the same file stem, so each HTML sits beside its own .svg/.png.
    """
    scatter_dir = Path(scatter_dir)
    out_dir = scatter_dir / "figures" / "percent_insoluble_reactivity"

    functions_csv = (
        scatter_dir.parent / UNIPROT_FUNCTIONS_CSV if functions_csv is None
        else Path(functions_csv)
    )
    assert functions_csv.exists(), (
        f"{functions_csv} not found - run write_percent_insoluble_excel first, which builds it"
    )
    functions = pd.read_csv(functions_csv, keep_default_na=False)
    detail = reactivity_hover_detail()
    reactivity = (
        load_reactivity_changes() if reactivity_proteins is None else reactivity_proteins
    )

    paths = []
    for x_comparison, y_comparison, kind in SCATTER_PAIRS:
        if kind != "state_vs_state":
            continue
        joined, _ = run_scatter(
            long_df, x_comparison, y_comparison, kind, scatter_dir,
            reactivity_proteins=reactivity,
        )
        annotated = joined["uniprot"].isin(
            functions.loc[functions["uniprot_function"] != "", "uniprot"]
        ).sum()
        path = write_reactivity_html(joined, out_dir, functions, detail)
        paths.append(path)
        print(f"{path.name}   {len(joined)} proteins, {annotated} with a UniProt function")
    return paths


# --------------------------------------------------------------------------------------
# Fraction shares - the barplot input
# --------------------------------------------------------------------------------------


def fraction_shares(df, cfg=REVISION, out_dir=None, protein_sets=None):
    """Each protein's soluble and insoluble share of its state's total, per channel.

    One share per replicate channel (two donors x two TMT duplicates); soluble duplicate `i`
    pairs with insoluble duplicate `i` of the same donor so every point's partner sums to 100.

    CAVEAT: equal protein MASS was labelled per fraction with no Final dilution factor, so
    shares are a relative display, not a literal "X% of the protein is insoluble".

    `protein_sets` names the barplot panels: one `fraction_shares_<key>.csv` per entry. A panel
    whose genes are not all quantified fails here rather than silently drawing short.
    """
    sets = BARPLOT_PROTEIN_SETS if protein_sets is None else protein_sets
    proteins = sets["proteins_of_interest"]

    records = []
    for state in STATES:
        for replicate in cfg.replicates:
            soluble = df[f"{state}-Soluble_{replicate}"]
            insoluble = df[f"{state}-Insoluble_{replicate}"]
            total = soluble + insoluble
            for fraction, value in [("Soluble", soluble), ("Insoluble", insoluble)]:
                block = df[ID_COLS].copy()
                block["state"] = state
                block["donor"] = cfg.donor_of(replicate)
                block["replicate"] = replicate
                block["duplicate"] = replicate.split("_")[-1]
                block["fraction"] = fraction
                block["group"] = f"{state}-{fraction}"
                block["percent"] = 100 * value / total
                records.append(block)

    shares = pd.concat(records, ignore_index=True)
    shares["protein_of_interest"] = shares["protein"].isin(proteins)

    # Soluble + insoluble must total 100 for every protein x state x channel.
    totals = shares.dropna(subset=["percent"]).groupby(
        ["uniprot", "state", "replicate"], sort=False
    )["percent"].sum()
    assert np.allclose(totals, 100), "fraction shares do not sum to 100"

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        shares.to_csv(out_dir / "fraction_shares.csv", index=False)
        print(f"fraction shares: {shares['uniprot'].nunique()} proteins")
        for name, panel in sets.items():
            subset = shares[shares["protein"].isin(panel)]
            missing = sorted(set(panel) - set(subset["protein"]))
            assert not missing, f"{name}: not quantified in this run: {missing}"
            subset.to_csv(out_dir / f"fraction_shares_{name}.csv", index=False)
            print(f"  {name:22s} {subset['protein'].nunique()} proteins, {len(subset)} rows")

    return shares


def check_percent_agreement(shares, long_df, out_dir=None, tolerance_pp=10.0):
    """Reconcile the barplot's percent with the scatter's percent.

    A barplot bar is the MEAN OF THE PER-CHANNEL SHARES (`fraction_shares`); a scatter point
    is the SHARE OF THE CHANNEL MEANS (`percent_insoluble`). Both legitimate, not identical -
    asserts the worst disagreement stays under `tolerance_pp`, so a real divergence fails here.
    """
    bars = (
        shares[shares["fraction"] == "Insoluble"]
        .groupby(["uniprot", "state"], as_index=False)["percent"]
        .mean()
        .rename(columns={"percent": "percent_mean_of_shares"})
    )

    contrasts = long_df[long_df["contrast_type"] == "fraction"].copy()
    contrasts["state"] = contrasts["condition"].str.split("-").str[0]
    contrasts["percent_share_of_means"] = percent_insoluble(contrasts["log2_FC"])

    merged = bars.merge(
        contrasts[["uniprot", "state", "percent_share_of_means"]],
        on=["uniprot", "state"],
        how="inner",
    ).dropna(subset=["percent_mean_of_shares", "percent_share_of_means"])
    merged["difference_pp"] = (
        merged["percent_share_of_means"] - merged["percent_mean_of_shares"]
    )

    worst = merged["difference_pp"].abs().max()
    print(
        f"percent agreement (barplot vs scatter): n={len(merged)}, "
        f"median |diff| = {merged['difference_pp'].abs().median():.3f} pp, "
        f"max |diff| = {worst:.3f} pp"
    )

    # Written before the assertion so the diagnostic survives a failure.
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        merged.to_csv(out_dir / "percent_agreement.csv", index=False)

    assert worst < tolerance_pp, (
        f"mean-of-shares and share-of-means disagree by {worst:.2f} pp, "
        f"above the {tolerance_pp} pp tolerance"
    )

    return merged


# --------------------------------------------------------------------------------------
# PCA - copied from bin/analysis_utils.py::get_pca_plot, minus the throwaway plotly figure
# --------------------------------------------------------------------------------------


def get_pca_plot(df, index_cols, out_dir):
    """Principal component analysis on log2 channel ratios.

    Writes the three CSVs the R layer expects: pca_results.csv, loadings_results.csv and
    percent_explained.csv. Only proteins quantified in every channel are used.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = df.dropna()
    df = df.set_index(index_cols)

    df_log = np.log2(df)
    df_log = df_log.replace(-np.inf, np.nan)
    df_log = df_log.dropna()

    number_of_proteins_in_common = df_log.shape[0]

    df_log_t = df_log.transpose()
    df_log_t.reset_index(inplace=True)
    df_log_t = df_log_t.rename(columns={"index": "channel_name"})

    X = df_log_t.drop("channel_name", axis=1)

    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(X)
    principalDf = pd.DataFrame(
        data=principalComponents, columns=["PC1", "PC2", "PC3"]
    )
    finalDf = pd.concat([principalDf, df_log_t[["channel_name"]]], axis=1)

    percent_df = pd.DataFrame(
        [
            {
                "principal_component": f"PC{i + 1}",
                "percent_explained": round(pca.explained_variance_ratio_[i] * 100, 2),
            }
            for i in range(3)
        ]
    )

    # `<state>-<fraction>_<replicate>` -> `<state>-<fraction>`
    finalDf["condition"] = finalDf["channel_name"].str.split("_").str[0]
    finalDf["state"] = finalDf["condition"].str.split("-").str[0]
    finalDf["fraction"] = finalDf["condition"].str.split("-").str[1]

    loadings_df = pd.DataFrame(
        data=np.transpose(pca.components_), columns=["PC1", "PC2", "PC3"]
    )
    loadings_df["variable"] = df_log.index.tolist()
    loadings_df = loadings_df.set_index("variable")

    finalDf.to_csv(out_dir / "pca_results.csv", index=False)
    loadings_df.to_csv(out_dir / "loadings_results.csv")
    percent_df.to_csv(out_dir / "percent_explained.csv", index=False)

    return loadings_df, finalDf, percent_df, number_of_proteins_in_common


# --------------------------------------------------------------------------------------
# Gene set enrichment analysis
# --------------------------------------------------------------------------------------

# MSigDB C5 GO biological process, the collection used in notebooks/01_whole_proteome.
GSEA_GMT = "c5.go.bp"


def perform_gsea_proteomics(ranked_proteins, gene_sets, output_dir):
    """Pre-ranked GSEA: 10000 permutations, seed 42, no plots, BH `FDR q-val`."""
    import gseapy as gp

    ranked_proteins = ranked_proteins.sort_values(
        ascending=True, by=ranked_proteins.columns[0]
    )

    return gp.prerank(
        rnk=ranked_proteins,
        gene_sets=gene_sets,
        threads=4,
        permutation_num=10000,
        outdir=str(output_dir),
        seed=42,
        no_plot=True,
    )


def gsea_ranking(long_df, condition, control):
    """log2FC ranking for one comparison, keyed on gene symbol.

    GSEA is rank-based, so the large global offset of the fraction contrasts (median
    log2(I/S) is about -3 to -4) does not affect the result: only the relative order of
    proteins matters. The three duplicated gene symbols are collapsed by mean.
    """
    sub = long_df[
        (long_df["condition"] == condition) & (long_df["control_condition"] == control)
    ]
    ranked = (
        sub.groupby("protein", as_index=True)["log2_FC"].mean().to_frame("log2_FC")
    )
    return ranked


def run_fraction_gsea(long_df, out_dir, gene_sets=None):
    """GSEA on each of the three insoluble/soluble fraction contrasts.

    Positive NES means the set is enriched among proteins with the *highest* log2(I/S),
    i.e. proteins gaining intensity in the insoluble fraction.

    Writes `<state>/gseapy.gene_set.prerank.report.csv` per state, the file the Rmd reads.
    """
    from gseapy import Msigdb

    out_dir = Path(out_dir)
    if gene_sets is None:
        gene_sets = Msigdb().get_gmt(GSEA_GMT)

    results = {}
    for condition, control, contrast_type in COMPARISONS:
        if contrast_type != "fraction":
            continue
        state = condition.split("-")[0]
        ranked = gsea_ranking(long_df, condition, control)
        res = perform_gsea_proteomics(ranked, gene_sets, out_dir / state)
        report = res.res2d.copy()
        report["state"] = state
        report["condition"] = condition
        report["control_condition"] = control
        results[state] = report
        print(
            f"{state:4s} {len(ranked)} proteins ranked, {len(report)} sets tested, "
            f"{(report['FDR q-val'] < 0.01).sum()} with FDR < 0.01"
        )

    combined = pd.concat(results.values(), ignore_index=True)
    combined.to_csv(out_dir / "gsea_all_fractions.csv", index=False)
    return results


# --------------------------------------------------------------------------------------
# Supplementary data export
# --------------------------------------------------------------------------------------
# `create_entry_cache` and `get_function` are copied verbatim from
# notebooks/05_low_input_proteomics/low_input.py so this module depends only on unipressed.


def get_function(x, cache):
    ret = set()
    if x in cache and "comments" in cache[x]:
        for i in cache[x]["comments"]:
            if "texts" in i and i["commentType"] == "FUNCTION":
                for j in i["texts"]:
                    ret.add(j["value"])
    return "|".join(list(ret))


def create_entry_cache(df):
    """Cache of UniProt entries, keyed by accession, for each unique protein in `df`."""
    from unipressed import UniprotkbClient

    chunk_size = 500
    entry_dict = dict()
    uniprots = set()  # use set so each identifier is unique
    for uniprot in set(df["uniprot"]):
        if uniprot is not None and uniprot != "None":
            uniprots.add(uniprot)

    # Break the list of ids into smaller lists to not overwhelm uniprot
    i = 0
    uniprots = list(uniprots)  # convert to list so we can subscript
    chunks = [uniprots[x : x + chunk_size] for x in range(0, len(uniprots), chunk_size)]
    print("Querying UniProt...")
    for chunk in chunks:
        chunk_num = i * chunk_size
        print("Retrieved " + str(chunk_num) + " entries out of " + str(len(uniprots)))
        i = i + 1
        entries = UniprotkbClient.fetch_many(chunk)
        for entry in entries:
            entry_dict[entry["primaryAccession"]] = entry
    print("Done querying UniProt.")
    return entry_dict


def drop_contaminants(df):
    """Contaminant and keratin rows, excluded from the supplementary tables.

    Matches the filter in low_input.write_percent_control_to_excel. These rows stay in the
    analysis outputs, so the supplementary counts are 2 lower than the volcano counts.
    """
    keep = ~(
        df["uniprot"].str.contains("contaminant", na=False)
        | df["description"].str.contains("Keratin", na=False)
    )
    return df[keep].copy()


def uniprot_functions(df, cache_csv=None):
    """Map each uniprot accession to its UniProt FUNCTION comment.

    Querying 5330 accessions takes a few minutes, so the result is cached to `cache_csv`
    and reused on later runs. Delete that file to re-query.
    """
    if cache_csv is not None:
        cache_csv = Path(cache_csv)
        if cache_csv.exists():
            cached = pd.read_csv(cache_csv, keep_default_na=False)
            missing = set(df["uniprot"]) - set(cached["uniprot"])
            if not missing:
                return cached

    cache = create_entry_cache(df)
    functions = pd.DataFrame(
        {
            "uniprot": sorted(set(df["uniprot"])),
        }
    )
    functions["uniprot_function"] = [
        get_function(u, cache) for u in functions["uniprot"]
    ]

    if cache_csv is not None:
        cache_csv.parent.mkdir(parents=True, exist_ok=True)
        functions.to_csv(cache_csv, index=False)
    return functions


# The Data S2 sheet this table ships as. `S2_6_TITLE` becomes row 1 of the sheet and the
# `contents` entry; `S2_6_DESCRIPTION` becomes row 2. Both are placeholder wording - edit
# them here and nowhere else.
S2_6_SHEET = "S2-6 Detergent solubility"
S2_6_TITLE = (
    "S2-6 Detergent solubility proteomics (TMT-exp) of activated, acutely, and "
    "chronically stimulated human T cells"
)
S2_6_DESCRIPTION = (
    "Percent of each protein recovered in the detergent-insoluble fraction, per replicate "
    "channel, with the median per cell state. Shares are relative: equal protein mass was "
    "labelled from each fraction, so compare across states rather than against 50."
)

# Front columns, in the order the S2-8..S2-11 low-input sheets use - `protein` first, not
# `uniprot`. Deliberately not ID_COLS order, so S2-6 matches its neighbours in Data S2.
SUPP_FRONT_COLS = ["protein", "uniprot", "description", "uniprot_function"]


def percent_insoluble_table(shares, cfg=REVISION, cache_csv=None):
    """Percent insoluble per replicate channel, plus the median per cell state.

    `shares` is what `fraction_shares` returns, so the percent itself is never recomputed
    here - this only selects the insoluble half, widens it and adds the medians.

    Replicate columns are `<state>_<replicate>` (`D2_d1_1`, ...), matching the
    `<condition>_<donor>_<duplicate>` naming of the low-input sheets alongside it in Data S2.

    CAVEAT: equal protein MASS was labelled per fraction and the revision run has no Final
    dilution factor, so a share is a relative position, not "X% of the protein is insoluble".
    """
    insoluble = drop_contaminants(shares[shares["fraction"] == "Insoluble"]).copy()

    # Pivot on `uniprot` alone: it is one row per protein here, and a multi-column index
    # would make pandas build the cartesian product of all three id columns.
    ids = insoluble[ID_COLS].drop_duplicates(subset="uniprot").set_index("uniprot")

    insoluble["channel"] = insoluble["state"] + "_" + insoluble["replicate"]
    replicate_cols = [f"{state}_{rep}" for state in STATES for rep in cfg.replicates]
    wide = insoluble.pivot(index="uniprot", columns="channel", values="percent").reindex(
        columns=replicate_cols
    )

    medians = (
        insoluble.groupby(["uniprot", "state"])["percent"]
        .median()
        .unstack("state")
        .reindex(columns=STATES)
        .rename(columns={state: f"{state} median" for state in STATES})
    )

    table = ids.join(wide).join(medians).reset_index()

    functions = uniprot_functions(table, cache_csv=cache_csv)
    table = table.merge(functions, on="uniprot", how="left")
    table = table[SUPP_FRONT_COLS + [c for c in table.columns if c not in SUPP_FRONT_COLS]]

    # Each median must sit inside its own state's replicate range; a median landing outside
    # would mean the two pivots disagreed about which channels belong to which state.
    for state in STATES:
        reps = table[[f"{state}_{rep}" for rep in cfg.replicates]]
        median = table[f"{state} median"]
        within = (median >= reps.min(axis=1) - 1e-9) & (median <= reps.max(axis=1) + 1e-9)
        assert within[median.notna()].all(), f"{state} median outside its replicate range"

    return table


def write_percent_insoluble_excel(shares, out_path, cfg=REVISION, cache_csv=None):
    """Write the Data S2-6 sheet: percent insoluble per channel and per-state medians.

    One sheet, styled as in low_input.write_percent_control_to_excel. Contaminant and
    keratin rows are dropped, and every row carries the UniProt FUNCTION comment.
    """
    import polars as pl
    import xlsxwriter

    table = percent_insoluble_table(shares, cfg=cfg, cache_csv=cache_csv)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    workbook = xlsxwriter.Workbook(out_path, {"nan_inf_to_errors": True})
    pl.from_pandas(table).write_excel(
        workbook=workbook,
        worksheet=S2_6_SHEET,
        table_style="Table Style Light 8",
        column_widths=150,
        freeze_panes=(1, 1),
    )
    workbook.close()

    return table
