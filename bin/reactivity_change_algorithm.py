import pandas as pd
import pathlib
import matplotlib.pyplot as plt
from functools import reduce

PEPTIDE_REQUIREMENT = 2

def median_or_ratio_filter(r, median_filter=True):
    """For proteins with PEPTIDE_REQUIREMENT or more quantified peptides, require that the maximum 
    peptide R ratio differ more than 2-fold from the minimum peptide R ratio. 
    If median filter is true, cysteines with five or more peptides must be more 
    than two fold from the median of all quantified cysteines for that protein"""
    r_values = r["percent_control"]
    residues = r["residue"]
    median = pd.Series(r_values).median()
    num_quantified_peptides = len(r_values)
    passing_residues = []
    passing_values = []
    if median_filter and (num_quantified_peptides >= 5):
        for i, value in enumerate(r_values):
            ratio_to_median = value / median
            if ratio_to_median > 2 or ratio_to_median < 0.5:
                passing_residues.append(residues[i])
                passing_values.append(value)
    elif num_quantified_peptides >= PEPTIDE_REQUIREMENT:
        for i, value in enumerate(r_values):
            ratio = max(r_values) / min(r_values)
            if ratio > 2 or ratio < 0.5:
                passing_values.append(value)
                passing_residues.append(residues[i])
    r["num_quantified_peptides"] = num_quantified_peptides
    r["TMT_ABPP_median_or_ratio_filter"] = passing_values
    r["residues_median_or_ratio_filter"] = passing_residues
    r["num_median_or_ratio_filter"] = len(passing_residues)
    return r


def expression_filter(r):
    """A cysteine passes the expression filter if its peptide R value differs more than 
    two-fold from the protein expression level measured by TMT-exp data"""
    r_values = r["percent_control"]
    residues = r["residue"]
    wp_percent_control = r["wp_percent_control"]
    passing_residues = []
    passing_values = []
    if pd.isnull(wp_percent_control):
        passing_residues = residues
        passing_values = r_values
    else:
        for i, value in enumerate(r_values):
            ratio_to_median = value / wp_percent_control
            if ratio_to_median > 2 or ratio_to_median < 0.5:
                passing_residues.append(residues[i])
                passing_values.append(value)
    r["TMT_ABPP_expression_filter"] = passing_values
    r["residues_expression_filter"] = passing_residues
    r["num_expression_filter"] = len(passing_residues)
    return r

def bio_replicate_variation_filter(r):
    """proteins are required to have at least one peptide R ratio within 1.5-fold
    of the protein expression level measured in TMT-exp experiments (if available)
    and are excluded if all peptide R ratios were greater than 2.0 or less than 0.5
    """
    passes = False
    if all(x > 200 for x in r["percent_control"]) or all(
        x < 50 for x in r["percent_control"]
    ):
        passes = False
    elif pd.isnull(r["wp_percent_control"]):
        passes = True
    else:
        passes = any(x / r["wp_percent_control"] < 1.5 for x in r["percent_control"])
    r["passes_bio_replicate_expression_variation_filter"] = passes
    return r

def manual_curation(r, curation):
    uniprot = r["uniprot"]
    protein_curation = curation[curation["uniprot"] == uniprot]
    if len(protein_curation) > 1:
        raise Exception("Redundant annotation for {}".format(uniprot))
    if len(protein_curation) < 1:
        return r
    curation_r = protein_curation.iloc[0]
    for curation_field in curation:
        r[curation_field] = curation_r[curation_field]
    return r


def reactivity_changes(r):
    reactivity_changes_vals = []
    reactivity_changes_residues = []
    r["two_peptide_no_expression_case"] = r[
        "num_quantified_peptides"
    ] == 2 and pd.isnull(r["wp_percent_control"])
    if (
        r["passes_bio_replicate_expression_variation_filter"]
        and r["curation_color"] != "red"
    ):
        reactivity_changes_vals = [
            x
            for x in r["TMT_ABPP_expression_filter"]
            if x in r["TMT_ABPP_median_or_ratio_filter"]
        ]
        reactivity_changes_residues = [
            x
            for x in r["residues_expression_filter"]
            if x in r["residues_median_or_ratio_filter"]
        ]
    curarted_residue = r["{}_curated_residue".format(r["condition"])]
    curated_vals = []
    curated_residues = []
    if type(curarted_residue) == float:
        curarted_residue = ""
    curation_residues = curarted_residue.split("|")
    if r["curation_color"] == "pink":
        for residue in curation_residues:
            if residue in reactivity_changes_residues:
                residue_index = reactivity_changes_residues.index(residue)
                curated_vals.append(reactivity_changes_vals[residue_index])
                curated_residues.append(residue)
        reactivity_changes_vals = curated_vals
        reactivity_changes_residues = curated_residues
    # green proteins are proteins not passing original filters which are manually
    # annotated as reactivity changes
    elif r["curation_color"] == "green":
        all_residues = r["residue"]
        all_rc_vals = r["percent_control"]
        for residue in curation_residues:
            if residue in all_residues:
                residue_index = all_residues.index(residue)
                curated_vals.append(all_rc_vals[residue_index])
                curated_residues.append(residue)
        reactivity_changes_vals = curated_vals
        reactivity_changes_residues = curated_residues

    r["TMT_ABPP_reactivity_changes"] = reactivity_changes_vals
    r["residues_reactivity_changes"] = reactivity_changes_residues
    r["num_reactivity_changes"] = len(reactivity_changes_residues)
    return r


def assign_directionality(r):
    """Label reactivity changes labeled as higher/lower by comparison of the fold change of the cysteine to
    matched whole proteome expression data if available, or otherwise by comparison to
    the median of cysteine fold change in reactivity profiling."""
    if pd.isnull(r["wp_percent_control"]):
        reference_percent_control = r["median_TMT_ABPP_percent_control"]
    else:
        reference_percent_control = r["wp_percent_control"]
    if r["TMT_ABPP_reactivity_changes"] > reference_percent_control:
        r["direction_of_reactivity_change"] = "Higher"
    else:
        r["direction_of_reactivity_change"] = "Lower"
    return r


def reactivity_change_algorithm(abpp_agg_df, curation_df):
    """run complete reactivity change algorithm
    ~4 minutes for normal dataset"""
    return (
        abpp_agg_df.apply(bio_replicate_variation_filter, axis=1)
        .apply(median_or_ratio_filter, axis=1)
        .apply(expression_filter, axis=1)
        .apply(manual_curation, curation=curation_df, axis=1)
        .apply(reactivity_changes, axis=1)
    )
