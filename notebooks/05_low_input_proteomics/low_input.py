import requests
from time import sleep
import scipy.stats as stat
import math
import matplotlib.pyplot as plt
import seaborn as sns
import polars as pl
import polars.selectors as cs
import re
from pathlib import Path
import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import plotly.express as px
import gseapy as gp
from gseapy import Msigdb
from unipressed import UniprotkbClient

# for cleaning channel names
ORIGINAL_LABELS = [
    "_Stim",
    "_CK",
    "_TG",
    "_BSJ",
    "_ISRIB",
    "_Thapsigargin",
    "_Gossypetin",
    "D8C_1",
    "D8C_2",
    "D8A_1",
    "D8A_2",
]
CLEANED_LABELS = [
    " + stim",
    " + CK",
    " + Thapsigargin",
    " + BSJ",
    " + ISRIB",
    " + Thapsigargin",
    " + Gossypetin",
    "D8C + DMSO_1",
    "D8C + DMSO_2",
    "D8A + DMSO_1",
    "D8A + DMSO_2",
]

# defining TPs list not subject to 2 peptide and 2 replicate filters
TRUE_POSITIVE_LIST = ["PDCD1","LAG3","TIGIT","HAVCR2","CD39","NR4A1","NR4A3","GZMB","GSR","CD62L","CTLA4","PRKCH","PRKCQ","SAMHD1"]
MOUSE_TRUE_POSITIVE_LIST = ["Pdcd1","Lag3","Tigit","Havcr2","Cd39","Nr4a1","Nr4a3","Gzmb","Gsr","Cd62l","Ctla4","Prkch","Prkcq","Samhd1"]
CURATED_FILTERS = True # switch to turn curated filtering on/off
curated_filter_expression = (pl.col("protein").is_in(TRUE_POSITIVE_LIST) & CURATED_FILTERS) 

def get_pca_plot(df, index_cols, title_string, out_dir):
    """Principal component analysis.
    Scales data with log. Runs and plots PCA with scikit learn
    """
    out_dir.mkdir(exist_ok=True)
    # drop na
    df = df.dropna()
    df = df.set_index(index_cols)
    # log2 transform
    df_log = np.log2(df)

    df_log = df_log.replace(-np.inf, np.nan)
    df_log = df_log.dropna()

    number_of_proteins_in_common = df_log.shape[0]
    print(number_of_proteins_in_common)

    # transpose table
    df_log_t = df_log.transpose()
    df_log_t.reset_index(inplace=True)
    df_log_t = df_log_t.rename(columns={"index": "channel_name"})

    # get X and y
    X = df_log_t.drop("channel_name", axis=1)

    # get PCA 2, fit transform, get df
    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(X)
    principalDf = pd.DataFrame(
        data=principalComponents,
        columns=[
            "principal component 1",
            "principal component 2",
            "principal component 3",
        ],
    )
    # add channel name
    finalDf = pd.concat([principalDf, df_log_t[["channel_name"]]], axis=1)

    # get PCA %
    pca_1_percent = round((pca.explained_variance_ratio_[0] * 100), 2)
    pca_2_percent = round((pca.explained_variance_ratio_[1] * 100), 2)
    pca_3_percent = round((pca.explained_variance_ratio_[2] * 100), 2)

    percent_df = pd.DataFrame(
        [
            {"principal_component": "PC1", "percent_explained": pca_1_percent},
            {"principal_component": "PC2", "percent_explained": pca_2_percent},
            {"principal_component": "PC3", "percent_explained": pca_3_percent},
        ]
    )

    # add condition
    finalDf["condition"] = finalDf["channel_name"]

    finalDf["condition"] = finalDf["condition"].str.split("_").str[0]
    title_info = title_string + " (" + str(number_of_proteins_in_common) + " Proteins)"

    fig = px.scatter(
        finalDf,
        x="principal component 1",
        y="principal component 2",
        hover_data=["channel_name"],
        color=finalDf["condition"],
        labels={
            "principal component 1": "PC1 ({}%)".format(pca_1_percent),
            "principal component 2": "PC2 ({}%)".format(pca_2_percent),
            "condition": "Condition",
        },
        title=title_info,
        template="plotly_white",
    )

    fig.update_traces(marker=dict(size=12), selector=dict(mode="markers"))
    fig.update_layout(height=600, width=600, showlegend=True, legend_title_text="")

    loadings_df = pd.DataFrame(
        data=np.transpose(pca.components_), columns=["PC1", "PC2", "PC3"]
    )
    loadings_df["variable"] = df_log.index.tolist()

    loadings_df = loadings_df.set_index("variable")

    finalDf.to_csv(out_dir / "pca_results.csv")
    loadings_df.to_csv(out_dir / "loadings_results.csv")
    percent_df.to_csv(out_dir / "percent_explained.csv")

    return (fig, loadings_df, finalDf, percent_df)


def filter_to_two_reps(
    curated_filter_expression, id_columns, allow_one_rep, channel_ratio_df
):
    if "identifier" in id_columns:
        channel_ratio_df = channel_ratio_df.with_columns(
            pl.concat_str(pl.col("protein"), pl.col("residue")).alias("key")
        )
    else:
        channel_ratio_df = channel_ratio_df.with_columns(pl.col("protein").alias("key"))

    print(channel_ratio_df.columns)
    # filter to min two reps or TP proteins
    two_rep_proteins = set(
        channel_ratio_df.group_by(id_columns + ["key"])
        .agg(pl.col("donor").n_unique())
        .filter(pl.col("donor") > 1)["key"]
    )

    channel_ratio_df = (
        channel_ratio_df.filter(
            (pl.col("key").is_in(two_rep_proteins))
            | curated_filter_expression
            | allow_one_rep,
        )
        .with_columns(
            pl.concat_str(
                ["condition", "donor", "technical_replicate"], separator="_"
            ).alias("clean_channel_name")
        )
        .drop("key")
    )

    return channel_ratio_df


# uniprot function annotation
def get_function(x, cache):
    ret = set()
    if x in cache and "comments" in cache[x]:
        for i in cache[x]["comments"]:
            if "texts" in i and i["commentType"] == "FUNCTION":
                for j in i["texts"]:
                    ret.add(j["value"])
    return "|".join(list(ret))


def create_entry_cache(df):
    """Given a dataframe of proteins,
    create a cache of UniProt entries with each unique
    protein. Remember to remove contaminants
    """
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


def channel_ratio_to_percent_control(
    channel_ratio_df, control_condition, grouping_columns
):
    control_condition_medians = (
        channel_ratio_df.filter(condition=control_condition)
        .group_by(grouping_columns)
        .agg(pl.col("value").median().alias(f"{control_condition}_median"))
    )

    percent_controls = channel_ratio_df.join(
        other=control_condition_medians, on=grouping_columns
    ).with_columns(
        ((pl.col("value") / pl.col(f"{control_condition}_median")) * 100).alias(
            "percent_control"
        )
    )

    return percent_controls


def write_percent_control_to_excel(
    df, file_name, id_columns=["protein", "uniprot", "description"]
):
    """Filter out contaminants and keratins, add uniprot function,
    and write to excel file"""
    df = df.select(id_columns + ["percent_control", "clean_channel_name"]).pivot(
        index=id_columns,
        values="percent_control",
        on="clean_channel_name",
    )
    # filter out contaminants and keratins
    df = df.filter(
        pl.col("uniprot").str.contains("contaminant").not_(),
        pl.col("description").str.contains("Keratin").not_(),
    )

    cache = create_entry_cache(df)

    formatted_df = df.with_columns(
        pl.col("uniprot")
        .map_elements(lambda x: get_function(x, cache), return_dtype=pl.String)
        .alias("uniprot_function")
    ).select(
        pl.col(id_columns),
        "uniprot_function",
        cs.contains("D"),
    )

    formatted_df.write_excel(
        workbook=file_name,
        table_style=f"Table Style Light 8",
        column_widths=150,
        freeze_panes=(1, 1),
    )


def channel_ratio_boxplot(df, ax, title, color_palette):
    sns.boxplot(
        data=df,
        x="value",
        y="clean_channel_name",
        log_scale=True,
        showfliers=False,
        hue="condition",
        ax=ax,
        legend=False,
        palette=color_palette,
    )
    ax.set_xlabel("log10(channel ratio)")
    ax.set_ylabel("")
    ax.set_title(title)


def median_normalize_data(df, height=5, color_palette=None):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, height), sharey=True, sharex=True)

    channel_ratio_boxplot(df, ax1, "Before normalization", color_palette)

    channel_medians = df.group_by(["donor", "clean_channel_name"]).agg(
        pl.col("value").median().alias("channel_median")
    )

    overall_median = df.select(pl.col("value").median()).item()

    df_normalized = (
        df.join(channel_medians, on=["donor", "clean_channel_name"])
        .with_columns(
            (overall_median / pl.col("channel_median")).alias("normalization_factor")
        )
        .with_columns((pl.col("value") * pl.col("normalization_factor")).alias("value"))
        .drop("channel_median")
    )

    channel_ratio_boxplot(
        df_normalized, ax2, "After median normalization", color_palette
    )

    plt.tight_layout()

    plt.show()
    return df_normalized


def read_channel_ratios(
    data_dir,
    donors,
    curated_filter_expression,
    id_columns=["uniprot", "protein", "description"],
    fn="03_combfiles_forpca_channelratio_or_rawsignal_wp.csv",
    folder_suffix="_WP_one peptide/03_combined_files",
    allow_one_rep=False,
    field_order=["condition", "technical_replicate", "drop"],
):

    channel_ratio_dfs = []

    for donor in donors.keys():
        file_name = data_dir / (donor + folder_suffix) / donor / fn
        print(file_name)
        donor_channel_ratios = pl.read_csv(file_name).drop("")

        pepnum_columns = donor_channel_ratios.select(cs.contains("pepNum")).columns
        if pepnum_columns:
            peptide_number_column = donor_channel_ratios.select(
                cs.contains("pepNum")
            ).columns[0]
            donor_channel_ratios = (
                donor_channel_ratios.rename(
                    mapping={peptide_number_column: "peptide_number"}
                )
                .unpivot(index=id_columns + ["peptide_number"])
                .filter(  # filter to two peptide measurements or TP proteins
                    (pl.col("peptide_number") > 1) | curated_filter_expression
                )
            )
        else:
            donor_channel_ratios = donor_channel_ratios.unpivot(index=id_columns)
        donor_channel_ratios = donor_channel_ratios.with_columns(
            # rename conditions
            pl.col("variable").str.replace_many(
                ORIGINAL_LABELS,
                CLEANED_LABELS,
            )
        )

        if "donor" not in field_order:
            donor_channel_ratios = donor_channel_ratios.with_columns(
                pl.lit(donors[donor]).alias("donor"),
            )

        channel_ratio_dfs.append(donor_channel_ratios)


    channel_ratio_df = (
        pl.concat(channel_ratio_dfs)
        .with_columns(
            pl.col("variable")
            .str.splitn("_", len(field_order))
            .struct.rename_fields(field_order)
        )
        .unnest("variable")
        .select(~cs.contains("drop"))
        .filter(
            # filter out unused TMT channels
            pl.col("condition")
            != "No",
        )
    )

    channel_ratio_df = filter_to_two_reps(
        curated_filter_expression, id_columns, allow_one_rep, channel_ratio_df
    )

    sns.boxplot(
        data=channel_ratio_df,
        x="value",
        y="clean_channel_name",
        log_scale=True,
        fliersize=1,
    )

    return channel_ratio_df


def get_p_value(row, cond_1, cond_2):
    ttest_result = stat.ttest_ind(row[cond_1], row[cond_2], nan_policy="omit")
    return ttest_result[1]


def get_expr(row):
    p_value_column = "-log10_pval"
    p_value_cutoff = -1 * math.log10(0.05)
    fc_cutoff = math.log2(1.5)
    if (row["log2_FC"] > fc_cutoff) & (row[p_value_column] > p_value_cutoff):
        return "Significant Up"
    if (row["log2_FC"] > fc_cutoff) & (row[p_value_column] < p_value_cutoff):
        return "Not Significant Up"
    if (row["log2_FC"] < -fc_cutoff) & (row[p_value_column] > p_value_cutoff):
        return "Significant Down"
    if (row["log2_FC"] < -fc_cutoff) & (row[p_value_column] < p_value_cutoff):
        return "Not Significant Down"
    if (
        (row["log2_FC"] > -fc_cutoff)
        & (row["log2_FC"] < fc_cutoff)
        & (row[p_value_column] > p_value_cutoff)
    ):
        return "Significant but <1.5 FC"
    else:
        return "Not Significant"


def get_volcano_plot_treatment_vs_control(
    conditions_list, control_labelling, df, file_name, folder_path
):
    volcano_plot_list = []
    volcano_df_list = []
    long_format_list = []
    index_cols = ["uniprot", "protein", "description"]
    df = df.set_index(index_cols)

    for condition in conditions_list:
        copy_df = df.copy()
        list_cond_1 = copy_df.filter(like=condition).columns.tolist()
        list_cond_2 = copy_df.filter(like=control_labelling).columns.tolist()
        if len(list_cond_1) <= 1 or len(list_cond_2) <= 1:
            print(
                f"Condition {condition} does not have enough replicates to be shown in volcano plot!"
            )
            print(list_cond_1)
            print(list_cond_2)
            continue

        for labelling in [control_labelling, condition]:
            list_columns_labelling = copy_df.filter(like=labelling).columns.tolist()
            labelling_df = copy_df[list_columns_labelling]
            copy_df["mean_" + labelling] = labelling_df.mean(axis=1)

        copy_df["FC"] = (
            copy_df["mean_" + condition] / copy_df["mean_" + control_labelling]
        )

        idx_cond_1 = copy_df.columns.get_indexer(list_cond_1)
        idx_cond_2 = copy_df.columns.get_indexer(list_cond_2)
        copy_df["p_value"] = copy_df.apply(
            get_p_value,
            axis=1,
            args=(idx_cond_1, idx_cond_2),
        )

        copy_df["log2_FC"] = np.log2(copy_df["FC"])
        volcano_df = copy_df[["p_value", "log2_FC"]]

        volcano_df = volcano_df.dropna()
        volcano_df["-log10_pval"] = -1 * np.log10(volcano_df["p_value"])
        # Using BH correction here (previous versions used Boneferroni)
        volcano_df["-log10_pval_adj"] = -1 * np.log10(
            stat.false_discovery_control(volcano_df["p_value"])
        )
        volcano_df["Regulation"] = volcano_df.apply(get_expr, axis=1)

        volcano_df["Regulation"] = volcano_df["Regulation"].astype("category")

        volcano_df = volcano_df.reset_index()
        volcano_df["name"] = volcano_df["protein"]

        title_name = (
            file_name
            + " - "
            + condition
            + " vs. "
            + control_labelling
            + " ("
            + str(len(volcano_df))
            + " Proteins)"
        )
        long_df = volcano_df.copy()
        long_df["condition"] = condition
        long_df["control_condition"] = control_labelling
        long_format_list.append(long_df)
    long_df = pd.concat(long_format_list, axis=0)

    return volcano_df, copy_df, long_df


UNIPROT_API = "https://rest.uniprot.org"


def convert_to_mouse_uniprot(
    human_uniprots, batch_size=100, source_organism="9606", destination_organism="10090"
):
    """Convert human UniProt IDs to mouse orthologs efficiently via batched API calls"""

    if not human_uniprots:
        return []

    # Filter out NaN values and ensure all are strings
    uniprot_list = [str(uid).strip() for uid in human_uniprots if pd.notna(uid)]

    if not uniprot_list:
        return []

    print(f"Converting {len(uniprot_list)} human proteins to mouse orthologs...")

    gene_to_human = {}  # Track which human IDs map to which genes
    mouse_proteins = set()

    # Process in batches
    for i in range(0, len(uniprot_list), batch_size):
        batch = uniprot_list[i : i + batch_size]

        try:
            # Batch query for human proteins to get gene names
            human_query = " OR ".join([f"accession:{uid}" for uid in batch])
            search_url = f"{UNIPROT_API}/uniprotkb/search"
            params = {
                "query": f"({human_query}) AND organism_id:{source_organism}",
                "fields": "accession,gene_primary",
                "format": "tsv",
                "size": batch_size,
            }

            response = requests.get(search_url, params=params)
            response.raise_for_status()

            # Parse gene names from batch
            gene_names = set()
            for line in response.text.strip().split("\n")[1:]:  # Skip header
                if line:
                    parts = line.split("\t")
                    if len(parts) >= 2 and parts[1]:
                        gene_name = parts[1].strip()
                        if gene_name:
                            gene_names.add(gene_name)
                            gene_to_human[gene_name] = parts[0]

            # Batch query for mouse orthologs if we have gene names
            if gene_names:
                sleep(0.3)  # Rate limiting

                mouse_query = " OR ".join([f"gene:{gene}" for gene in gene_names])
                params = {
                    "query": f"({mouse_query}) AND organism_id:{destination_organism} AND reviewed:true",
                    "fields": "accession",
                    "format": "tsv",
                    "size": 500,
                }

                response = requests.get(search_url, params=params)
                response.raise_for_status()

                # Collect mouse proteins
                for line in response.text.strip().split("\n")[1:]:  # Skip header
                    if line:
                        mouse_id = line.strip()
                        if mouse_id:
                            mouse_proteins.add(mouse_id)

            print(
                f"  Processed batch {i//batch_size + 1}/{(len(uniprot_list)-1)//batch_size + 1}: {len(mouse_proteins)} proteins found"
            )
            sleep(0.3)  # Rate limiting

        except Exception as e:
            print(f"  Error processing batch {i//batch_size + 1}: {e}")
            continue

    result = list(mouse_proteins)

    return result


def annotate_functional_categories(volcano_df):
    """Annotate the volcano_data plot with functional categories. See 'protein_lists/README.md' for
    further documentation
    """
    protein_list_dir = Path(
        "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Vinogradova Laboratory/Henry_data processing/01_data_analysis_folders/03_Exhaustion/exhaustion/exhaustion_notebook/01_wp_analysis/"
    )

    fun_groups = pd.read_csv(
        protein_list_dir / "protein_lists/group_annotations_final.csv"
    )
    for i, row in fun_groups.iterrows():
        column_name = row["column_name"]
        file_name = row["file_name"]
        print(f"Processing {column_name}...")

        protein_table = pd.read_csv(protein_list_dir / file_name)
        functional_proteins = set(protein_table["uniprot"])

        print(f"  Converting {len(functional_proteins)} proteins to mouse UniProts...")
        mouse_functional_proteins = convert_to_mouse_uniprot(functional_proteins)
        print(mouse_functional_proteins)
        print(f"  Found {len(mouse_functional_proteins)} mouse orthologs")

        volcano_df[column_name] = volcano_df.index.isin(mouse_functional_proteins)

    return volcano_df
