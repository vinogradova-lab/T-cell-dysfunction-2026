"""Assemble Supplementary Table 2 with the solubility table inserted as S2-6.

The previous submission has no solubility sheet and its S2-6 slot is taken by the low-input
mouse table. Inserting solubility at 6 shifts every sheet from there on up by one:

    S2-1 .. S2-5   unchanged
    S2-6           Detergent solubility          NEW
    S2-7 .. S2-10  the four low-input sheets     were S2-6 .. S2-9

That lands every sheet on the number the current manuscript cites it by: `Data S2-1` through
`Data S2-10`, with solubility at `Data S2-6`. The numbering is the previous submission's, not
this script's - inserting one sheet is the only change it makes to the sequence.

Every sheet's figure citation is also rewritten, both in Nature style and against the current
manuscript: the previous submission cites figures that no longer exist. See `NATURE_CITATIONS`.

`SOURCE` is the previous submission, `Data S2_20260528.xlsx`, which is read and never written.
An earlier build of this script used `Data S2_EV2.xlsx` from Dropbox, a divergent version
carrying an extra PCA-pathway-enrichment sheet and the low-input sheets in the other order;
neither matches the manuscript, so it is not the source. The new sheet's tab name and its
row-1/row-2 text come from `S2_6_SHEET` / `S2_6_TITLE` / `S2_6_DESCRIPTION` in `solubility.py`,
so the wording lives in one place. Tab names keep their `S2-<n>` form: Excel caps a tab name at
31 characters, which a spelled-out "Supplementary Table 2-<n> ..." would exceed.

Run in the `polars` env:

    conda run -n polars python bin/renumber_data_s2.py
"""

import re
import shutil
import sys
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.worksheet.table import Table, TableStyleInfo

REPO = Path(__file__).resolve().parent.parent
sys.path.append(str(REPO / "notebooks" / "06_solubility_proteomics"))

from solubility import (  # noqa: E402
    REVISION,
    S2_6_DESCRIPTION,
    S2_6_SHEET,
    S2_6_TITLE,
    STATES,
    fraction_percent_table,
)

SOURCE = REPO / "supp_data" / "Data S2_20260528.xlsx"
OUT = REPO / "supp_data" / "Supplementary_Table_2.xlsx"

SHARES_CSV = REPO / "data/solubility/revision/03_results/barplots/fraction_shares.csv"
FUNCTIONS_CSV = REPO / "data/solubility/revision/03_results/uniprot_functions.csv"

# The sheet solubility displaces, and everything below it. Renamed from the bottom up so a
# new name never collides with a name still in use.
FIRST_SHIFTED = 6
LAST_SHIFTED = 9

# Sheets above the shift point keep their tab name AND their row-1 title verbatim - including
# the stray "Data " prefix the source's S2-2 tab carries.
UNTOUCHED = [
    "contents",
    "S2-1 Unenriched proteomics",
    "Data S2-2 Bulk RNA vs WP",
    "S2-3 Reactivity changes",
    "S2-4 Cysteine reactivity vs WP",
    "S2-5 ATP add-back",
]

EXPECTED_ORDER = UNTOUCHED + [
    S2_6_SHEET,
    "S2-7 Low-input WP, mouse",
    "S2-8 Low-input WP, TILs",
    "S2-9 Low-input RC, inhibitors",
    "S2-10 Low-input WP, inhibitors",
]


# Every citation, re-derived from `Manuscript_Main_NI_20260904_clean_yellow.docx`, which has
# five main figures and eight Extended Data figures. The source workbook's numbers track an
# older draft with seven main figures - it cites `Figure 5`, `Figure 6` and `Figure S7`, none
# of which mean what they used to - so these are not a restyling of what it holds but a
# replacement of it. Nature style follows: `Fig.`/`Figs.`, supplementary as `Extended Data
# Fig. <n>`.
#
# Each is taken from the manuscript's own `Data S2-<n>` callouts plus the panels that
# unambiguously plot that dataset; `notebooks/06_solubility_proteomics/README.md` records which
# sentence supports which. The keys below are SOURCE tab names, i.e. before the shift, so a
# key's number is one lower than the sheet's final number from S2-6 on.
#
# The new solubility sheet is absent because its citation comes from `S2_6_TITLE`, which is
# already in this style.
NATURE_CITATIONS = {
    "S2-1 Unenriched proteomics": "(related to Fig. 1 and Extended Data Figs. 2, 6)",
    "Data S2-2 Bulk RNA vs WP": "(related to Extended Data Fig. 2)",
    "S2-3 Reactivity changes": "(related to Figs. 1, 2 and Extended Data Fig. 2)",
    "S2-4 Cysteine reactivity vs WP": "(related to Fig. 1 and Extended Data Fig. 2)",
    "S2-5 ATP add-back": "(related to Figs. 2, 3 and Extended Data Figs. 3, 5)",
    "S2-6 Low-input WP, mouse": "(related to Fig. 4 and Extended Data Fig. 7)",
    "S2-7 Low-input WP, TILs": "(related to Fig. 4 and Extended Data Fig. 7)",
    "S2-8 Low-input RC, inhibitors": "(related to Fig. 5 and Extended Data Fig. 8)",
    "S2-9 Low-input WP, inhibitors": "(related to Fig. 5 and Extended Data Fig. 8)",
}

# The trailing `(related to ...)` clause, which is the only part of a title this rewrites -
# the body wording is left alone, including where a sheet's row 1 and its contents row word
# it differently (S2-1 says "quiescent" in one and "acutely" in the other).
CITATION_RX = re.compile(r"\s*\(related to [^)]*\)\s*$", re.I)


def restyle(text, citation):
    """Swap a title's trailing citation clause for `citation`, appending if it had none."""
    return f"{CITATION_RX.sub('', text).rstrip()} {citation}"


def restyle_citations(wb):
    """Rewrite every sheet's row-1 citation, and its contents row, into Nature style."""
    contents = wb["contents"]
    # Keyed by the `S2-<n>` token, not the whole tab name: the S2-2 tab carries a stray
    # "Data " prefix that its contents row does not, so the two strings differ.
    rows = {
        contents.cell(row=r, column=2).value.split()[0]: r
        for r in range(2, contents.max_row + 1)
        if contents.cell(row=r, column=2).value
    }

    for source_name, citation in NATURE_CITATIONS.items():
        name = source_name if source_name in UNTOUCHED else bump_name(source_name)
        ws = wb[name]
        ws["A1"] = restyle(ws["A1"].value, citation)

        row = rows[name.removeprefix("Data ").split()[0]]
        cell = contents.cell(row=row, column=3)
        cell.value = restyle(cell.value, citation)
        print(f"  {name:32s} {citation}")


def bump(text, old, new):
    """Rewrite a leading `S2-<old>` token, leaving the rest of the string untouched."""
    bumped, n = re.subn(rf"^S2-{old}\b", f"S2-{new}", text)
    assert n == 1, f"no leading S2-{old} in {text!r}"
    return bumped


def bump_name(name):
    """The shifted tab name a source tab becomes: `S2-6 ATP add-back` -> `S2-7 ...`."""
    old = int(re.match(r"^S2-(\d+)\b", name).group(1))
    return bump(name, old, old + 1)


def cell_values(table):
    """Rows of plain Python scalars - openpyxl cannot write numpy types or NaN."""
    frame = table.astype(object).where(table.notna(), None)
    return frame.itertuples(index=False, name=None)


def build_sheet(wb, table, index):
    """Add the S2-6 sheet in the workbook's house style: title, description, header, data."""
    ws = wb.create_sheet(title=S2_6_SHEET, index=index)
    ws["A1"] = S2_6_TITLE
    ws["A2"] = S2_6_DESCRIPTION
    # Writing A2 leaves openpyxl's append cursor on row 2, so the header lands on row 3 -
    # the same layout as every other sheet, where the two title rows precede the header.
    ws.append(list(table.columns))
    for row in cell_values(table):
        ws.append(list(row))

    # Match the low-input sheets alongside it, which carry a striped table over the data.
    last_col = ws.cell(row=3, column=table.shape[1]).column_letter
    ws.add_table(
        Table(
            displayName="Frame07",
            ref=f"A3:{last_col}{3 + len(table)}",
            tableStyleInfo=TableStyleInfo(showRowStripes=True),
        )
    )
    return ws


def rewrite_contents(ws):
    """Insert the S2-6 row and renumber every entry below it.

    Row 1 is the header; entry `n` sits on row `n + 1`. The `#` column is written back in
    whatever type the sheet already uses, so a text-formatted column stays text.
    """
    number_type = type(ws.cell(row=2, column=1).value)
    row = FIRST_SHIFTED + 1

    ws.insert_rows(row)
    ws.cell(row=row, column=1).value = number_type(FIRST_SHIFTED)
    ws.cell(row=row, column=2).value = S2_6_SHEET
    # The description column repeats the sheet's row-1 title without its `S2-6 ` prefix,
    # which is what every other entry holds - not the row-2 methods sentence.
    ws.cell(row=row, column=3).value = S2_6_TITLE.removeprefix(f"S2-{FIRST_SHIFTED} ")

    for old in range(LAST_SHIFTED, FIRST_SHIFTED - 1, -1):
        title = ws.cell(row=old + 2, column=2)
        ws.cell(row=old + 2, column=1).value = number_type(old + 1)
        title.value = bump(title.value, old, old + 1)


def sample(ws, first_row, n=8):
    """The first `n` rows from `first_row`, with trailing blanks trimmed.

    Two of the low-input sheets declare 16384 columns, and the round trip materialises the
    empty ones, so a raw row-tuple comparison would fail on padding alone.
    """
    rows = []
    for values in ws.iter_rows(min_row=first_row, max_row=first_row + n - 1, values_only=True):
        row = list(values)
        while row and row[-1] is None:
            row.pop()
        rows.append(row)
    return rows


def verify(path, source_path, source_rows):
    """Reopen the written workbook and check the numbering and the carried-over data."""
    wb = load_workbook(path)
    try:
        assert wb.sheetnames == EXPECTED_ORDER, wb.sheetnames

        for name in wb.sheetnames:
            if name == "contents":
                continue
            title = wb[name]["A1"].value
            number = name.removeprefix("Data ").split()[0]
            assert title.startswith(number), f"{name}: row 1 reads {title[:40]!r}"

        # Every sheet carried over must still hold the data it held in the source, compared
        # from row 3 - below the title and description rows this script rewrites.
        for old_name, rows in source_rows.items():
            new_name = old_name if old_name in UNTOUCHED else bump_name(old_name)
            assert sample(wb[new_name], 3) == rows, f"{old_name}: data changed"
        print(f"  {len(source_rows)} sheets verified against {source_path.name}")

        # No `Figure`/`Figures` spelled out and no `S<n>` figure token may survive inside a
        # citation clause, in either place it appears. This is what proves the conversion was
        # complete rather than merely attempted.
        stale = re.compile(r"\bFigures?\b|\bS\d+\b")
        for name in wb.sheetnames:
            if name == "contents":
                continue
            found = CITATION_RX.search(wb[name]["A1"].value or "")
            assert found, f"{name}: row 1 has no citation clause"
            assert not stale.search(found.group(0)), f"{name} row 1: {found.group(0).strip()}"

        for r in range(2, wb["contents"].max_row + 1):
            text = wb["contents"].cell(row=r, column=3).value or ""
            found = CITATION_RX.search(text)
            assert found, f"contents row {r} has no citation clause"
            assert not stale.search(found.group(0)), (
                f"contents row {r}: {found.group(0).strip()}"
            )
        print("  citations verified: no `Figure`/`S<n>` left in any clause")

        contents = wb["contents"]
        # Entry n sits on row n + 1, and there is one entry per sheet but the `contents` tab.
        n_entries = len(EXPECTED_ORDER) - 1
        entries = [
            (contents.cell(row=r, column=1).value, contents.cell(row=r, column=2).value)
            for r in range(2, n_entries + 2)
        ]
        assert [int(n) for n, _ in entries] == list(range(1, n_entries + 1)), entries
        for (_, title), tab in zip(entries, EXPECTED_ORDER[1:]):
            assert title.split()[0] == tab.removeprefix("Data ").split()[0], (title, tab)
    finally:
        wb.close()

    print(
        f"verified: {len(EXPECTED_ORDER)} sheets, "
        f"contents 1-{len(EXPECTED_ORDER) - 1} consistent with the tabs"
    )


def main():
    shares = pd.read_csv(SHARES_CSV)
    table = fraction_percent_table(shares, cfg=REVISION, cache_csv=FUNCTIONS_CSV)
    print(f"{S2_6_SHEET}: {table.shape[0]} proteins x {table.shape[1]} cols")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(SOURCE, OUT)
    print(f"copied {SOURCE.name} -> {OUT}")

    wb = load_workbook(OUT)

    # Recorded before the round trip so `verify` can prove no sheet lost its data. Sampled
    # from row 3 for every sheet, below the two title rows: row 1 is now rewritten on all of
    # them, not just the shifted ones, since the citation restyle touches every sheet.
    source_rows = {
        name: sample(wb[name], 3) for name in wb.sheetnames if name != "contents"
    }

    for old in range(LAST_SHIFTED, FIRST_SHIFTED - 1, -1):
        ws = next(wb[n] for n in wb.sheetnames if n.startswith(f"S2-{old} "))
        renamed = bump(ws.title, old, old + 1)
        ws["A1"] = bump(ws["A1"].value, old, old + 1)
        print(f"  {ws.title!r} -> {renamed!r}")
        ws.title = renamed

    build_sheet(wb, table, index=FIRST_SHIFTED)
    rewrite_contents(wb["contents"])
    restyle_citations(wb)

    wb.save(OUT)
    wb.close()
    print(f"wrote {OUT}")

    verify(OUT, SOURCE, source_rows)


if __name__ == "__main__":
    main()
