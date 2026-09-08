"""Add the two detergent-solubility revision runs to the PRIDE submission Metadata sheet.

`Metadata.xlsx` is the PRIDE submission's per-run description: one row per MS run, giving the
TMT plex and the channel -> sample map. The copy downloaded from the current submission
(`SOURCE`) carries 24 rows - whole proteome, reactivity, metabolite add-back and the low-input
runs - and nothing at all from the detergent solubility assay. The revision run's two files,
`HK3-129` (12-plex, acquired 2026-08-22) and `HK3-130` (16-plex, 2026-08-23), are being added
to the deposition, so this script appends their two rows.

Only those two. `HK3-124`, the original solubility run, is deliberately left out.

The channel names are read from `METADATA_CSV`, the upstream pipeline's own channel -> sample
map, rather than retyped here. That file has one trap: HK3-129's `Rename channels` row lists
16 names when the run has only 12 channels, the last four being `No_d1_*` padding for a plex
the run does not have. Zipping the two rows truncates them, which is right - the header of
`03_combfiles_forpca_channelratio_or_rawsignal_wp.csv` confirms HK3-129 has no `No_*` columns
while HK3-130 has four.

Edited with openpyxl rather than rebuilt with polars `write_excel`: appending in place keeps
the sheet's Excel Table, its style, the frozen header row and the column widths exactly as
PRIDE served them, where a round trip through polars would restyle the whole sheet.

`SOURCE` is read and never written; the edit lands on `OUT`, which is what gets uploaded.
`*.xlsx` is gitignored repo-wide, so only this script is tracked.

Run in the `polars` env:

    conda run -n polars python bin/add_solubility_to_pride_metadata.py
"""

import csv
import shutil
from pathlib import Path

from openpyxl import load_workbook

REPO = Path(__file__).resolve().parent.parent

SOURCE = Path.home() / "Downloads" / "Metadata.xlsx"
OUT = REPO / "supp_data" / "Metadata.xlsx"

METADATA_CSV = Path(
    "/Users/henrysanford/dev/test_data/texh_revisions/detergent_proteomics/input"
    "/01_processed_files/metadata.csv"
)

SHEET = "Sheet1"
TABLE = "Frame0"

# The processed-file stem each run is keyed by in `METADATA_CSV`, and the row it becomes.
# `experiment` follows the census filename, as the `_WP` / `_Reactivity` rows already do.
RUNS = [
    ("processed_census-out_20260822_HK3-129_WP", "20260822_HK3-129_WP"),
    ("processed_census-out_20260823_HK3-130_WP", "20260823_HK3-130_WP"),
]

EXPERIMENT_TYPE = "detergent solubility proteomics"
SPECIES = "human"
IN_VITRO_OR_VIVO = "in vitro"
CONDITIONS = (
    "D2-soluble, D2-insoluble, D8A-soluble, D8A-insoluble, D8C-soluble, D8C-insoluble"
)

FIRST_TAG_COL = 7  # column G, `tag_126.127726`; A-F are the six description columns.


def read_channel_maps():
    """`{file stem: (plex, {tag: sample name})}` from the pipeline's metadata.csv.

    Each run is a block of four lines: `<stem>,<plex>`, the control replicates, the channel
    tags, and the renamed samples. Zipping tags to names drops HK3-129's four padding names.
    """
    with METADATA_CSV.open() as fh:
        rows = [row for row in csv.reader(fh) if any(field.strip() for field in row)]

    maps = {}
    for i, row in enumerate(rows):
        stem = row[0]
        if not stem.startswith("processed_census-out_"):
            continue
        plex = int(row[1])
        tags = [f for f in rows[i + 2][1:] if f]
        names = [f for f in rows[i + 3][1:] if f]
        assert rows[i + 2][0] == "Channels", rows[i + 2][0]
        assert rows[i + 3][0] == "Rename channels", rows[i + 3][0]
        assert len(tags) == plex, (stem, len(tags), plex)
        assert len(names) >= plex, (stem, len(names), plex)
        maps[stem] = (plex, dict(zip(tags, names)))
    return maps


def build_row(header, experiment, plex, channels):
    """One sheet row, in the header's column order; unused tags stay empty."""
    values = {
        "experiment": experiment,
        "experiment type": EXPERIMENT_TYPE,
        "species": SPECIES,
        "in-vitro or in-vivo": IN_VITRO_OR_VIVO,
        "conditions": CONDITIONS,
        "plex": plex,
    }
    unknown = set(channels) - set(header)
    assert not unknown, f"{experiment}: tags absent from the sheet: {sorted(unknown)}"
    return [values.get(col, channels.get(col)) for col in header]


def verify(source_rows, new_rows):
    """Reopen the written workbook and prove nothing but the two new rows changed."""
    wb = load_workbook(OUT)
    try:
        ws = wb[SHEET]
        n = len(source_rows) + len(new_rows)
        assert ws.max_row == n, ws.max_row
        assert ws.max_column == len(source_rows[0]), ws.max_column

        table = ws.tables[TABLE]
        assert table.ref == f"A1:V{n}", table.ref
        assert table.tableStyleInfo.name == "TableStyleLight8", table.tableStyleInfo.name
        assert ws.freeze_panes == "A2", ws.freeze_panes

        written = [[c.value for c in row] for row in ws.iter_rows()]
        assert written[: len(source_rows)] == source_rows, "an existing row changed"
        assert written[len(source_rows) :] == new_rows, "a new row did not round-trip"

        for row, expected in zip(written[len(source_rows) :], new_rows):
            plex = row[header_index(written[0], "plex")]
            tags = [v for v in row[FIRST_TAG_COL - 1 :] if v is not None]
            assert isinstance(plex, int), (row[0], type(plex))
            assert len(tags) == plex, (row[0], len(tags), plex)
    finally:
        wb.close()

    print(f"verified: {len(source_rows) - 1} existing rows untouched, {len(new_rows)} added")


def header_index(header, name):
    return header.index(name)


def main():
    maps = read_channel_maps()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(SOURCE, OUT)
    print(f"copied {SOURCE} -> {OUT}")

    wb = load_workbook(OUT)
    try:
        ws = wb[SHEET]
        header = [c.value for c in ws[1]]
        source_rows = [[c.value for c in row] for row in ws.iter_rows()]

        new_rows = []
        for stem, experiment in RUNS:
            plex, channels = maps[stem]
            row = build_row(header, experiment, plex, channels)
            ws.append(row)
            new_rows.append(row)
            named = [v for v in row[FIRST_TAG_COL - 1 :] if v is not None]
            print(f"  {experiment}: {plex}-plex, {len(named)} channels")
            print(f"    {', '.join(named)}")

        # openpyxl does not grow a table's ref on append, and Excel calls a workbook whose
        # table ref falls short of its data corrupt.
        n = ws.max_row
        ws.tables[TABLE].ref = f"A1:V{n}"

        wb.save(OUT)
    finally:
        wb.close()

    print(f"wrote {OUT}")
    verify(source_rows, new_rows)


if __name__ == "__main__":
    main()
