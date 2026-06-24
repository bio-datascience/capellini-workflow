"""Taxonomy helpers: NCBI name lookup, index sanitization, bacteria taxonomy cleaning."""

# NOTE: This file is duplicated across packages by design (to avoid a shared-utils dependency).
# Sister copies live at:
#   - progenomes_harmonizer/taxonomy.py   (full copy including NCBI lookup functions)
#   - capellini_workflow/taxonomy.py       (this file — subset: sanitize_*, clean_*, load_bacteria_taxonomy)
# When fixing a bug or adding a feature here, update all copies.

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from typing import Any

import pandas as pd

RANKS_FOR_TAXID: list[str] = ["Genus", "Family", "Order", "Class", "Phylum", "Kingdom"]

CUSTOM_RENAMES: dict[str, str] = {
    "Clostridium sensu stricto": "Clostridium sensu stricto 1",
}


def build_name_to_ncbi(names_dmp_path: str) -> dict[str, int]:
    """Parse an NCBI names.dmp file into a scientific-name → taxid mapping.

    Args:
        names_dmp_path: Path to names.dmp extracted from taxdmp.zip.

    Returns:
        Dictionary mapping scientific name strings to integer NCBI taxids.
    """
    name_to_taxid: dict[str, int] = {}
    with open(names_dmp_path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            if parts[3].strip() != "scientific name":
                continue
            try:
                taxid = int(parts[0].strip())
            except ValueError:
                continue
            name_to_taxid[parts[1].strip()] = taxid
    return name_to_taxid


def lookup_ncbi_taxid(
    row: Any,
    name_to_taxid: dict[str, int],
    ranks: list[str] = RANKS_FOR_TAXID,
) -> tuple[Any, Any]:
    """Look up an NCBI taxid for a taxonomy row, trying ranks from finest to coarsest.

    Args:
        row: A dict-like row with taxonomy rank keys.
        name_to_taxid: Mapping from scientific name to taxid.
        ranks: Ordered list of ranks to try (finest first).

    Returns:
        Tuple of (taxid, matched_rank) or (pd.NA, None) if no match found.
    """
    for rank in ranks:
        val = row.get(rank, None)
        if pd.isna(val) or str(val).strip() == "":
            continue
        taxid = name_to_taxid.get(str(val).strip())
        if taxid is not None:
            return taxid, rank
    return pd.NA, None


def assign_ncbi_taxids(
    taxonomy_table: pd.DataFrame,
    name_to_ncbi: dict[str, int],
) -> pd.DataFrame:
    """Assign NCBI taxids to every row of a taxonomy table and print a summary.

    Args:
        taxonomy_table: DataFrame with ranks as columns.
        name_to_ncbi: Scientific-name → taxid mapping from build_name_to_ncbi.

    Returns:
        DataFrame with added NCBI_taxid and taxid_matched_rank columns.
    """
    ncbi_taxids = []
    matched_ranks = []
    for _, row in taxonomy_table.iterrows():
        taxid, rank = lookup_ncbi_taxid(row, name_to_ncbi)
        ncbi_taxids.append(taxid)
        matched_ranks.append(rank)

    taxonomy_table = taxonomy_table.copy()
    taxonomy_table["NCBI_taxid"] = pd.array(ncbi_taxids, dtype="Int64")
    taxonomy_table["taxid_matched_rank"] = matched_ranks

    n_total = len(taxonomy_table)
    n_matched = taxonomy_table["NCBI_taxid"].notna().sum()
    print(f"\nNCBI taxids assigned: {n_matched} / {n_total} ({100*n_matched/n_total:.1f}%)")
    rank_counts = taxonomy_table["taxid_matched_rank"].value_counts(dropna=True)
    for rank, count in rank_counts.items():
        print(f"  {rank:<12}: {count:>6}  ({100*count/n_total:.1f}%)")
    return taxonomy_table


def build_rank_to_taxids(df_all_ncbis: pd.DataFrame, rank_col: str) -> dict[str, set]:
    """Build a mapping from rank name to the set of ProGenomes taxids in that rank.

    Args:
        df_all_ncbis: DataFrame with taxid and rank columns.
        rank_col: Column name for the rank (e.g., "genus", "family").

    Returns:
        Dictionary mapping rank name → set of integer taxids.
    """
    from collections import defaultdict

    m: dict[str, set] = defaultdict(set)
    sub = df_all_ncbis[["taxid", rank_col]].dropna()
    for taxid, name in sub.itertuples(index=False):
        if isinstance(name, str) and name.strip():
            m[name].add(int(taxid))
    return dict(m)


_SANITIZE_RE_SPACES = re.compile(r"\s+")
_SANITIZE_RE_X_PREFIX = re.compile(r"^X\.")
_SANITIZE_RE_DOTS = re.compile(r"\.+")


def sanitize_taxon_name(s: Any) -> str:
    """Normalize a taxon string consistently across studies.

    Strips, collapses whitespace, removes brackets/quotes, drops the R
    ``X.`` prefix, removes spaces/underscores, and collapses dot runs.
    """
    if s is None:
        return s  # type: ignore[return-value]
    s = str(s).strip()
    s = _SANITIZE_RE_SPACES.sub(" ", s)
    s = s.replace("[", "").replace("]", "")
    s = s.replace("'", "").replace('"', "")
    s = _SANITIZE_RE_X_PREFIX.sub("", s)
    s = s.replace(" ", "").replace("_", "")
    s = _SANITIZE_RE_DOTS.sub(".", s)
    return s


def sanitize_index(idx: Iterable) -> pd.Index:
    """Apply :func:`sanitize_taxon_name` to every element of an index."""
    return pd.Index([sanitize_taxon_name(x) for x in idx])


def clean_index_ids(idx: Iterable) -> list:
    """Strip trailing .0 float artefacts from string-cast integer IDs.

    Args:
        idx: Iterable of index values.

    Returns:
        List of cleaned string values.
    """
    out = []
    for x in idx:
        s = str(x).strip()
        if s.endswith(".0") and s[:-2].isdigit():
            s = s[:-2]
        out.append(s)
    return out


def clean_df_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Apply clean_index_ids to both the row index and column index of a DataFrame.

    Args:
        df: Input DataFrame.

    Returns:
        Copy of df with cleaned index and columns.
    """
    df = df.copy()
    df.index = clean_index_ids(df.index)
    df.columns = clean_index_ids(df.columns)
    return df


def parse_bool_series(s: pd.Series) -> pd.Series:
    """Robustly parse a boolean metadata column that may be stored as strings.

    Args:
        s: Series of bool or string values.

    Returns:
        Boolean Series.
    """
    if s.dtype == bool:
        return s
    return (
        s.astype(str)
        .str.strip()
        .str.lower()
        .map({"true": True, "1": True, "yes": True, "y": True,
              "false": False, "0": False, "no": False, "n": False})
        .fillna(False)
        .astype(bool)
    )


def load_bacteria_taxonomy(path: str) -> pd.DataFrame:
    """Load a bacteria taxonomy CSV, handling the old notebook's double-index convention.

    Args:
        path: Path to the taxonomy CSV file.

    Returns:
        DataFrame with ASV/OTU index.
    """
    import pandas as pd

    tax = pd.read_csv(path, index_col=0)
    if "Unnamed: 0" in tax.columns:
        tax = tax.set_index("Unnamed: 0")
    return tax


def clean_bacteria_taxonomy(
    tax: pd.DataFrame,
    cols_to_clean: Sequence[str] = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
    keep_cols: Sequence[str] = ("target_taxids",),
) -> pd.DataFrame:
    """Sanitize bacteria taxonomy columns and index.

    Args:
        tax: Taxonomy DataFrame.
        cols_to_clean: Rank columns to sanitize.
        keep_cols: Columns to copy through without sanitization.

    Returns:
        Cleaned taxonomy DataFrame.
    """
    out = tax.copy()
    available = [c for c in cols_to_clean if c in out.columns]
    out_clean = out.loc[:, available].apply(lambda col: col.astype(str).map(sanitize_taxon_name))
    for c in keep_cols:
        if c in out.columns:
            out_clean[c] = out[c]
    out_clean.index = sanitize_index(out_clean.index)
    return out_clean


def apply_custom_renames(df: pd.DataFrame, renames: dict[str, str] = CUSTOM_RENAMES) -> pd.DataFrame:
    """Apply a fixed set of column renames to a DataFrame.

    Args:
        df: Input DataFrame.
        renames: Mapping from old column name to new name.

    Returns:
        DataFrame with renamed columns.
    """
    present = {k: v for k, v in renames.items() if k in df.columns}
    return df.rename(columns=present)


def rename_clostridium_sensu_stricto(df: pd.DataFrame) -> pd.DataFrame:
    """Rename the Clostridium sensu stricto column to include the subspecies number.

    Args:
        df: DataFrame with genus-level abundance columns.

    Returns:
        Copy of df with corrected column name.
    """
    out = df.copy()
    out.columns = [
        "Clostridium sensu stricto 1" if c == "Clostridium sensu stricto" else c
        for c in out.columns
    ]
    return out
