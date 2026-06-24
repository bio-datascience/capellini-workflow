"""I/O helpers: file reading, writing, downloading, subprocess execution."""

# NOTE: This file is duplicated across packages by design (to avoid a shared-utils dependency).
# Sister copies live at:
#   - progenomes_harmonizer/io.py   (full copy)
#   - capellini_workflow/io.py       (this file — full copy)
#   - procs_maker/io.py              (subset: download_if_missing, sh, atomic helpers)
# When fixing a bug or adding a feature here, update all copies.

from __future__ import annotations

import gzip
import logging
import os
import shutil
import subprocess
import urllib.request
from contextlib import contextmanager
from pathlib import Path
from typing import IO

logger = logging.getLogger(__name__)


# ── Atomic-write primitives ────────────────────────────────────────────────
#
# Snakemake only checks file *existence* — it doesn't validate content. A
# crash mid-write of a non-declared (side-effect) file would leave a partial
# file on disk that the next run silently reuses. These helpers guarantee
# that ``final_path`` either does not exist or contains the complete output:
# we always write to ``final_path + ".part"`` and ``os.replace()`` to the
# canonical name only on success.

def _part_path(final_path: Path) -> Path:
    return final_path.with_name(final_path.name + ".part")


@contextmanager
def atomic_write(final_path: str | Path, mode: str = "w", **open_kwargs):
    """Context manager: write to ``<final>.part`` then atomically rename."""
    final_path = Path(final_path)
    final_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = _part_path(final_path)
    fh = open(tmp_path, mode, **open_kwargs)
    try:
        yield fh
        fh.flush()
        try:
            os.fsync(fh.fileno())
        except (OSError, AttributeError):
            pass
    except BaseException:
        try:
            fh.close()
        except Exception:
            pass
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass
        raise
    else:
        fh.close()
        os.replace(tmp_path, final_path)


def atomic_download(url: str, dest: str | Path, *, label: str | None = None) -> Path:
    """Download ``url`` to ``dest`` via a ``.part`` file, rename on success."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = _part_path(dest)
    name = label or dest.name
    print(f" • Download of {name} - this could take a while")
    try:
        urllib.request.urlretrieve(url, str(tmp_path))
    except BaseException:
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass
        raise
    os.replace(tmp_path, dest)
    return dest


def download_if_missing(url: str, dest: str | Path, *, label: str | None = None) -> Path:
    """Download ``url`` to ``dest`` if ``dest`` is missing.

    Uses ``atomic_download`` so a Ctrl-C'd download never leaves a partial
    file masquerading as a complete one. Also auto-cleans any stale
    ``.part`` from a prior crashed run.
    """
    dest = Path(dest)
    name = label or dest.name
    if dest.exists():
        print(f"{name} found - skipping download")
        return dest
    if not url:
        raise ValueError(
            f"Cannot download {name}: no URL configured "
            f"(target path: {dest})"
        )
    stale = _part_path(dest)
    if stale.exists():
        print(f" • Removing stale partial download: {stale.name}")
        stale.unlink()
    return atomic_download(url, dest, label=label)


def open_maybe_gzip(path: str | Path, mode: str = "rt", encoding: str = "utf-8") -> IO:
    """Open a file, transparently decompressing if it ends with .gz.

    Args:
        path: Path to the file.
        mode: File open mode.
        encoding: Text encoding.

    Returns:
        File-like object.
    """
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, mode, encoding=encoding, errors="replace")
    return open(p, mode, encoding=encoding, errors="replace")


def read_table(path: str | Path, index_col: int = 0, **kwargs):
    """Read a CSV, TSV, or Excel file into a DataFrame, including gzip-compressed variants.

    Args:
        path: Path to the tabular file.
        index_col: Column to use as the row index.
        **kwargs: Additional keyword arguments forwarded to pandas.

    Returns:
        pd.DataFrame with the file contents.
    """
    import pandas as pd

    path = Path(path)
    suffixes = "".join(path.suffixes).lower()

    if suffixes.endswith(".xlsx") or suffixes.endswith(".xls"):
        return pd.read_excel(path, index_col=index_col, **kwargs)
    if suffixes.endswith(".tsv") or suffixes.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t", index_col=index_col, **kwargs)
    if suffixes.endswith(".csv") or suffixes.endswith(".csv.gz") or suffixes.endswith(".gz"):
        return pd.read_csv(path, index_col=index_col, **kwargs)
    raise ValueError(f"Unsupported file type: {path}")


def write_df(df, path: str | Path, *, overwrite: bool = True, verbose: bool = False, **to_csv_kwargs) -> Path:
    """Write a DataFrame to CSV, creating parent directories as needed.

    Args:
        df: DataFrame to write.
        path: Destination file path.
        overwrite: Skip if file exists and overwrite is False.
        verbose: Print a message on skip or save.
        **to_csv_kwargs: Forwarded to DataFrame.to_csv.

    Returns:
        Path to the written file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        if verbose:
            print("skip existing:", path)
        return path
    df.to_csv(path, **to_csv_kwargs)
    if verbose:
        print("saved:", path, getattr(df, "shape", ""))
    return path


def sh(cmd: str, desc: str = "") -> subprocess.CompletedProcess:
    """Run a shell command, printing description and command, raising on failure.

    Args:
        cmd: Shell command string.
        desc: Human-readable description printed before execution.

    Returns:
        CompletedProcess result.

    Raises:
        RuntimeError: If the command exits with a non-zero return code,
            including full stdout and stderr in the message.
    """
    if desc:
        print(desc)
    print(f"Executing command: {cmd}")
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        msg = (
            f"Command failed with code {e.returncode}\n"
            f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
        )
        raise RuntimeError(msg) from e
    if r.stdout:
        print(r.stdout)
    if r.stderr.strip():
        print("STDERR (tool messages):\n" + r.stderr)
    return r
