"""Download the bundled spacers_CompleteCollection.fasta from the GitHub release.

This file is too large to ship inside the source repo, so it is hosted as
a release asset on GitHub and downloaded on demand.

Override the source tag with the ``CAPELLINI_REFERENCES_TAG`` env var if you
need to pin a different release.
"""

from __future__ import annotations

import importlib.resources as pkg_resources
import logging
import os
import shutil
import sys
import urllib.request
from pathlib import Path
from urllib.error import HTTPError, URLError

logger = logging.getLogger(__name__)

DEFAULT_TAG = "v0.1.0"
RELEASES_BASE = "https://github.com/AlexDellOrti/Capellini/releases/download"

# (asset_filename, destination relative to capellini_workflow/data/)
ASSETS: tuple[tuple[str, str], ...] = (
    ("spacers_CompleteCollection.fasta", "references/spacers/spacers_CompleteCollection.fasta"),
)


def _data_root() -> Path:
    """Return the absolute path to the installed ``capellini_workflow/data`` directory."""
    ref = pkg_resources.files("capellini_workflow").joinpath("data")
    return Path(str(ref))


def _release_tag() -> str:
    return os.environ.get("CAPELLINI_REFERENCES_TAG", DEFAULT_TAG).strip() or DEFAULT_TAG


def _asset_url(tag: str, filename: str) -> str:
    return f"{RELEASES_BASE}/{tag}/{filename}"


def _human_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def _download_one(url: str, dest: Path) -> None:
    """Stream-download ``url`` to ``dest`` with a simple progress bar."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")

    def _progress(block_num: int, block_size: int, total: int) -> None:
        if total <= 0:
            return
        downloaded = block_num * block_size
        pct = min(100.0, downloaded * 100.0 / total)
        bar = "█" * int(pct // 2) + "·" * (50 - int(pct // 2))
        sys.stdout.write(
            f"\r  [{bar}] {pct:5.1f}%  {_human_size(downloaded)} / {_human_size(total)}"
        )
        sys.stdout.flush()

    print(f"  -> {url}")
    try:
        urllib.request.urlretrieve(url, str(tmp), reporthook=_progress)
    except (HTTPError, URLError) as exc:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
        raise RuntimeError(f"Download failed for {url}: {exc}") from exc
    print()  # newline after progress bar
    shutil.move(str(tmp), str(dest))


def fetch_references(
    tag: str | None = None,
    overwrite: bool = False,
) -> list[Path]:
    """Download every bundled reference the workflow needs.

    This fetches:
      * ``spacers_CompleteCollection.fasta`` into ``capellini_workflow/data/``
        (owned by this package).
      * ``progenome16S.fasta`` into ``progenomes_harmonizer/data/`` by
        delegating to that tool's own ``fetch_references``. Each tool owns
        its own bundle, but from the user's perspective one "Fetch
        references" action should populate both — that's what this wrapper
        guarantees.

    Args:
        tag: Release tag to pull from. Defaults to env var or DEFAULT_TAG.
        overwrite: If True, re-download even when files already exist.

    Returns:
        List of paths to downloaded (or already-present) files, from both
        packages.
    """
    tag = tag or _release_tag()
    data_root = _data_root()

    print(f"capellini-workflow references - release tag: {tag}")
    print(f"Target directory:  {data_root}\n")

    paths: list[Path] = []
    for filename, rel in ASSETS:
        dest = data_root / rel
        if dest.exists() and not overwrite:
            print(f"  {rel}  (already present, skipping)")
            paths.append(dest)
            continue
        print(f"  {rel}")
        _download_one(_asset_url(tag, filename), dest)
        size = dest.stat().st_size
        print(f"  saved: {dest}  ({_human_size(size)})\n")
        paths.append(dest)

    # Chain the harmonizer's fetch so progenome16S.fasta also ends up bundled.
    # We import lazily so the workflow can still be installed without the
    # harmonizer present (e.g. for development on the workflow alone).
    try:
        from progenomes_harmonizer.fetch_references import (
            fetch_references as _fetch_harmonizer_refs,
        )
    except ImportError:
        print(
            "\n[warning] progenomes-harmonizer is not installed in this "
            "environment — skipping progenome16S.fasta fetch.\n"
            "Install it with `pip install progenomes-harmonizer` and re-run "
            "fetch-references to get the 16S bundle."
        )
    else:
        print("\n--- delegating to progenomes-harmonizer fetch-references ---")
        try:
            paths.extend(_fetch_harmonizer_refs(tag=tag, overwrite=overwrite))
        except RuntimeError as exc:
            # Let the caller decide how to surface this; print and continue
            # so a transient 16S failure doesn't roll back the spacers.
            print(f"[warning] harmonizer fetch failed: {exc}")

    print("Done.")
    return paths
