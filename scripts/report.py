#!/usr/bin/env python3
"""CAPELLINI post-run report: per-stage statistics + figures.

Runs after a pipeline execution and, for **every** stage, reports whether it was
performed (its outputs exist) or skipped, together with summary statistics. It
then reproduces the figures from the original analysis notebook
(``1_IBD/SpacePHARER_IBD.ipynb``): degree distributions, spacer support,
p-value distribution, top-20 host/virus bars, the bipartite host-phage network,
a presence/absence heatmap, and the per-genome protein-count distribution.

Design goals:
  * self-contained — reads only the run's output files (no in-memory pipeline
    state), so it works after a full run, a partial run, or a resumed run;
  * v4-aware — proGenomes4 ids carry the versioned GCA, and taxids are resolved
    through the GCA->taxid map (not parsed from the header as in the v3 notebook);
  * never fatal — a missing input for any single stat/figure is reported as
    "skipped", not raised, so the report always completes.

Usage:
    python report.py --configfile config/ibd.yaml
    # or, programmatically:
    from report import generate_report; generate_report(cfg_dict)
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless: figures are saved, never shown
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

# ── styling (mirrors the notebook) ───────────────────────────────────────────
PALETTE = {"host": "#E07B39", "virus": "#4A90C4", "edge": "#AAAAAA"}
HIGH_DEGREE_PERCENTILE = 0.95
_RULE = "═" * 68


def _apply_style() -> None:
    plt.rcParams.update({
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight",
    })


# ── small helpers ────────────────────────────────────────────────────────────

def _fmt(n) -> str:
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return str(n)


def _fasta_count(path: Path) -> int:
    """Count records in a (optionally gzipped) FASTA by scanning for '>' lines."""
    op = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def _count_data_lines(path: Path) -> int:
    """Fast newline count minus a header line (for CSV row counts)."""
    op = gzip.open if str(path).endswith(".gz") else open
    total = 0
    with op(path, "rb") as fh:
        while True:
            buf = fh.read(1 << 20)
            if not buf:
                break
            total += buf.count(b"\n")
    return max(0, total - 1)


def _csv_dims(path: Path) -> tuple[int, int]:
    """(rows, cols) of a CSV without loading it — header for cols, stream for rows."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        header = fh.readline()
    cols = header.count(",")  # data cols = commas (index col is the +0th field)
    return _count_data_lines(path), cols


def _first_gca(text: str) -> str | None:
    m = re.search(r"(GC[AF]_\d+\.\d+)", str(text))
    return m.group(1) if m else None


# ── report printer: tees to stdout and to a saved summary file ───────────────

class _Report:
    def __init__(self) -> None:
        self._buf = io.StringIO()

    def __call__(self, msg: str = "") -> None:
        print(msg)
        self._buf.write(msg + "\n")

    def status(self, performed: bool, label: str, note: str = "") -> None:
        mark = "✓ performed" if performed else "○ skipped  "
        line = f"  [{mark}]  {label}"
        if note:
            line += f"   ({note})"
        self(line)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self._buf.getvalue())


# ── GCA -> taxid (lazy, v4) ──────────────────────────────────────────────────

def _load_gca_to_taxid(download_path: Path, cfg: dict) -> dict:
    """Best-effort load of the proGenomes4 GCA->taxid map (empty dict on failure)."""
    candidates = []
    if cfg.get("pg4_ncbi_taxonomy"):
        candidates.append(Path(cfg["pg4_ncbi_taxonomy"]))
    candidates += [
        download_path / "pg4_ncbi_taxonomy.tsv.gz",
        download_path / "pg4_ncbi_taxonomy.tsv",
    ]
    path = next((p for p in candidates if p and p.exists()), None)
    if path is None:
        return {}
    try:
        from progenomes_harmonizer.reference import load_gca_to_taxid
        return load_gca_to_taxid(path)
    except Exception:
        return {}


# ═════════════════════════════════════════════════════════════════════════════
# STAGE STATISTICS
# ═════════════════════════════════════════════════════════════════════════════

def _stage_dada2(p: _Report, dada2: Path, suffix: str) -> None:
    otu = dada2 / f"OTU_table_{suffix}.csv"
    tax = dada2 / f"taxonomy_table_{suffix}.csv"
    asv = dada2 / f"ASV_sequences_{suffix}.fasta"
    performed = otu.exists() and tax.exists()
    p.status(performed, "DADA2 — 16S amplicon processing")
    if not performed:
        return
    n_samples, _ = _csv_dims(otu)
    n_asvs = _count_data_lines(tax)
    n_asv_fa = _fasta_count(asv) if asv.exists() else None
    p(f"        samples: {_fmt(n_samples)}   ASVs: {_fmt(n_asvs)}"
      + (f"   (ASV FASTA: {_fmt(n_asv_fa)})" if n_asv_fa is not None else ""))


def _stage_harmonize(p: _Report, dada2: Path, mmseq: Path, suffix: str,
                     ref16s: Path | None, gca2tax: dict, cfg: dict) -> dict:
    """Returns the target taxid set (used later by SpacePHARER stats)."""
    withids = dada2 / f"taxonomy_table_{suffix}_withIDs.csv"
    performed = withids.exists()
    p.status(performed, "Harmonize — NCBI / GCA assignment")
    out: dict = {"target_taxids": set()}
    if not performed:
        return out

    tax = pd.read_csv(withids, index_col=0, low_memory=False)
    n = len(tax)
    sp = tax["progenomes_taxid_species"].notna().sum() if "progenomes_taxid_species" in tax else 0
    ge = tax["progenomes_taxid_genus"].notna().sum() if "progenomes_taxid_genus" in tax else 0
    fa = tax["progenomes_taxid_family"].notna().sum() if "progenomes_taxid_family" in tax else 0
    cols = [c for c in ("progenomes_taxid_species", "progenomes_taxid_genus",
                        "progenomes_taxid_family") if c in tax]
    anya = tax[cols].notna().any(axis=1).sum() if cols else 0
    p(f"        ASVs assigned (any layer): {_fmt(anya)} / {_fmt(n)} ({100*anya/max(n,1):.1f}%)")
    p(f"          species {_fmt(sp)} ({100*sp/max(n,1):.1f}%) · "
      f"genus {_fmt(ge)} ({100*ge/max(n,1):.1f}%) · "
      f"family {_fmt(fa)} ({100*fa/max(n,1):.1f}%)")

    # Target taxid set used to scope the SpacePHARER spacers. Mirrors
    # _target_taxid_set() in spacepharer_wrapper.py: species_level only ADDS the
    # species layer on top of genus, it never replaces it (dropping genus would
    # coarsen the set, since addSpecies matches only exactly).
    cols = ["progenomes_taxid_genus", "progenomes_taxid_family"]
    if cfg.get("species_level"):
        cols.append("progenomes_taxid_species")
    target = set()
    for col in cols:
        if col in tax:
            target |= set(tax[col].dropna().astype(int))
    out["target_taxids"] = target
    p(f"        unique NCBI target IDs: {_fmt(len(target))}")

    # MMseqs2 match quality from output.m8 (top hit per query)
    m8 = mmseq / "output.m8"
    if m8.exists():
        try:
            df = pd.read_csv(m8, sep="\t", header=None, usecols=[0, 1, 2, 11],
                             names=["query", "target", "pident", "bitscore"],
                             dtype={"query": str, "target": str})
            top = df.loc[df.groupby("query")["bitscore"].idxmax()]
            min_bs = cfg.get("min_bitscore", 50)
            matched = top[top["bitscore"] >= min_bs]
            scale = 1.0 if top["pident"].max() > 1.0 else 100.0
            p(f"        MMseqs2 16S hits (bitscore ≥ {min_bs}): "
              f"{_fmt(matched['query'].nunique())} ASVs "
              f"({100*matched['query'].nunique()/max(n,1):.1f}%)")
            p(f"          bitscore med/mean {top['bitscore'].median():.0f}/"
              f"{top['bitscore'].mean():.0f} · "
              f"identity med/mean {top['pident'].median()*scale:.1f}%/"
              f"{top['pident'].mean()*scale:.1f}%")
            if gca2tax and target:
                tgt_str = target
                top = top.assign(gca=top["target"].map(_first_gca))
                top = top.assign(taxid=top["gca"].map(gca2tax))
                n_gca = top[top["taxid"].isin(tgt_str)]["gca"].nunique()
                p(f"        bacterial genomes in target set (unique GCAs): {_fmt(n_gca)}")
        except Exception as exc:  # noqa: BLE001
            p(f"        [MMseqs2 match-quality stats unavailable: {exc}]")

    # 16S reference used
    if ref16s and ref16s.exists():
        n_rec = _fasta_count(ref16s)
        tsv = ref16s.with_suffix(".tsv")
        n_gen = None
        if tsv.exists():
            genomes = set()
            with open(tsv) as fh:
                next(fh, None)
                for line in fh:
                    genomes.add(line.split("\t", 1)[0])
            n_gen = len(genomes)
        p(f"        16S reference: {_fmt(n_rec)} copies"
          + (f" from {_fmt(n_gen)} genomes" if n_gen is not None else "")
          + f"   [{ref16s.name}]")
    return out


def _load_spacepharer(pred_tsv: Path, gca2tax: dict) -> dict:
    """Parse phage_host_predictions.tsv into the pairs/degree structures."""
    cols = ["spacer_id", "phage_id", "pvalue", "aln_len", "mismatch",
            "qstart", "qend", "qprot_aln", "sprot_aln"]
    sp = pd.read_csv(pred_tsv, sep="\t", comment="#", header=None, names=cols,
                     engine="python", dtype=str, on_bad_lines="skip"
                     ).dropna(how="all").reset_index(drop=True)
    for c in ("pvalue", "aln_len", "mismatch", "qstart", "qend"):
        sp[c] = pd.to_numeric(sp[c], errors="coerce")
    # SpacePHARER '>' hit lines keep the '>' on the spacer id — drop it so host
    # ids render cleanly (the taxid is parsed from the GCA regardless).
    sp["spacer_id"] = sp["spacer_id"].str.replace(r"^>", "", regex=True)
    sp["host_id"] = sp["spacer_id"].str.replace(r"_spacer_\d+$", "", regex=True)
    sp["gca"] = sp["spacer_id"].map(_first_gca)
    sp["ncbi_taxid"] = sp["gca"].map(gca2tax) if gca2tax else None

    pairs = sp[["host_id", "phage_id"]].drop_duplicates().reset_index(drop=True)
    pairs = pairs.merge(
        sp.groupby(["host_id", "phage_id"])["spacer_id"].nunique()
          .rename("n_spacers").reset_index(), on=["host_id", "phage_id"])
    pairs = pairs.merge(
        sp.groupby(["host_id", "phage_id"])["pvalue"].min()
          .rename("best_pval").reset_index(), on=["host_id", "phage_id"])
    host_deg = pairs.groupby("host_id")["phage_id"].nunique().rename("degree")
    virus_deg = pairs.groupby("phage_id")["host_id"].nunique().rename("degree")

    # connected components over the bipartite graph
    edges = pairs[["host_id", "phage_id"]].dropna().drop_duplicates()
    h2p, p2h = defaultdict(set), defaultdict(set)
    for h, ph in edges.values:
        h2p[h].add(ph); p2h[ph].add(h)
    seen_h, seen_p, comps = set(), set(), []
    for h0 in edges["host_id"].unique():
        if h0 in seen_h:
            continue
        stack, ch, cp = [h0], set(), set()
        while stack:
            h = stack.pop()
            if h in seen_h:
                continue
            seen_h.add(h); ch.add(h)
            for ph in h2p.get(h, ()):
                if ph not in seen_p:
                    seen_p.add(ph); cp.add(ph)
                    stack.extend(h2 for h2 in p2h.get(ph, ()) if h2 not in seen_h)
        comps.append((len(ch), len(cp)))
    return {"sp": sp, "pairs": pairs, "host_deg": host_deg,
            "virus_deg": virus_deg, "components": comps}


def _stage_spacepharer(p: _Report, sp_dir: Path, cfg: dict, data: dict | None) -> None:
    pred = sp_dir / "output" / "phage_host_predictions.tsv"
    performed = pred.exists() and pred.stat().st_size > 0
    p.status(performed, "SpacePHARER — CRISPR host-phage predictions")
    if not performed or data is None:
        return
    pairs, host_deg, virus_deg = data["pairs"], data["host_deg"], data["virus_deg"]
    comps = data["components"]
    n_hosts, n_vir, n_int = host_deg.shape[0], virus_deg.shape[0], pairs.shape[0]
    density = n_int / (n_hosts * n_vir) if n_hosts * n_vir else 0.0
    largest = max(comps, key=lambda x: sum(x)) if comps else (0, 0)
    n_gca = data["sp"]["gca"].nunique()
    p(f"        interactions: {_fmt(n_int)}  ·  hosts: {_fmt(n_hosts)} "
      f"(genomes: {_fmt(n_gca)})  ·  phages: {_fmt(n_vir)}")
    p(f"        spacer-hit rows: {_fmt(len(data['sp']))}  ·  density: {density:.5f}"
      f"  ·  FDR: {cfg.get('fdr', 0.05)}")
    p(f"        components: {_fmt(len(comps))}  ·  largest: "
      f"{_fmt(largest[0])} hosts / {_fmt(largest[1])} phages")
    p(f"        host degree med/mean/max: {host_deg.median():.0f}/"
      f"{host_deg.mean():.1f}/{int(host_deg.max())}  ·  "
      f"virus degree med/mean/max: {virus_deg.median():.0f}/"
      f"{virus_deg.mean():.1f}/{int(virus_deg.max())}")


def _stage_procs(p: _Report, procs: Path) -> None:
    pc = procs / "pc_matrix.csv"
    pb = procs / "pb_matrix.csv"
    performed = pc.exists() or pb.exists()
    p.status(performed, "ProCs Maker — protein-cluster matrices")
    if not performed:
        return
    bac = procs / "Targets Proteins Extraction" / "BacterialProteinsCollection.fasta"
    vir = procs / "Targets Proteins Extraction" / "ViralProteinsCollection.fasta"
    if bac.exists():
        p(f"        bacterial proteins: {_fmt(_fasta_count(bac))}")
    if vir.exists():
        p(f"        viral proteins:     {_fmt(_fasta_count(vir))}")
    clust = procs / "Clustering" / "clusterRes_cluster.tsv"
    if clust.exists():
        members = clusters = 0
        seen = set()
        with open(clust) as fh:
            for line in fh:
                members += 1
                rep = line.split("\t", 1)[0]
                if rep not in seen:
                    seen.add(rep); clusters += 1
        p(f"        protein clusters: {_fmt(clusters)}  (members: {_fmt(members)}, "
          f"mean size {members/max(clusters,1):.1f})")
    for m in (pc, pb):
        if m.exists():
            r, c = _csv_dims(m)
            p(f"        {m.name}: {_fmt(r)} × {_fmt(c)}")


def _stage_network(p: _Report, net: Path, cfg: dict) -> None:
    subs = [
        ("common_abundance",        "run_common_abundance",       net / "common" / "B_processed.csv"),
        ("shrinkage_correlations",  "run_shrinkage_correlations", net / "shrinkage" / "shrinkage_corr_BV.csv.gz"),
        ("raw_crispr_network",      "run_raw_crispr_networks",    net / "crispr_raw" / "crispr_net.csv.gz"),
        ("smooth_crispr",           "run_smooth_crispr",          net / "crispr_smooth" / "crispr_smooth_vir_bac.csv.gz"),
        ("xstar",                   "run_xstar",                  net / "xstar" / "X_star.csv.gz"),
    ]
    for name, flag, out in subs:
        performed = out.exists()
        note = "" if performed else ("disabled in config" if not cfg.get(flag, True) else "no output")
        p.status(performed, f"Network — {name}", note)
        if performed:
            try:
                r, c = _csv_dims(out)
                shape = "(empty)" if r == 0 and c == 0 else f"{_fmt(r)} × {_fmt(c)}"
                p(f"        {out.name}: {shape}")
            except Exception:  # noqa: BLE001
                pass
    common = net / "common"
    if (common / "B_processed.csv").exists() and (common / "V_processed.csv").exists():
        br, bc = _csv_dims(common / "B_processed.csv")
        vr, vc = _csv_dims(common / "V_processed.csv")
        p(f"        aligned network: {_fmt(br)} samples · "
          f"{_fmt(bc)} bacterial × {_fmt(vc)} viral features")


# ═════════════════════════════════════════════════════════════════════════════
# FIGURES  (reproduce the notebook)
# ═════════════════════════════════════════════════════════════════════════════

def _save(fig_dir: Path, name: str) -> None:
    for ext in ("png", "pdf"):
        plt.savefig(fig_dir / f"{name}.{ext}")
    plt.close()


def _figures_spacepharer(p: _Report, data: dict, fig_dir: Path) -> int:
    import seaborn as sns
    import networkx as nx
    pairs = data["pairs"]; host_deg = data["host_deg"]; virus_deg = data["virus_deg"]
    made = 0

    # 4a — degree distributions (log-binned)
    fig, axes = plt.subplots(1, 2, figsize=(6, 2.2))
    for ax, d, lab, col in [
        (axes[0], host_deg, "Host degree (phages / bacterium)", PALETTE["host"]),
        (axes[1], virus_deg, "Virus degree (bacteria / phage)", PALETTE["virus"]),
    ]:
        if len(d):
            bins = np.logspace(np.log10(max(d.min(), 1)), np.log10(max(d.max(), 2)), 30)
            ax.hist(d, bins=bins, color=col, edgecolor="white", linewidth=0.3)
            ax.set_xscale("log")
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
        ax.set_xlabel(lab); ax.set_ylabel("Count")
    fig.suptitle("Degree distributions", y=1.02)
    plt.tight_layout(); _save(fig_dir, "01_degree_distributions"); made += 1

    # 4b — spacer support
    fig, ax = plt.subplots(figsize=(2.2, 1.7))
    ax.hist(pairs["n_spacers"], bins=30, color=PALETTE["host"], edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Unique spacers per interaction"); ax.set_ylabel("Count")
    ax.set_title("Spacer support per host–phage pair")
    _save(fig_dir, "02_spacer_support"); made += 1

    # 4c — significance
    neg_log = -np.log10(pairs["best_pval"].replace(0, np.nan).dropna())
    if len(neg_log):
        fig, ax = plt.subplots(figsize=(2.2, 1.7))
        ax.hist(neg_log, bins=40, color=PALETTE["virus"], edgecolor="white", linewidth=0.3)
        ax.set_xlabel("–log10(best p-value)"); ax.set_ylabel("Count")
        ax.set_title("Interaction significance")
        _save(fig_dir, "03_pvalue_distribution"); made += 1

    # 4d — top-20 bars
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2))
    for ax, d, col, title in [
        (axes[0], host_deg.sort_values(ascending=False).head(20), PALETTE["host"], "Top 20 hosts by virus range"),
        (axes[1], virus_deg.sort_values(ascending=False).head(20), PALETTE["virus"], "Top 20 viruses by host range"),
    ]:
        ax.barh(range(len(d)), d.values, color=col)
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels([x[:30] + "…" if len(x) > 30 else x for x in d.index], fontsize=6)
        ax.invert_yaxis(); ax.set_xlabel("Degree"); ax.set_title(title)
    plt.tight_layout(); _save(fig_dir, "04_top20_degree"); made += 1

    # 4e — bipartite network (cap node count)
    MAX_NODES = 500
    if host_deg.shape[0] + virus_deg.shape[0] > MAX_NODES:
        top_h = host_deg.sort_values(ascending=False).head(MAX_NODES // 2).index
        sub = pairs[pairs["host_id"].isin(top_h)]
        subtitle = f"(top {MAX_NODES // 2} hosts by degree)"
    else:
        sub, subtitle = pairs, ""
    B = nx.Graph()
    B.add_nodes_from(sub["host_id"].unique(), bipartite=0)
    B.add_nodes_from(sub["phage_id"].unique(), bipartite=1)
    B.add_edges_from(sub[["host_id", "phage_id"]].itertuples(index=False, name=None))
    if B.number_of_nodes():
        hosts_in = [n for n, d in B.nodes(data=True) if d["bipartite"] == 0]
        virus_in = [n for n, d in B.nodes(data=True) if d["bipartite"] == 1]
        pos = nx.spring_layout(B, k=0.5, seed=42)
        plt.figure(figsize=(7, 7))
        nx.draw_networkx_nodes(B, pos, nodelist=hosts_in, node_size=25,
                               node_color=PALETTE["host"], alpha=0.6, label="Bacterial hosts")
        nx.draw_networkx_nodes(B, pos, nodelist=virus_in, node_size=25,
                               node_color=PALETTE["virus"], alpha=0.6, label="Viruses")
        nx.draw_networkx_edges(B, pos, width=0.6, alpha=0.35)
        plt.legend(scatterpoints=1); plt.axis("off")
        plt.title(f"Host–phage interaction network {subtitle}")
        _save(fig_dir, "05_bipartite_network"); made += 1

    # 4f — presence/absence heatmap
    pivot = (pairs.groupby(["host_id", "phage_id"])["n_spacers"].sum()
                  .unstack(fill_value=0).astype(bool).astype(int))
    pivot = pivot.loc[pivot.sum(axis=1) >= 3, pivot.sum(axis=0) >= 3]
    if pivot.shape[0] and pivot.shape[1]:
        if pivot.shape[0] > 80:
            pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).head(80).index]
        if pivot.shape[1] > 80:
            pivot = pivot[pivot.sum(axis=0).sort_values(ascending=False).head(80).index]
        fig, ax = plt.subplots(figsize=(max(2, pivot.shape[1] * 0.06),
                                        max(2, pivot.shape[0] * 0.06)))
        sns.heatmap(pivot, ax=ax, cmap="Oranges", xticklabels=False, yticklabels=False,
                    cbar_kws={"label": "Interaction (0/1)", "shrink": 0.5})
        ax.set_xlabel("Viruses"); ax.set_ylabel("Bacterial hosts")
        ax.set_title(f"Host–phage presence/absence ({pivot.shape[0]}×{pivot.shape[1]})")
        _save(fig_dir, "06_presence_absence_heatmap"); made += 1
    return made


def _figure_protein_counts(procs: Path, fig_dir: Path) -> bool:
    bac = procs / "Targets Proteins Extraction" / "BacterialProteinsCollection.fasta"
    if not bac.exists():
        return False
    counter: Counter = Counter()
    with open(bac) as fh:
        for line in fh:
            if line.startswith(">"):
                g = _first_gca(line)
                if g:
                    counter[g] += 1
    if not counter:
        return False
    values = np.array(sorted(counter.values(), reverse=True))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(range(len(values)), values, color="#4F75FF", width=1.0)
    mean_v = float(values.mean())
    ax.axhline(mean_v, color="red", linestyle="--", linewidth=1.5, label=f"Mean: {mean_v:.1f}")
    ax.set_xticks([]); ax.set_xlabel("Genomes"); ax.set_ylabel("Protein counts")
    ax.set_title("Distribution of protein counts among genomes")
    ax.legend()
    plt.tight_layout(); _save(fig_dir, "07_protein_count_distribution")
    return True


def _high_degree_hosts(p: _Report, data: dict) -> None:
    host_deg = data["host_deg"]; sp = data["sp"]
    if host_deg.empty:
        return
    df = host_deg.reset_index().rename(columns={"degree": "host_degree"})
    taxid_map = (sp[["host_id", "ncbi_taxid"]].dropna().drop_duplicates()
                 .groupby("host_id")["ncbi_taxid"].first()) if "ncbi_taxid" in sp else {}
    thr = df["host_degree"].quantile(HIGH_DEGREE_PERCENTILE)
    high = df[df["host_degree"] >= thr].sort_values("host_degree", ascending=False)
    p("")
    p(f"  HIGH-DEGREE HOSTS  (top {int((1-HIGH_DEGREE_PERCENTILE)*100)}% — degree ≥ {thr:.0f})")
    p(f"    hosts in network: {_fmt(df.shape[0])}   high-degree: {_fmt(high.shape[0])}")
    for _, row in high.head(15).iterrows():
        hid = row["host_id"]
        tid = taxid_map.get(hid) if hasattr(taxid_map, "get") else None
        p(f"    · {hid[:48]:<48}  degree {int(row['host_degree']):>4}"
          + (f"  taxid {int(tid)}" if pd.notna(tid) else ""))


# ═════════════════════════════════════════════════════════════════════════════
# ORCHESTRATION
# ═════════════════════════════════════════════════════════════════════════════

def generate_report(cfg: dict) -> Path | None:
    """Print per-stage statistics and write figures for a completed run.

    Returns the figures directory (or None if the run has no ``base``).
    """
    _apply_style()
    base = cfg.get("base")
    if not base:
        print("[report] no 'base' in config — nothing to report")
        return None
    base = Path(base)
    suffix = {"forward": "F", "reverse": "R", "paired": "P"}.get(cfg.get("direction", "forward"), "F")
    download_path = Path(cfg.get("download_path", "") or base)
    dada2 = base / "DADA2 output"
    mmseq = base / "MMSeqs2 Output"
    sp_dir = base / "SpacePHARER output"
    procs = base / "Procs Estimations"
    net = base / "Enhanced Networks"
    report_dir = base / "report"
    fig_dir = report_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # resolve the 16S reference (per-run copy, else shared cache)
    ref16s = next((c for c in (dada2 / "pg4_16s.fasta",
                               download_path / "pg4_16s.fasta") if c.exists()), None)

    p = _Report()
    import datetime
    p(_RULE)
    p("  CAPELLINI — RUN REPORT")
    p(f"  base : {base}")
    p(f"  time : {datetime.datetime.now():%Y-%m-%d %H:%M}")
    p(_RULE)

    p("")
    p("STAGES")
    gca2tax = _load_gca_to_taxid(download_path, cfg)
    _stage_dada2(p, dada2, suffix)
    harmonize = _stage_harmonize(p, dada2, mmseq, suffix, ref16s, gca2tax, cfg)

    # SpacePHARER: load once, reuse for stats + figures
    pred = sp_dir / "output" / "phage_host_predictions.tsv"
    sp_data = None
    if pred.exists() and pred.stat().st_size > 0:
        try:
            sp_data = _load_spacepharer(pred, gca2tax)
        except Exception as exc:  # noqa: BLE001
            p(f"  [SpacePHARER parse failed: {exc}]")
    _stage_spacepharer(p, sp_dir, cfg, sp_data)
    _stage_procs(p, procs)
    _stage_network(p, net, cfg)

    # ── figures ──────────────────────────────────────────────────────────────
    p("")
    p("FIGURES")
    n_fig = 0
    if sp_data is not None:
        try:
            n_fig += _figures_spacepharer(p, sp_data, fig_dir)
        except Exception as exc:  # noqa: BLE001
            p(f"  [SpacePHARER figures skipped: {exc}]")
    try:
        if _figure_protein_counts(procs, fig_dir):
            n_fig += 1
    except Exception as exc:  # noqa: BLE001
        p(f"  [protein-count figure skipped: {exc}]")
    p(f"  {n_fig} figure(s) written to: {fig_dir}")

    # ── high-degree hosts ─────────────────────────────────────────────────────
    if sp_data is not None:
        try:
            _high_degree_hosts(p, sp_data)
        except Exception as exc:  # noqa: BLE001
            p(f"  [high-degree host analysis skipped: {exc}]")

    p(_RULE)
    p.save(report_dir / "summary.txt")
    print(f"[report] summary saved to: {report_dir / 'summary.txt'}")
    return fig_dir


def _load_yaml(path: Path) -> dict:
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="CAPELLINI post-run report (stats + figures)")
    ap.add_argument("--configfile", required=True, help="Path to the run's YAML config")
    args = ap.parse_args(argv)
    cfg = _load_yaml(Path(args.configfile))
    generate_report(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
