"""Interactive terminal UI for CAPELLINI workflow.

Provides a Rich + questionary interactive menu (``capellini`` with no args)
as well as direct Snakemake dispatch (``capellini run --configfile …``).
"""

from __future__ import annotations

import importlib.resources as pkg_resources
import os
import shutil
import subprocess
import sys
import threading
from collections import deque
from pathlib import Path

import questionary
import yaml
from rich.console import Console, Group
from rich.panel import Panel
from rich.live import Live
from rich.table import Table
from rich.text import Text

CONSOLE = Console()
LAST_CONFIG_FILE = Path.home() / ".capellini" / "last_config"

# Locate the Snakefile relative to this package
_WORKFLOW_DIR = Path(__file__).resolve().parent.parent
_SNAKEFILE = _WORKFLOW_DIR / "Snakefile"

LOGO = r"""
              /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /
             /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   
            /   /    ▗▄▄▖ ▗▄▖ ▗▄▄▖ ▗▄▄▄▖▗▖   ▗▖   ▗▄▄▄▖▗▖  ▗▖▗▄▄▄▖      /   /   /
           /   /    ▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌   ▐▌   ▐▌     █  ▐▛▚▖▐▌  █       /   /   /
          /   /     ▐▌   ▐▛▀▜▌▐▛▀▘ ▐▛▀▀▘▐▌   ▐▌     █  ▐▌ ▝▜▌  █      /   /   /
         /   /      ▝▚▄▄▖▐▌ ▐▌▐▌   ▐▙▄▄▖▐▙▄▄▖▐▙▄▄▖▗▄█▄▖▐▌  ▐▌▗▄█▄▖   /   /   /
        /   /   /                                                   /   /   /
       /   /   /       CRISPR-Abundance Phage-Evidence Linkage     /   /   /
      /   /   /                for Network Inference              /   /   /
     /   /   /   /   /      Snakemake Workflow Edition   /   /   /   /   /            
    /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   / 
"""




# ── Pipeline stages (Snakemake rule names) ───────────────────────────────────

STAGE_ORDER: list[str] = [
    "dada2",
    "harmonize",
    "spacepharer",
    "protein_clusters",
    "common_abundance",
    "shrinkage_correlations",
    "raw_crispr_network",
    "smooth_crispr",
    "xstar",
]

STAGE_LABELS: dict[str, str] = {
    "dada2":                   "DADA2 — 16S amplicon processing",
    "harmonize":               "Harmonize — NCBI / GCA assignment",
    "spacepharer":             "SpacePHARER — CRISPR spacer predictions",
    "protein_clusters":        "ProCs Maker — protein-cluster count + presence matrices",
    "common_abundance":        "Common abundance — align V & B matrices",
    "shrinkage_correlations":  "Shrinkage correlations — Schäfer-Strimmer",
    "raw_crispr_network":      "Raw CRISPR network — aggregate predictions",
    "smooth_crispr":           "Smooth CRISPR — taxonomy kernel smoothing",
    "xstar":                   "X* — residual propagation + correlations",
}


# ── Config section grouping for "Show config" ───────────────────────────────

CONFIG_SECTIONS: list[tuple[str, list[str]]] = [
    ("Execution", ["cores"]),
    ("Paths", [
        "base", "download_path",
        "silva_ref_path", "silva_taxmap_path",
        "bacterial_raw_fasta_folder", "virus_fasta_name",
        "metadata_path",
    ]),
    ("Resources (download URLs)", [
        "ncbi_taxdmp_url", "genes_reference_url",
        "bacContigs_reference_url", "protein_reference_url",
    ]),
    ("Global settings", ["species_level", "fresh_start", "ref_removal"]),
    ("DADA2", ["direction", "bacteria_fasta_name", "fasta_generation", "chimera_removal"]),
    ("Harmonizer (MMSeqs2)", [
        "isolate_ref_16S", "mapping_saving", "min_bitscore", "max_matches",
        "add_taxonomy", "extend_taxonomy",
    ]),
    ("SpacePHARER", [
        "min_n_spacers", "min_length", "max_length", "fdr",
        "keep_spacers_collection", "remove_decomp_fasta",
    ]),
    ("Protein clustering", [
        "matrix_type", "batch_size", "filter_1bac_1vir",
    ]),
    ("Network — run flags", [
        "run_common_abundance", "run_shrinkage_correlations",
        "run_raw_crispr_networks", "run_smooth_crispr", "run_xstar",
    ]),
    ("Network — common abundance", [
        "prevalence", "keep_column", "bacteria_taxonomy_rank",
    ]),
    ("Network — CRISPR smoothing", [
        "bacterial_ranks", "bacterial_weights",
        "viral_ranks", "viral_weights",
        "crispr_smooth_alpha", "transpose_raw_crispr_after_load",
        "aggregate_viral_rank",
    ]),
    ("Network — X*", ["pseudocount", "lam", "n_steps", "preserve_scale"]),
    ("Network — study inputs", [
        "virus_abundance_raw", "bacteria_otu", "bacteria_taxonomy",
        "phage_host_predictions", "tax_bac_for_smoothing", "tax_vir",
    ]),
]


# ── Config persistence ───────────────────────────────────────────────────────

def _read_last_config_path() -> Path | None:
    if LAST_CONFIG_FILE.exists():
        try:
            p = LAST_CONFIG_FILE.read_text().strip()
            if p:
                return Path(p)
        except OSError:
            pass
    return None


def _write_last_config_path(path: Path) -> None:
    LAST_CONFIG_FILE.parent.mkdir(parents=True, exist_ok=True)
    LAST_CONFIG_FILE.write_text(str(path))


# ── Config loading ───────────────────────────────────────────────────────────

def _load_yaml(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def _resolve_config_path(raw: str) -> Path | None:
    """Try CWD first, then fall back to the workflow config/ directory."""
    p = Path(raw)
    if p.exists():
        return p.resolve()
    wf = _WORKFLOW_DIR / raw
    if wf.exists():
        return wf.resolve()
    return None


# ── Editor ───────────────────────────────────────────────────────────────────

def _edit_file(path: Path) -> None:
    if shutil.which("micro") is None:
        CONSOLE.print(
            "[red]The 'micro' editor is required but was not found on PATH.[/red]\n"
            "Install it: https://github.com/zyedidia/micro  (e.g. `brew install micro`)."
        )
        return
    CONSOLE.print(
        "[dim]micro shortcuts:[/dim] "
        "[bold]Ctrl+S[/bold] save  ·  "
        "[bold]Ctrl+Q[/bold] quit  ·  "
        "[bold]Ctrl+G[/bold] full help"
    )
    _pause("Press Enter to open the editor…")
    subprocess.call(["micro", str(path)])


# ── Display helpers ──────────────────────────────────────────────────────────

def _show_logo() -> None:
    CONSOLE.print(f"[cyan]{LOGO}[/cyan]")


def _refresh_screen() -> None:
    CONSOLE.clear()
    _show_logo()


def _show_config(cfg: dict) -> None:
    for section, fields in CONFIG_SECTIONS:
        table = Table(title=section, show_header=True, header_style="bold cyan")
        table.add_column("Field", style="cyan", no_wrap=True)
        table.add_column("Value", style="white")
        for f in fields:
            if f in cfg:
                table.add_row(f, str(cfg[f]))
        if table.row_count:
            CONSOLE.print(table)


def _bundled_data_dir() -> Path:
    ref = pkg_resources.files("capellini_workflow").joinpath("data")
    return Path(str(ref))


def _bundled_reference_paths() -> tuple[Path, Path]:
    data = _bundled_data_dir()
    return (
        data / "references" / "progenome16S.fasta",
        data / "references" / "spacers" / "spacers_CompleteCollection.fasta",
    )


# ── Snakemake runner ─────────────────────────────────────────────────────────

def _build_snakemake_cmd(
    configfile: Path,
    cores: int = 4,
    targets: list[str] | None = None,
    dry_run: bool = False,
    extra: list[str] | None = None,
) -> list[str]:
    cmd = [
        "snakemake",
        "--snakefile", str(_SNAKEFILE),
        "--configfile", str(configfile),
        "--cores", str(cores),
        # If a previous run was interrupted (Ctrl-C, kill, power loss) while
        # a rule was mid-write, Snakemake flags the partial outputs in
        # .snakemake/incomplete_files.txt. Default behavior is to refuse the
        # next run with an IncompleteFilesException — we want it to just
        # silently re-run those rules so an interruption is never a foot-gun.
        "--rerun-incomplete",
    ]
    if dry_run:
        cmd.append("--dry-run")
    if targets:
        cmd.extend(targets)
    if extra:
        cmd.extend(extra)
    return cmd


def _run_snakemake_rule(
    configfile: Path,
    rule: str,
    cores: int = 4,
    dry_run: bool = False,
) -> subprocess.CompletedProcess:
    """Run a single Snakemake rule, streaming output to stdout/stderr."""
    cmd = _build_snakemake_cmd(configfile, cores, targets=[rule], dry_run=dry_run)
    return subprocess.run(cmd)


def _unlock_snakemake(configfile: Path) -> None:
    """Clear any stale Snakemake lock from an interrupted previous run.

    Snakemake locks the working directory while it runs; if the process is
    killed (Ctrl-C, OOM, power loss) the lock file is left behind and
    subsequent runs error out with a ``LockException``. This helper invokes
    ``snakemake --unlock`` as a best-effort cleanup. It is safe to call when
    no lock is present — Snakemake exits cleanly in that case.
    """
    cmd = [
        "snakemake",
        "--snakefile", str(_SNAKEFILE),
        "--configfile", str(configfile),
        "--unlock",
    ]
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except Exception:
        # Best-effort: never let unlock failures mask the real error.
        pass


def _run_dry_run_summary(configfile: Path, cores: int) -> None:
    """Run ``snakemake --dry-run`` and print only the Job Stats + Reasons sections.

    Snakemake's default dry-run output includes the full per-rule descriptor
    for every job (input/output/jobid/reason), plus a duplicated Job Stats
    block at the end. This helper captures the full output and emits only:

      1. The first ``Job stats:`` table (one occurrence).
      2. The ``Reasons:`` section, including the closing dry-run footer.
    """
    cmd = _build_snakemake_cmd(configfile, cores, dry_run=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    lines = (result.stdout + result.stderr).splitlines()

    # State machine: pre → stats → middle → reasons → done
    state = "pre"
    for line in lines:
        if state == "pre":
            if line.startswith("Job stats:"):
                state = "stats"
                CONSOLE.print(line)
        elif state == "stats":
            # First per-rule header marks the end of the stats table.
            if line.startswith("[") or line.startswith("rule "):
                state = "middle"
            else:
                CONSOLE.print(line)
        elif state == "middle":
            if line.startswith("Reasons:"):
                state = "reasons"
                CONSOLE.print()
                CONSOLE.print(line)
        elif state == "reasons":
            CONSOLE.print(line)

    if state == "pre":
        # Snakemake never produced a Job stats block — fall back to raw output
        # so the user can see whatever error message came through.
        CONSOLE.print("[yellow]Snakemake produced no Job stats block.[/yellow]")
        if result.stdout:
            CONSOLE.print(result.stdout)
        if result.stderr:
            CONSOLE.print(f"[red]{result.stderr}[/red]")

    if result.returncode != 0:
        CONSOLE.print(f"\n[red]Snakemake exited with code {result.returncode}[/red]")


# ── Validation ───────────────────────────────────────────────────────────────

# Labels of bundled references that can be downloaded via fetch_references().
# Only these block preflight and trigger the "Download now?" prompt.
_REFERENCE_LABELS = {
    "Bundled spacers_CompleteCollection.fasta",
}

# Labels of references the workflow regenerates automatically at runtime
# (e.g. progenome16S.fasta is filtered from the genes reference by the
# harmonizer the first time it's needed). These are shown in "Validate
# inputs" for visibility but never block the pipeline.
_AUTO_REGENERABLE_LABELS = {
    "Bundled progenome16S.fasta",
}


def _check_inputs(cfg: dict) -> tuple[list[tuple[str, bool, str]], list[tuple[str, bool, str]]]:
    """Run all input/dependency checks."""
    base = cfg.get("base", "")
    direction = cfg.get("direction", "forward")
    suffix = {"forward": "F", "reverse": "R", "paired": "P"}.get(direction, "F")
    fasta_folder = f"{base}/Inputs/Fasta Collection"

    paths: list[tuple[str, bool, str]] = []

    virus = Path(fasta_folder) / cfg.get("virus_fasta_name", "")
    paths.append(("Virus FASTA", virus.exists(), str(virus)))

    silva_ref = cfg.get("silva_ref_path", "")
    paths.append(("SILVA reference", bool(silva_ref) and Path(silva_ref).exists(), silva_ref))

    silva_tax = cfg.get("silva_taxmap_path", "")
    paths.append(("SILVA taxmap", bool(silva_tax) and Path(silva_tax).exists(), silva_tax))

    bac_raw = cfg.get("bacterial_raw_fasta_folder", "")
    paths.append(("Bacterial raw FASTA folder",
                  bool(bac_raw) and Path(bac_raw).is_dir(), bac_raw))

    meta = cfg.get("metadata_path", "")
    paths.append(("Metadata", bool(meta) and Path(meta).exists(), meta))

    vir_abund = cfg.get("virus_abundance_raw", "")
    paths.append(("Virus abundance raw",
                  bool(vir_abund) and Path(vir_abund).exists(), vir_abund))

    tax_vir = cfg.get("tax_vir", "")
    paths.append(("Viral taxonomy", bool(tax_vir) and Path(tax_vir).exists(), tax_vir))

    try:
        bundled_16s, bundled_spacers = _bundled_reference_paths()
        paths.append(("Bundled progenome16S.fasta", bundled_16s.exists(), str(bundled_16s)))
        paths.append(("Bundled spacers_CompleteCollection.fasta",
                      bundled_spacers.exists(), str(bundled_spacers)))
    except (ModuleNotFoundError, FileNotFoundError):
        paths.append(("Bundled progenome16S.fasta", False, "<not installed>"))
        paths.append(("Bundled spacers_CompleteCollection.fasta", False, "<not installed>"))

    deps: list[tuple[str, bool, str]] = []
    for tool in ("spacepharer", "minced", "Rscript", "prodigal", "micro"):
        deps.append((tool, shutil.which(tool) is not None, ""))
    mmseqs_ok = shutil.which("mmseqs") is not None or shutil.which("mmseqs2") is not None
    deps.append(("mmseqs/mmseqs2", mmseqs_ok, ""))
    deps.append(("snakemake", shutil.which("snakemake") is not None, ""))

    return paths, deps


def _validate_inputs(cfg: dict) -> None:
    paths, deps = _check_inputs(cfg)

    def make_table(title: str) -> Table:
        t = Table(title=title)
        t.add_column("Check", style="cyan", no_wrap=True)
        t.add_column("Result")
        return t

    def add_row(t: Table, name: str, ok: bool, detail: str = "") -> None:
        marker = "[green]OK[/green]" if ok else "[red]MISSING[/red]"
        t.add_row(name, f"{marker} {detail}".strip())

    paths_table = make_table("Input paths")
    for n, ok, d in paths:
        add_row(paths_table, n, ok, d)

    deps_table = make_table("External dependencies")
    for n, ok, d in deps:
        add_row(deps_table, n, ok, d)

    CONSOLE.print(paths_table)
    CONSOLE.print()
    CONSOLE.print(deps_table)


# ── Preflight ────────────────────────────────────────────────────────────────

def _preflight_full_pipeline(session: "_Session") -> bool:
    """Validate inputs before launching. Returns True if the pipeline can proceed."""
    if not session.ensure_loaded():
        _pause()
        return False
    cfg = session.config
    paths, deps = _check_inputs(cfg)

    missing_refs = [(n, d) for (n, ok, d) in paths if not ok and n in _REFERENCE_LABELS]
    # Real input paths that block the pipeline — exclude downloadable references
    # AND auto-regenerable references (regenerated by the workflow at runtime).
    missing_paths = [
        (n, d) for (n, ok, d) in paths
        if not ok
        and n not in _REFERENCE_LABELS
        and n not in _AUTO_REGENERABLE_LABELS
    ]
    missing_deps = [(n, d) for (n, ok, d) in deps if not ok]

    # Only references missing → offer to download
    if missing_refs and not missing_paths and not missing_deps:
        _refresh_screen()
        CONSOLE.print("[yellow]References not found:[/yellow]")
        for n, _ in missing_refs:
            CONSOLE.print(f"  • {n}")
        CONSOLE.print()
        if _confirm("Download the missing references now?", default=True):
            _refresh_screen()
            from capellini_workflow.fetch_references import fetch_references
            try:
                fetch_references(overwrite=False)
            except RuntimeError as exc:
                CONSOLE.print(f"[red]{exc}[/red]")
                _pause()
                return False
            _pause()
            return _preflight_full_pipeline(session)
        return False

    # Anything else missing → abort with summary
    if missing_refs or missing_paths or missing_deps:
        _refresh_screen()
        CONSOLE.print("[red]Cannot start the full pipeline — missing items:[/red]\n")
        if missing_refs:
            CONSOLE.print("[bold]References:[/bold]")
            for n, d in missing_refs:
                CONSOLE.print(f"  • {n}  [dim]{d}[/dim]")
            CONSOLE.print()
        if missing_paths:
            CONSOLE.print("[bold]Input paths:[/bold]")
            for n, d in missing_paths:
                CONSOLE.print(f"  • {n}  [dim]{d or '<empty>'}[/dim]")
            CONSOLE.print()
        if missing_deps:
            CONSOLE.print("[bold]External dependencies:[/bold]")
            for n, _ in missing_deps:
                CONSOLE.print(f"  • {n}  [dim]not on PATH[/dim]")
            CONSOLE.print()
        _select(
            "",
            choices=[questionary.Choice(" » Back to main menu", "back")],
            default="back",
        )
        return False

    return True


# ── Live progress panel ──────────────────────────────────────────────────────

_TAIL_MAX_LINES = 60
_LOG_DIR = Path.home() / ".capellini" / "logs"


def _print_failure(exc: BaseException) -> None:
    """Print the failure message and, if attached, the full-log file path.

    ``_live_progress`` attaches ``failure_log_path`` to the raised exception
    so the caller can surface the saved log without re-running the stage.
    """
    CONSOLE.print(f"\n[red]{exc}[/red]")
    log_path = getattr(exc, "failure_log_path", None)
    if log_path is not None:
        CONSOLE.print(
            f"[dim]Full stage log saved to:[/dim] [bold]{log_path}[/bold]\n"
            f"[dim]Open it with[/dim] [bold]less '{log_path}'[/bold] "
            f"[dim]or[/dim] [bold]micro '{log_path}'[/bold] "
            f"[dim]to see the real error[/dim]"
        )


def _save_failure_log(stage: str, lines: list[str]) -> Path:
    """Dump the full captured output of a failed stage to ``~/.capellini/logs``.

    Returns the path to the written log file so the caller can surface it to
    the user. The bounded tail buffer only keeps the last ``_TAIL_MAX_LINES``
    lines, but the real root cause is usually buried earlier (e.g. a Python
    traceback that scrolled off). The full log preserves everything.
    """
    _LOG_DIR.mkdir(parents=True, exist_ok=True)
    import time
    ts = time.strftime("%Y%m%d-%H%M%S")
    log_path = _LOG_DIR / f"{stage}_{ts}.log"
    log_path.write_text("\n".join(lines) + "\n")
    return log_path


def _live_progress(stages: list[str], session: "_Session", dry_run: bool = False) -> None:
    """Run the given Snakemake rules sequentially with a live Rich progress panel."""
    statuses = {s: "waiting" for s in stages}
    tail: deque[str] = deque(maxlen=_TAIL_MAX_LINES)
    # Unbounded copy of the current stage's output, used to dump a full log on
    # failure. Cleared at the start of every stage so we only persist the
    # output of the stage that actually failed.
    full_log: list[str] = []
    current_stage: dict[str, str] = {"name": ""}
    lock = threading.Lock()

    def render_table() -> Table:
        t = Table(title="CAPELLINI workflow")
        t.add_column("Stage", style="cyan")
        t.add_column("Status")
        for s in stages:
            mark = {
                "waiting": "[grey50]waiting[/grey50]",
                "running": "[yellow]running[/yellow]",
                "done":    "[green]✓ done[/green]",
                "skipped": "[blue]● up to date[/blue]",
                "failed":  "[red]✗ failed[/red]",
            }[statuses[s]]
            t.add_row(STAGE_LABELS.get(s, s), mark)
        return t

    def render() -> Group:
        with lock:
            body_text = "\n".join(tail) if tail else "[dim](no output yet)[/dim]"
            title = STAGE_LABELS.get(current_stage["name"], current_stage["name"]) or "output"
        panel = Panel(
            Text.from_markup(body_text),
            title=f"[bold]{title}[/bold] (last {_TAIL_MAX_LINES} lines)",
            border_style="grey50",
            height=_TAIL_MAX_LINES + 2,
        )
        return Group(render_table(), Text(""), panel)

    class _Refreshable:
        def __rich__(self):
            return render()

    # Save real terminal fds before redirecting
    saved_stdout_fd = os.dup(1)
    saved_stderr_fd = os.dup(2)
    real_term = os.fdopen(saved_stdout_fd, "w", buffering=1)
    live_console = Console(file=real_term, force_terminal=True)

    r_fd, w_fd = os.pipe()
    saved_py_stdout = sys.stdout
    saved_py_stderr = sys.stderr

    def _is_progress_bar(line: str) -> bool:
        """Heuristic: MMseqs2-style progress lines like ``[==== ] 3.00K 0s 100ms``.

        Tools that print these to a pipe emit one full line per percent step
        (instead of using ``\\r`` on a TTY), which floods the live panel and
        causes flicker. We collapse consecutive progress lines into one.
        """
        s = line.lstrip()
        return s.startswith("[") and "]" in s and ("=" in s or " " in s.split("]", 1)[0])

    def reader_loop() -> None:
        buf = bytearray()
        while True:
            try:
                chunk = os.read(r_fd, 4096)
            except OSError:
                return
            if not chunk:
                return
            buf.extend(chunk)
            while True:
                nl = buf.find(b"\n")
                if nl < 0:
                    break
                line = bytes(buf[:nl]).decode("utf-8", errors="replace").rstrip("\r")
                del buf[: nl + 1]
                with lock:
                    # Collapse consecutive progress-bar lines: if the previous
                    # tail entry was also a progress bar, overwrite it instead
                    # of appending. The full_log keeps every line for the
                    # failure dump so nothing is lost on disk.
                    if (
                        tail
                        and _is_progress_bar(line)
                        and _is_progress_bar(tail[-1])
                    ):
                        tail[-1] = line
                    else:
                        tail.append(line)
                    full_log.append(line)

    reader = threading.Thread(target=reader_loop, daemon=True)

    configfile = session.config_path
    cores = session.cores

    # Proactive unlock: clear any stale lock left by a previous interrupted run
    # so the user never has to manually run `snakemake --unlock`.
    _unlock_snakemake(configfile)

    try:
        os.dup2(w_fd, 1)
        os.dup2(w_fd, 2)
        os.close(w_fd)
        sys.stdout = os.fdopen(1, "w", buffering=1, closefd=False)
        sys.stderr = os.fdopen(2, "w", buffering=1, closefd=False)
        reader.start()

        failure_log_path: Path | None = None
        # 4 fps is enough to feel "live" without flickering. The reader_loop
        # collapses consecutive MMseqs2-style progress bars (see _is_progress_bar)
        # so the bounded tail doesn't churn at hundreds of lines per second.
        with Live(_Refreshable(), console=live_console, refresh_per_second=4):
            for s in stages:
                with lock:
                    tail.clear()
                    full_log.clear()
                    current_stage["name"] = s
                statuses[s] = "running"
                try:
                    cmd = _build_snakemake_cmd(
                        configfile, cores, targets=[s], dry_run=dry_run,
                    )
                    result = subprocess.run(cmd)
                    if result.returncode != 0:
                        statuses[s] = "failed"
                        raise RuntimeError(
                            f"Snakemake rule '{s}' failed (exit {result.returncode})"
                        )
                except KeyboardInterrupt:
                    # User pressed Ctrl-C — Snakemake child died holding the lock.
                    statuses[s] = "failed"
                    with lock:
                        failure_log_path = _save_failure_log(s, list(full_log))
                    _unlock_snakemake(configfile)
                    raise
                except RuntimeError as exc:
                    # Stage exited non-zero — release the lock and persist
                    # the full output before bubbling up. Attach the log path
                    # to the exception so the menu wrapper can show it.
                    with lock:
                        failure_log_path = _save_failure_log(s, list(full_log))
                    exc.failure_log_path = failure_log_path  # type: ignore[attr-defined]
                    _unlock_snakemake(configfile)
                    raise
                except Exception as exc:
                    statuses[s] = "failed"
                    with lock:
                        failure_log_path = _save_failure_log(s, list(full_log))
                    exc.failure_log_path = failure_log_path  # type: ignore[attr-defined]
                    _unlock_snakemake(configfile)
                    raise
                statuses[s] = "done"
    finally:
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        except Exception:
            pass
        os.dup2(saved_stdout_fd, 1)
        os.dup2(saved_stderr_fd, 2)
        sys.stdout = saved_py_stdout
        sys.stderr = saved_py_stderr
        os.close(saved_stderr_fd)
        reader.join(timeout=2)
        try:
            os.close(r_fd)
        except OSError:
            pass
        try:
            real_term.close()
        except Exception:
            pass


# ── Config session ───────────────────────────────────────────────────────────

class _Session:
    def __init__(self) -> None:
        self.config_path: Path | None = None
        self.config: dict | None = None
        self.cores: int = 4

    def ensure_loaded(self) -> bool:
        if self.config is not None:
            return True
        path = self._resolve_initial_path()
        if path is None:
            return False
        return self._load(path)

    def _resolve_initial_path(self) -> Path | None:
        path = _read_last_config_path()
        if path is not None and path.exists():
            return path
        CONSOLE.print(
            "[yellow]No previous config remembered.[/yellow] "
            "Use Settings → Load config to point CAPELLINI at one."
        )
        return None

    def _load(self, path: Path) -> bool:
        if not path.exists():
            CONSOLE.print(f"[red]Config not found:[/red] {path}")
            return False
        try:
            self.config = _load_yaml(path)
        except Exception as exc:
            CONSOLE.print(f"[red]Failed to load config:[/red] {exc}")
            return False
        self.config_path = path.resolve()
        # Pull cores from YAML so it persists across sessions. The Settings
        # menu can still override it for the current session.
        yaml_cores = self.config.get("cores")
        if isinstance(yaml_cores, int) and yaml_cores > 0:
            self.cores = yaml_cores
        CONSOLE.print(f"[green]Loaded config:[/green] {path}")
        return True

    def reload(self) -> None:
        if self.config_path is not None:
            self._load(self.config_path)

    def switch(self, path: Path) -> bool:
        return self._load(path)


# ── questionary wrappers ─────────────────────────────────────────────────────

def _select(message: str, choices: list, default=None):
    try:
        return questionary.select(message, choices=choices, default=default, qmark="»").ask()
    except KeyboardInterrupt:
        return None


def _ask_path(message: str, default: str) -> Path | None:
    answer = questionary.path(message, default=default).ask()
    if answer is None or answer.strip() == "":
        return None
    return Path(answer).expanduser()


def _confirm(message: str, default: bool = True) -> bool:
    return bool(questionary.confirm(message, default=default).ask())


def _pause(message: str = "Press Enter to return…") -> None:
    try:
        input(f"\n{message}")
    except (KeyboardInterrupt, EOFError):
        pass


# ── Settings sub-menu ────────────────────────────────────────────────────────

def _settings_menu(session: _Session) -> None:
    while True:
        _refresh_screen()
        current = session.config_path
        CONSOLE.print(
            f"[dim]Current config:[/dim] {current if current else '[red]<none>[/red]'}\n"
            f"[dim]Cores:[/dim] {session.cores}\n"
        )
        choice = _select(
            "Settings",
            choices=[
                questionary.Choice("Load config (set / change current)", "load"),
                questionary.Choice("Edit current config", "edit"),
                questionary.Choice("Show current config", "show"),
                questionary.Choice("Set number of cores", "cores"),
                questionary.Choice("Validate inputs", "validate"),
                questionary.Separator(),
                questionary.Choice(" » Back to main menu", "back"),
            ],
            default="back",
        )
        if choice in (None, "back"):
            return

        if choice == "load":
            default_str = str(current) if current else str(Path.cwd())
            target = _ask_path("Config path", default=default_str)
            if target is None:
                continue
            if not target.exists():
                CONSOLE.print(f"[red]Not found:[/red] {target}")
                _pause()
                continue
            if session.switch(target):
                _write_last_config_path(target)
                CONSOLE.print("[green]Config remembered for next time.[/green]")

        elif choice == "edit":
            if current is None:
                CONSOLE.print(
                    "[yellow]No config loaded.[/yellow] Use 'Load config' first."
                )
                _pause()
                continue
            _edit_file(current)
            session.reload()

        elif choice == "show":
            if not session.ensure_loaded():
                _pause()
                continue
            _refresh_screen()
            _show_config(session.config)

        elif choice == "cores":
            try:
                available = os.cpu_count() or 1
                CONSOLE.print(
                    f"[dim]Detected logical CPU cores on this machine: "
                    f"[bold]{available}[/bold]   (current: {session.cores})[/dim]"
                )
                answer = questionary.text(
                    f"Number of cores (1–{available})",
                    default=str(session.cores),
                ).ask()
                if answer and answer.strip().isdigit():
                    n = int(answer.strip())
                    if n < 1:
                        CONSOLE.print("[red]Must be at least 1.[/red]")
                    else:
                        if n > available:
                            CONSOLE.print(
                                f"[yellow]Warning: {n} exceeds the detected {available} cores. "
                                f"Snakemake will still run but may oversubscribe the CPU.[/yellow]"
                            )
                        session.cores = n
                        CONSOLE.print(f"[green]Cores set to {session.cores}[/green]")
            except (KeyboardInterrupt, ValueError):
                pass

        elif choice == "validate":
            if not session.ensure_loaded():
                _pause()
                continue
            _refresh_screen()
            _validate_inputs(session.config)

        _pause()


# ── Stage sub-menus ──────────────────────────────────────────────────────────

def _selected_stages_menu(session: _Session) -> None:
    if not session.ensure_loaded():
        return
    while True:
        _refresh_screen()
        CONSOLE.print(
            "[dim]Toggle a stage with [space], confirm with [enter].\n"
            "To go back, do not select any stage and press enter.[/dim]\n"
        )
        try:
            picked = questionary.checkbox(
                "Select stages to run",
                choices=[questionary.Choice(STAGE_LABELS.get(s, s), s) for s in STAGE_ORDER],
                qmark="»",
            ).ask()
        except KeyboardInterrupt:
            return
        if picked is None or not picked:
            return
        ordered = [s for s in STAGE_ORDER if s in picked]
        _refresh_screen()
        try:
            _live_progress(ordered, session)
        except RuntimeError as exc:
            _print_failure(exc)
        _pause()


def _single_stage_menu(session: _Session) -> None:
    if not session.ensure_loaded():
        return
    while True:
        _refresh_screen()
        options = [questionary.Choice(STAGE_LABELS.get(s, s), s) for s in STAGE_ORDER]
        options.extend([questionary.Separator(), questionary.Choice(" » Back to main menu", "back")])

        choice = _select("Run a single stage", choices=options, default="back")
        if choice in (None, "back"):
            return
        _refresh_screen()
        try:
            _live_progress([choice], session)
        except RuntimeError as exc:
            _print_failure(exc)
        _pause()


def _run_pipeline_menu(session: _Session) -> None:
    while True:
        _refresh_screen()
        choice = _select(
            "Run pipeline",
            choices=[
                questionary.Choice("Run full pipeline", "run_all"),
                questionary.Choice("Run selected stages", "run_selected"),
                questionary.Choice("Run single stage", "run_one"),
                questionary.Choice("Dry run (show what would be done)", "dry_run"),
                questionary.Separator(),
                questionary.Choice(" » Back to main menu", "back"),
            ],
            default="run_all",
        )
        if choice in (None, "back"):
            return
        if choice == "run_selected":
            _selected_stages_menu(session)
            continue
        if choice == "run_one":
            _single_stage_menu(session)
            continue
        if choice == "dry_run":
            if not session.ensure_loaded():
                _pause()
                continue
            _refresh_screen()
            _run_dry_run_summary(session.config_path, session.cores)
            _pause()
            continue
        if choice == "run_all":
            if not _preflight_full_pipeline(session):
                continue
            _refresh_screen()
            try:
                _live_progress(list(STAGE_ORDER), session)
            except RuntimeError as exc:
                _print_failure(exc)
            _pause()


# ── Direct CLI mode (non-interactive) ────────────────────────────────────────

def _cli_run(argv: list[str]) -> int:
    """``capellini run --configfile …`` — direct Snakemake dispatch."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="capellini run",
        description="Run CAPELLINI workflow via Snakemake (non-interactive).",
    )
    parser.add_argument(
        "--configfile", required=True,
        help="Path to YAML config file (e.g. config/ibd.yaml)",
    )
    parser.add_argument(
        "--cores", "-j", type=int, default=None,
        help="Number of cores for Snakemake "
             "(default: value of `cores` in YAML, falling back to 4)",
    )
    parser.add_argument(
        "--target", nargs="*", default=None,
        help="Specific Snakemake target(s) to build (default: run all)",
    )
    parser.add_argument(
        "--dry-run", "-n", action="store_true",
        help="Show what would be done without executing",
    )
    parser.add_argument(
        "--use-conda", action="store_true",
        help="Use conda environments for external tools",
    )
    parser.add_argument(
        "--snakemake-args", nargs=argparse.REMAINDER, default=[],
        help="Additional arguments passed directly to Snakemake",
    )

    args = parser.parse_args(argv)

    configfile = _resolve_config_path(args.configfile)
    if configfile is None:
        print(f"ERROR: Config file not found: {args.configfile}", file=sys.stderr)
        return 1

    if not _SNAKEFILE.exists():
        print(f"ERROR: Snakefile not found: {_SNAKEFILE}", file=sys.stderr)
        return 1

    # Resolve cores: CLI flag wins, else YAML `cores`, else 4.
    if args.cores is None:
        try:
            yaml_cores = _load_yaml(configfile).get("cores")
            args.cores = yaml_cores if isinstance(yaml_cores, int) and yaml_cores > 0 else 4
        except Exception:
            args.cores = 4

    extra = []
    if args.use_conda:
        extra.append("--use-conda")
    if args.snakemake_args:
        extra.extend(args.snakemake_args)

    cmd = _build_snakemake_cmd(
        configfile, args.cores,
        targets=args.target, dry_run=args.dry_run, extra=extra,
    )
    print(f"Running: {' '.join(cmd)}\n")

    # Proactive unlock so a Ctrl-C from a previous run doesn't block this one.
    # Skip for dry-run since dry-run never acquires the lock.
    if not args.dry_run:
        _unlock_snakemake(configfile)

    try:
        result = subprocess.run(cmd)
    except KeyboardInterrupt:
        # Snakemake child got SIGINT — release the lock before exiting.
        if not args.dry_run:
            _unlock_snakemake(configfile)
        return 130  # standard SIGINT exit code

    if result.returncode != 0 and not args.dry_run:
        _unlock_snakemake(configfile)
    return result.returncode


# ── Main entry point ─────────────────────────────────────────────────────────

def main() -> None:
    """Entry point for the ``capellini`` CLI.

    - ``capellini``            → interactive menu
    - ``capellini run …``      → direct Snakemake dispatch
    - ``capellini fetch-references`` → download bundled references
    """
    # Sub-command dispatch
    if len(sys.argv) > 1:
        sub = sys.argv[1]
        if sub in {"run", "--configfile"}:
            # "capellini run --configfile ..." or legacy "capellini --configfile ..."
            argv = sys.argv[2:] if sub == "run" else sys.argv[1:]
            sys.exit(_cli_run(argv))
        if sub in {"fetch-references", "fetch_references"}:
            from capellini_workflow.fetch_references import fetch_references
            try:
                fetch_references(overwrite="--overwrite" in sys.argv)
            except RuntimeError as exc:
                print(f"ERROR: {exc}", file=sys.stderr)
                sys.exit(1)
            sys.exit(0)

    # Interactive mode
    _show_logo()
    session = _Session()
    session.ensure_loaded()

    while True:
        _refresh_screen()

        choice = _select(
            "Main menu",
            choices=[
                questionary.Choice("Run pipeline", "run"),
                questionary.Choice("Settings", "settings"),
                questionary.Choice("Fetch/Update reference FASTAs from GitHub release", "fetch_refs"),
                questionary.Separator(),
                questionary.Choice("Quit", "quit"),
            ],
            default="run",
        )
        if choice in (None, "quit"):
            return
        if choice == "run":
            _run_pipeline_menu(session)
            continue
        if choice == "settings":
            _settings_menu(session)
            continue
        if choice == "fetch_refs":
            _refresh_screen()
            from capellini_workflow.fetch_references import fetch_references
            overwrite = _confirm(
                "Re-download even if the files already exist?", default=False
            )
            try:
                fetch_references(overwrite=overwrite)
            except RuntimeError as exc:
                CONSOLE.print(f"[red]{exc}[/red]")
            _pause()
            continue


if __name__ == "__main__":
    sys.exit(main() or 0)
