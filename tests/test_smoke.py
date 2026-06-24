"""Smoke tests for capellini-workflow."""

import subprocess
import sys


def test_import():
    """Package imports cleanly."""
    import capellini_workflow
    assert hasattr(capellini_workflow, "__version__")


def test_submodules():
    """All library submodules import."""
    from capellini_workflow import transforms, taxonomy, io, network_utils
    assert callable(transforms.clr)
    assert callable(transforms.schaefer_strimmer_corr)
    assert callable(taxonomy.load_bacteria_taxonomy)
    assert callable(io.read_table)
    assert callable(network_utils.build_xstar_from_smoothed_crispr)


def test_fetch_references_import():
    """fetch_references module imports."""
    from capellini_workflow import fetch_references
    assert callable(fetch_references.fetch_references)


def test_scripts_parse():
    """All scripts parse without syntax errors."""
    scripts = [
        "scripts/spacepharer_wrapper.py",
        "scripts/common_abundance.py",
        "scripts/shrinkage_correlations.py",
        "scripts/raw_crispr_network.py",
        "scripts/smooth_crispr.py",
        "scripts/xstar.py",
        "scripts/extract_gca_list.py",
    ]
    for script in scripts:
        result = subprocess.run(
            [sys.executable, "-c", f"import py_compile; py_compile.compile('{script}', doraise=True)"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"{script} has syntax errors: {result.stderr}"
