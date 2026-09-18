"""Figure-contract tests for the BIT figure set (Task 11), written before the generator.

They pin the frozen plotting-data inputs, the output names, the final size and typography, legend/text
clearance from data, embedded-font export settings, and deterministic regeneration.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
FIG_DIR = ROOT / "figures_bit"
MANIFEST = ROOT / "results" / "bit_revision" / "figures" / "manifest.json"

# SHA-256 of every frozen v3.4 data file the BIT generator is allowed to read (sealed archive bytes).
SOURCE_HASHES = {
    "results/jcp_figures_v34/figure1_exact_branches.csv": "2f81dd15a03ac5e05f5fff5908386448eb995584f6b8133cd3649af8bcf50be4",
    "results/jcp_figures_v34/figure1_conditioning.csv": "8c6673442bbefa5abcd00c6012f5bcaf56111fa4f74cbf0da01564f6d8e02c62",
    "results/jcp_figures_v34/figure2_candidate_regions.csv": "00e0c7bb28e3429513513a91a97ea71f606e0888868b4d5e49cd3aac7f27a292",
    "results/jcp_figures_v34/figure2_richardson_margin.csv": "1e21c8612fd7246935ffaeefa8a8af0ebe7ccf388c4b815effea35490a44702b",
    "results/jcp_figures_v34/figure3_noncommuting_endpoints.csv": "9e7d52aec53f681ee56bec687642f9e1c821761eb118670f91babcfe84050623",
    "results/jcp_figures_v34/figure4_conditioning.csv": "895c2e3af0e0536858b8b137ea781cfbbb6c1e55edefddcb21a4db2af0fd96d9",
    "results/jcp_figures_v34/figure5_tolerance_sweep.csv": "aa7ccd3a71de7f0acdd270680573c5f7baa0a0f91e4e7c1780054ff81f840738",
    "results/jcp_figures_v34/figure5_action_composition.csv": "bb40d3337fd37354032c16620d74069e9df9e84aaaae7059c6b72798235cc52e",
    "results/jcp_figures_v34/figureS2_operational.csv": "b2874b7e48a845ffe03dfe6f5d441bc37302c71929507284071bd107cb34705f",
    "results/v27/v27_verification.json": "ce81fa596dc3f6d0b6fd5089c13b7c406e92585ffa32d32cc08f34c2e9418030",
    "results/v34/v34_verification.json": "2a908280917f9cf1ff607fce19364cfcfb2691eb137a08fb530120a8dc7a14eb",
}

MAIN = ["Fig1", "Fig2", "Fig3"]
SUPP = ["FigS1", "FigS2", "FigS3", "FigS4"]
WIDTH_MM = 119.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def generator():
    import importlib

    return importlib.import_module("generate_bit_figures")


@pytest.fixture(scope="module")
def manifest(generator):
    assert MANIFEST.exists(), "run scripts/generate_bit_figures.py first"
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def test_generator_declares_only_frozen_sources(generator):
    declared = set(generator.SOURCE_FILES)
    assert declared == set(SOURCE_HASHES)


@pytest.mark.parametrize("rel", sorted(SOURCE_HASHES))
def test_source_data_unchanged(rel, manifest):
    """The generator records the hash of each input it read; the file must still hash the same."""
    path = ROOT / rel
    assert path.exists()
    assert sha256(path) == SOURCE_HASHES[rel], "frozen v3.4 plotting data changed"
    assert manifest["sources"][rel] == SOURCE_HASHES[rel]


def test_output_files_exist(manifest):
    for stem in MAIN + SUPP:
        for ext in ("eps", "pdf", "png"):
            assert (FIG_DIR / f"{stem}.{ext}").exists(), f"{stem}.{ext} missing"
    assert not any(FIG_DIR.glob("Graphical_Abstract*"))


def test_manifest_hashes_match_disk(manifest):
    for name, entry in manifest["outputs"].items():
        assert sha256(FIG_DIR / name) == entry["sha256"], name


@pytest.mark.parametrize("stem", MAIN + SUPP)
def test_final_size_and_typography(generator, stem):
    fig = generator.build(stem)
    w_in, h_in = fig.get_size_inches()
    assert abs(w_in * 25.4 - WIDTH_MM) < 0.05
    assert h_in * 25.4 <= 195.0
    sizes = [t.get_fontsize() for t in generator.all_text_artists(fig) if t.get_text().strip()]
    assert sizes, "no text"
    assert min(sizes) >= 8.0 - 1e-9 and max(sizes) <= 10.0 + 1e-9, (min(sizes), max(sizes))
    for ax in fig.get_axes():
        assert ax.get_title() == "", "titles must live in captions"
    generator.close(fig)


@pytest.mark.parametrize("stem", MAIN + SUPP)
def test_legends_and_texts_clear_the_data(generator, stem):
    fig = generator.build(stem)
    report = generator.overlap_audit(fig)
    generator.close(fig)
    assert report["legend_data_overlaps"] == [], report["legend_data_overlaps"]
    assert report["text_data_overlaps"] == [], report["text_data_overlaps"]
    assert report["text_text_overlaps"] == [], report["text_text_overlaps"]


def test_series_counts(generator):
    counts = generator.series_counts()
    assert counts["Fig1"]["a_fills"] == 2 and counts["Fig1"]["inset"] == 1
    assert counts["Fig2"]["a_intervals"] == 2 and counts["Fig2"]["b_lines"] == 1
    assert counts["Fig3"]["a_rows"] == 5 and counts["Fig3"]["b_marker_series"] == 3 and counts["Fig3"]["b_connecting_lines"] == 0
    assert counts["FigS4"]["a_series"] == 3 and counts["FigS4"]["b_categories"] == 4


def test_export_settings_embed_fonts(generator):
    rc = generator.EXPORT_RC
    assert rc["pdf.fonttype"] == 42
    assert rc["ps.fonttype"] in (3, 42)
    assert "Arial" in rc["font.sans-serif"][0]


def test_deterministic_regeneration(generator, tmp_path):
    os.environ["SOURCE_DATE_EPOCH"] = "1787184000"
    a = generator.render_to(tmp_path / "a", "Fig2")
    b = generator.render_to(tmp_path / "b", "Fig2")
    for ext in ("pdf", "eps", "png"):
        assert (a / f"Fig2.{ext}").read_bytes() == (b / f"Fig2.{ext}").read_bytes(), ext
