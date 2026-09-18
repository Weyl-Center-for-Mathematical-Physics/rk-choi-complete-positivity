from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "manuscript_contract"


def _read(relative: str) -> str:
    return (SOURCE / relative).read_text(encoding="utf-8")


def test_submission_manifest_uses_consistent_publication_names() -> None:
    entries = set(_read("SOURCE_MANIFEST.txt").splitlines())
    required = {
        "Abstract.txt",
        "Manuscript_JCP.tex",
        "Manuscript_JCP_clean.tex",
        "Supplementary_Material_JCP.tex",
        "figures/Fig2_robustness.pdf",
        "figures/Fig3_candidate_maps.pdf",
        "figures/FigS1_frequency_bands.pdf",
        "sections/05_robustness.tex",
        "sections/06_candidate_maps.tex",
        "sections/07_step_control.tex",
        "sections/08_operational_composition.tex",
        "sections/09_comparison.tex",
        "sections/10_conclusion.tex",
    }
    forbidden = {
        "Manuscript_JCP_review.tex",
        "Supplemental_Material_JCP.tex",
        "figures/Fig2_candidate_maps.pdf",
        "figures/Fig3_robustness.pdf",
        "figures/Fig6_frequency_bands.pdf",
        "sections/06_step_control.tex",
        "sections/06_operational_composition.tex",
        "sections/07_comparison.tex",
        "sections/08_conclusion.tex",
    }
    assert required <= entries
    assert not (forbidden & entries)


def test_article_order_abstract_and_submission_wrapper_are_consistent() -> None:
    article = _read("jcp_article.tex")
    ordered_inputs = [
        r"\input{sections/05_robustness}",
        r"\input{sections/06_candidate_maps}",
        r"\input{sections/07_step_control}",
        r"\input{sections/08_operational_composition}",
        r"\input{sections/09_comparison}",
        r"\input{sections/10_conclusion}",
    ]
    positions = [article.index(item) for item in ordered_inputs]
    assert positions == sorted(positions)

    abstract = article.split(r"\begin{abstract}", 1)[1].split(
        r"\end{abstract}", 1
    )[0]
    assert abstract.strip().replace("--", "–") == _read("Abstract.txt").strip()
    assert "coherent precession, bidirectional relaxation, pure dephasing" in abstract
    assert "precession-to-relaxation ratio" in abstract
    assert "executed binary64 step and candidate" in abstract
    for undefined_notation in (r"\theta", r"\nu", r"\varpi", r"\alpha_4", r"\Omega_x"):
        assert undefined_notation not in abstract

    wrapper = _read("Manuscript_JCP.tex")
    assert r"\documentclass[review,12pt]{elsarticle}" in wrapper
    assert r"\usepackage[switch]{lineno}" in wrapper
    assert r"\newcommand{\ReviewMode}{}" in wrapper
    assert r"\linenumbers" in article
    assert r"\usepackage[T1]{fontenc}" in _read("jcp_preamble.tex")


def test_highlights_are_separate_specific_and_within_elsevier_limit() -> None:
    highlights = [
        line.removeprefix("- ")
        for line in _read("Highlights.txt").splitlines()
        if line.strip()
    ]
    assert highlights == [
        "Exact CPTP conditions are derived for Runge--Kutta steps on a qubit generator.",
        "RK4 admissible step sets can be disconnected or reduce to isolated points.",
        "Full-step and two-half-step RK4 maps can have disjoint CPTP sets.",
        "RK4 Richardson extrapolation can fail although both inputs are channels.",
        "A tri-state binary64 guard certifies the executed candidate before acceptance.",
    ]
    assert all(len(line) <= 85 for line in highlights)


def test_cover_letter_contains_editorial_declarations_fit_and_reviewers() -> None:
    cover = _read("Cover_Letter_JCP.tex")
    for required in (
        "August 20, 2026",
        "original",
        "not under consideration by any other journal",
        "approved by all authors",
        "efficacy",
        "robustness",
        "computational complexity",
        "reproducibility",
        "Riesch and Jirauschek",
        "Havasi and Kazemi",
        "Appel",
        "Cao and Lu",
        "Christian Jirauschek",
        "jirauschek@tum.de",
        "Yingda Cheng",
        "yingda@vt.edu",
        "David I. Ketcheson",
        "david.ketcheson@kaust.edu.sa",
        "Reviewer exclusions: none",
    ):
        assert required in cover
    for removed in (
        "15 locator calls",
        "PI history",
        "work proxy",
        "graphical abstract",
    ):
        assert removed.lower() not in cover.lower()


def test_supplementary_material_uses_consistent_notation_and_metadata() -> None:
    supplement = _read("Supplementary_Material_JCP.tex")
    assert r"\usepackage[letterpaper,margin=18mm]{geometry}" in supplement
    assert "Supplementary Material" in supplement
    assert "Supplemental Material" not in supplement
    assert r"G. Blake Pierpoint$^{1,*}$, Olivier Bernard$^{2,3}$, and Yichen Liu$^{4,3}$" in supplement
    assert r"$\kappa=\gamma_\phi/\Gamma$" in supplement
    assert r"$\varpi=|\omega|/\Gamma$" in supplement
    assert "FigS1_frequency_bands.pdf" in supplement
    assert "Fig6_frequency_bands" not in supplement
    assert "share only $x=0$" in supplement
    assert "Thus $x=1$" in supplement
    assert "$H=0$" not in supplement
    assert "$H=1$" not in supplement
    assert "$H=2$" not in supplement
    assert r"$\Phi$ is completely positive if and only if $\mathcal U_h\circ\Phi$ is completely positive" in supplement


def test_credit_acknowledgment_data_and_related_work_are_submission_ready() -> None:
    article = _read("jcp_article.tex")
    expected_credit = (
        "G. Blake Pierpoint: Conceptualization, Methodology, Formal analysis, "
        "Investigation, Data curation, Visualization, Project administration, "
        "Software, Supervision, Writing -- original draft, Writing -- review and "
        "editing. Olivier Bernard: Conceptualization, Methodology, Formal analysis, "
        "Investigation, Software, Validation, Writing -- review and editing. "
        "Yichen Liu: Software, Methodology, Formal analysis, Investigation, "
        "Validation, Writing -- review and editing."
    )
    assert expected_credit in article
    assert article.count("Project administration") == 1
    assert "The Wahab cluster is supported in part by National Science Foundation" in article
    assert "CNS-1828593" in article
    assert "[ZENODO DOI TO BE INSERTED BEFORE SUBMISSION]" in article
    assert "MIT License" in article
    assert "Declaration of generative AI" not in article
    assert "OpenAI ChatGPT" not in article

    comparison = _read("sections/09_comparison.tex")
    assert "physics-based test equation" in comparison
    assert "finite-step CPTP" in comparison

    public_prose = "\n".join(
        (article, _read("Cover_Letter_JCP.tex"), _read("Supplementary_Material_JCP.tex"))
    )
    assert not re.search(r"(?<![A-Za-z])JCP(?![A-Za-z])", public_prose)
