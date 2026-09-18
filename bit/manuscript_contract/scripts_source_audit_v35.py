from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TEX_FILES = sorted(
    p for p in ROOT.rglob("*.tex")
    if not any(part.startswith(".") for part in p.relative_to(ROOT).parts)
)
TEXTS = {str(p.relative_to(ROOT)): p.read_text(errors="replace") for p in TEX_FILES}
ALL_TEXT = "\n".join(TEXTS.values())


def without_comments(text: str) -> str:
    lines: list[str] = []
    for line in text.splitlines():
        cut = len(line)
        for idx, char in enumerate(line):
            if char == "%" and (idx == 0 or line[idx - 1] != "\\"):
                cut = idx
                break
        lines.append(line[:cut])
    return "\n".join(lines)


ACTIVE_TEXT = without_comments(ALL_TEXT)
labels: list[tuple[str, str]] = []
refs: list[tuple[str, str]] = []
cites: list[tuple[str, str]] = []
for rel, text in TEXTS.items():
    active = without_comments(text)
    labels.extend((m.group(1), rel) for m in re.finditer(r"\\label\{([^}]+)\}", active))
    refs.extend((m.group(1), rel) for m in re.finditer(r"\\(?:ref|eqref|pageref|autoref)\{([^}]+)\}", active))
    for match in re.finditer(r"\\cite\w*\{([^}]+)\}", active):
        cites.extend((key.strip(), rel) for key in match.group(1).split(",") if key.strip())

bib_text = without_comments((ROOT / "references_jcp.tex").read_text(errors="replace"))
bibs = re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", bib_text)
label_names = [name for name, _ in labels]
cite_names = [name for name, _ in cites]
label_counts = Counter(label_names)
bib_counts = Counter(bibs)

article = (ROOT / "jcp_article.tex").read_text(errors="replace")
abstract_match = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", article, re.S)
if abstract_match is None:
    raise RuntimeError("abstract environment not found")
abstract = abstract_match.group(1)


def texcount_words(fragment: str) -> tuple[int, str]:
    texcount = shutil.which("texcount")
    if texcount:
        with tempfile.NamedTemporaryFile("w", suffix=".tex", delete=False) as handle:
            handle.write("\\documentclass{article}\n\\begin{document}\n")
            handle.write(fragment)
            handle.write("\n\\end{document}\n")
            temp_name = handle.name
        try:
            run = subprocess.run(
                [texcount, "-sum", "-1", temp_name],
                check=True,
                text=True,
                capture_output=True,
            )
            return int(run.stdout.strip()), "texcount"
        except (OSError, subprocess.CalledProcessError, ValueError):
            pass
        finally:
            Path(temp_name).unlink(missing_ok=True)
    plain = re.sub(r"\\[A-Za-z]+(?:\[[^\]]*\])?\{([^{}]*)\}", r"\1", fragment)
    plain = re.sub(r"\\[A-Za-z]+|[$~{}_^\\]", " ", plain)
    words = re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)*", plain)
    return len(words), "fallback-parser"


abstract_word_count, abstract_counter = texcount_words(abstract)
abstract_word_limit = 250

highlight_lines = [
    line.removeprefix("-").strip()
    for line in (ROOT / "Highlights.txt").read_text(errors="replace").splitlines()
    if line.strip().startswith("-")
]
highlight_lengths = [len(line) for line in highlight_lines]

missing_graphics: list[dict[str, str]] = []
for rel, text in TEXTS.items():
    source_path = ROOT / rel
    for match in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", without_comments(text)):
        requested = match.group(1)
        roots = [source_path.parent / requested, ROOT / requested, ROOT / "figures" / requested]
        candidates: list[Path] = []
        for root in roots:
            candidates.append(root)
            if not root.suffix:
                candidates.extend(root.with_suffix(ext) for ext in (".pdf", ".png", ".eps", ".jpg", ".jpeg"))
        if not any(path.exists() for path in candidates):
            missing_graphics.append({"source": rel, "target": requested})

missing_inputs: list[dict[str, str]] = []
for rel, text in TEXTS.items():
    source_path = ROOT / rel
    for match in re.finditer(r"\\(?:input|include)\{([^}]+)\}", without_comments(text)):
        requested = match.group(1)
        roots = [source_path.parent / requested, ROOT / requested]
        candidates: list[Path] = []
        for root in roots:
            candidates.append(root)
            if not root.suffix:
                candidates.append(root.with_suffix(".tex"))
        if not any(path.exists() for path in candidates):
            missing_inputs.append({"source": rel, "target": requested})

active_visible_paths = [
    "jcp_article.tex",
    "Cover_Letter_JCP.tex",
    "Supplementary_Material_JCP.tex",
    "Figure_Captions.txt",
    "Highlights.txt",
    "Title_Page.txt",
    "README.txt",
]
visible = "\n".join((ROOT / path).read_text(errors="replace") for path in active_visible_paths)

ambiguous_terms = sorted(
    set(
        match.group(0)
        for match in re.finditer(r"rotating slice|(?<!-frame )rotating RK4", visible, re.I)
    )
)
legacy_active_names = sorted(
    set(
        match.group(0)
        for match in re.finditer(
            r"benchmark_macros_v33|adaptive_benchmark_table_v33|Fig3_certified_guard|Fig4_tolerance_sweep|Fig5_robustness|Fig5_operational_refinement|Fig2_candidate_maps|Fig3_robustness|Fig6_frequency_bands|Supplemental_Material",
            ACTIVE_TEXT,
        )
    )
)

final_logs = [
    "Manuscript_JCP.log",
    "Manuscript_JCP_clean.log",
    "Supplementary_Material_JCP.log",
    "Cover_Letter_JCP.log",
]
log_failures: dict[str, list[str]] = {}
log_pattern = re.compile(
    r"undefined|multiply defined|destination with the same identifier|There were undefined|Missing character|Overfull",
    re.I,
)
for log_name in final_logs:
    log_path = ROOT / log_name
    if not log_path.exists():
        continue
    hits = [line for line in log_path.read_text(errors="replace").splitlines() if log_pattern.search(line)]
    if hits:
        log_failures[log_name] = hits

verification_text = (ROOT / "sections/09_comparison.tex").read_text(errors="replace")

result = {
    "version": "3.5.0",
    "tex_files": len(TEX_FILES),
    "labels": len(label_names),
    "unique_labels": len(set(label_names)),
    "duplicate_labels": sorted(name for name, count in label_counts.items() if count > 1),
    "reference_occurrences": len(refs),
    "undefined_references": sorted({name for name, _ in refs if name not in set(label_names)}),
    "citation_occurrences": len(cite_names),
    "unique_citations": len(set(cite_names)),
    "bibliography_entries": len(bibs),
    "undefined_citations": sorted(set(cite_names) - set(bibs)),
    "uncited_bibliography_entries": sorted(set(bibs) - set(cite_names)),
    "duplicate_bibliography_keys": sorted(name for name, count in bib_counts.items() if count > 1),
    "abstract_word_count": abstract_word_count,
    "abstract_counter": abstract_counter,
    "abstract_word_limit": abstract_word_limit,
    "highlights": [
        {"text": text, "characters_including_spaces": length, "within_85": length <= 85}
        for text, length in zip(highlight_lines, highlight_lengths, strict=True)
    ],
    "missing_graphics": missing_graphics,
    "missing_inputs": missing_inputs,
    "ambiguous_rotating_terminology": ambiguous_terms,
    "legacy_active_asset_names": legacy_active_names,
    "scope_contract": {
        "pointwise_transverse_values_present": all(token in visible for token in ("0.05", "0.10", "0.20")),
        "no_uniform_transverse_claim_present": bool(re.search(r"no uniform|does not certify the full interval", visible, re.I)),
        "controller_safety_scope_present": bool(re.search(r"fail-closed safety|not runtime efficiency|adaptive safety", visible, re.I)),
        "binary_recertification_failure_present": ("EXECUTED_BINARY_RECERTIFICATION_FAILED" in visible or "EXECUTED\\_BINARY\\_RECERTIFICATION\\_FAILED" in visible),
        "work_proxy_exclusion_present": bool(
            re.search(
                r"work proxy.*(?:exclud|reported separately|counts only)|excluded from (?:that|the) proxy",
                visible,
                re.I | re.S,
            )
        ),
    },
    "disclosures": {
        "no_ai_declaration_in_manuscript": (
            "Declaration of generative AI" not in article
            and "OpenAI ChatGPT" not in article
        ),
        "no_duplicate_ai_note_in_verification": "OpenAI ChatGPT" not in verification_text,
        "graphical_abstract_reproducible_provenance": (
            "regenerated deterministically with the versioned scripts" in verification_text
            and (
                "generated deterministically in Matplotlib" in visible
                or "generated deterministically by Matplotlib" in visible
            )
        ),
    },
    "final_latex_log_failures": log_failures,
}

critical = {
    "duplicate_labels": result["duplicate_labels"],
    "undefined_references": result["undefined_references"],
    "undefined_citations": result["undefined_citations"],
    "uncited_bibliography_entries": result["uncited_bibliography_entries"],
    "duplicate_bibliography_keys": result["duplicate_bibliography_keys"],
    "missing_graphics": result["missing_graphics"],
    "missing_inputs": result["missing_inputs"],
    "ambiguous_rotating_terminology": result["ambiguous_rotating_terminology"],
    "legacy_active_asset_names": result["legacy_active_asset_names"],
    "final_latex_log_failures": result["final_latex_log_failures"],
}
result["critical_failures"] = {key: value for key, value in critical.items() if value}
result["pass"] = (
    not result["critical_failures"]
    and 150 <= result["abstract_word_count"] <= abstract_word_limit
    and 3 <= len(highlight_lines) <= 5
    and all(length <= 85 for length in highlight_lengths)
    and all(result["scope_contract"].values())
    and all(result["disclosures"].values())
)

output = ROOT / "SOURCE_AUDIT_v3_5.json"
output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["pass"] else 1)
