#!/usr/bin/env python3
"""Deterministic BIT Numerical Mathematics source and package auditor (Task 16).

    python scripts/audit_bit_source.py --manuscript ../manuscript [--submission ../submission]
        [--build-dir ../build/main_build] [--code-dir .] --status PRE_SUBMISSION|FINAL [--json out.json]

Every rule is reported as PASS, FAIL, or SKIP with machine-readable detail.  The exit code is 0 only if no
rule failed.  A PASS is a mechanical check of the source/package; it is not an author approval.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

DELIBERATE_MARKERS = [
    "[PERSISTENT DOI OR REPOSITORY URL OF THE ARCHIVED RELEASE TO BE INSERTED BEFORE SUBMISSION]",
    "[THIS STATEMENT MUST BE REVIEWED AND APPROVED BY ALL AUTHORS BEFORE SUBMISSION.]",
    "[NOT-UNDER-CONSIDERATION-ELSEWHERE STATEMENT TO BE CONFIRMED BY THE AUTHORS BEFORE SUBMISSION]",
    "[DATE OF SUBMISSION]",
]
MARKER_WORDS = re.compile(r"\bTODO\b|\bTBD\b|PLACEHOLDER|TO BE INSERTED|ZENODO DOI TO BE|\bINSERT\b|\bXXX\b")
BRACKET_MARKER = re.compile(r"\[[A-Z][A-Z0-9 \-.,'/]{10,}\]")
PENDING_FINAL_STATUS = re.compile(
    r"final author approval(?: of this BIT revision)? is pending|"
    r"approval of (?:the )?final BIT revision remains pending|"
    r"new public deposit is pending|persistent location remains to be supplied|"
    r"prepared for public release|archive has not been deposited", re.I)
JCP_RESIDUE = re.compile(r"elsarticle|Elsevier|graphical abstract|\bHighlights\b|Journal of Computational Physics|submission to JCP", re.I)
STALE_TERMS = [
    "Supplemental Material", "Fig6_frequency_bands", "candidate record", "structural pass", "typed recipe",
    "PI history", "work proxy", "binary step", "Complete-positivity topology of Runge",
]
JCP_ONLY_FILES = re.compile(
    r"^(Highlights\.txt|Abstract\.txt|Title_Page\.txt|Figure_Captions\.txt|BUILD_INSTRUCTIONS\.txt|README\.txt|"
    r"Cover_Letter_JCP\.tex|Supplementary_Material_JCP\.tex|Manuscript_JCP(_clean)?\.tex|jcp_article\.tex|jcp_preamble\.tex|"
    r"references_jcp\.tex|scripts_source_audit_v35\.py|SOURCE_AUDIT_v3_5\.json|SOURCE_MANIFEST\.txt|SOURCE_SHA256SUMS\.txt|"
    r"Graphical_Abstract.*|elsarticle.*)$", re.I)
REQUIRED_DECLARATIONS = [
    "Funding", "Competing interests", "Ethics approval and consent to participate", "Consent for publication",
    "Data availability", "Code availability", "Author contributions", "Use of generative artificial intelligence",
]
TEXT_EXT = {".tex", ".bib", ".txt", ".md", ".cff", ".toml", ".csv", ".json", ".yaml", ".yml", ".cls", ".bst"}
JUNK_EXT = {".aux", ".log", ".out", ".fls", ".fdb_latexmk", ".synctex.gz", ".blg", ".toc", ".bak", ".pyc", ".DS_Store"}
SUBMISSION_FILES = [
    "Manuscript_BIT_review.pdf", "Manuscript_BIT_clean.pdf", "Cover_Letter_BIT.pdf", "Cover_Letter_BIT.txt", "ESM_1.pdf",
    "BIT_Manuscript_Source.zip", "BIT_ESM_1_Source.zip", "BIT_Code_and_Data.zip", "BIT_Submission_Metadata.md",
    "BIT_Declarations.md", "BIT_README_FIRST.md", "BIT_SHA256.txt",
]


# ---------------------------------------------------------------------------
# TeX helpers
# ---------------------------------------------------------------------------
def strip_comments(s: str) -> str:
    return re.sub(r"(?<!\\)%.*", "", s)


def braced_arg(s: str, start: int) -> tuple[str, int]:
    """Return the content of the {...} group starting at index `start` (which must be '{')."""
    assert s[start] == "{"
    depth, i = 0, start
    while i < len(s):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                return s[start + 1:i], i + 1
        i += 1
    raise ValueError("unbalanced braces")


def command_arg(s: str, name: str, optional=False) -> str | None:
    m = re.search(r"\\" + re.escape(name) + (r"\[[^\]]*\]" if optional else "") + r"\s*\{", s)
    if not m:
        return None
    return braced_arg(s, m.end() - 1)[0]


def detex(s: str) -> str:
    s = strip_comments(s)
    s = re.sub(r"\\(cite[tp]?|ref|eqref|label|autoref)\{[^}]*\}", " ", s)
    s = re.sub(r"\\href\{[^}]*\}\{([^}]*)\}", r"\1", s)
    for _ in range(3):
        s = re.sub(r"\\(textit|textbf|emph|texttt|textsc|textnormal|mbox)\{([^{}]*)\}", r"\2", s)
    s = re.sub(r"\$[^$]*\$", " EQN ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", s)
    s = s.replace("``", '"').replace("''", '"').replace("--", "-").replace("~", " ")
    s = re.sub(r"[{}]", "", s)
    return s


def word_count(s: str) -> int:
    return len([w for w in detex(s).split() if re.search(r"[A-Za-z0-9]", w)])


def read(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="replace")


def input_files(main: Path) -> list[Path]:
    files, seen = [main], set()

    def walk(p: Path):
        for m in re.finditer(r"\\input\{([^}]+)\}", strip_comments(read(p))):
            name = m.group(1)
            q = p.parent / (name if name.endswith(".tex") else name + ".tex")
            if q.exists() and q not in seen:
                seen.add(q)
                files.append(q)
                walk(q)

    walk(main)
    return files


def cite_keys(text: str) -> set[str]:
    keys = set()
    for m in re.finditer(r"\\cite[tp]?\*?(?:\[[^\]]*\])?\{([^}]*)\}", strip_comments(text)):
        keys.update(k.strip() for k in m.group(1).split(",") if k.strip())
    return keys


def bib_keys(bib: str) -> set[str]:
    return {m.group(1).strip() for m in re.finditer(r"^@\w+\{([^,]+),", bib, re.M)}


def labels(text: str) -> set[str]:
    return {m.group(1) for m in re.finditer(r"\\label\{([^}]+)\}", strip_comments(text))}


def refs(text: str) -> set[str]:
    out = set()
    for m in re.finditer(r"\\(?:ref|eqref|autoref|pageref)\{([^}]+)\}", strip_comments(text)):
        out.update(k.strip() for k in m.group(1).split(","))
    return out


def authors(main_text: str) -> list[dict]:
    out = []
    for m in re.finditer(r"\\author(\*?)\[([\d, ]+)\]\{\\fnm\{([^}]*)\}\s*\\sur\{([^}]*)\}\}\s*\\email\{([^}]*)\}", main_text):
        out.append({"corresponding": m.group(1) == "*", "affils": [a.strip() for a in m.group(2).split(",")],
                    "given": m.group(3).strip(), "family": m.group(4).strip(), "email": m.group(5).strip()})
    return out


# ---------------------------------------------------------------------------
# audit
# ---------------------------------------------------------------------------
class Audit:
    def __init__(self, manuscript: Path, status: str, submission: Path | None, build_dir: Path | None, code_dir: Path | None):
        self.m = manuscript
        self.status = status
        self.sub = submission
        self.build = build_dir
        self.code = code_dir
        self.rules: list[dict] = []
        self.main = manuscript / "main.tex"
        self.main_text = read(self.main) if self.main.exists() else ""
        self.article_files = input_files(self.main) if self.main.exists() else []
        self.article_text = "\n".join(read(p) for p in self.article_files)
        self.esm = manuscript / "ESM_1.tex"
        self.esm_text = read(self.esm) if self.esm.exists() else ""
        self.cover = manuscript / "cover_letter.tex"
        self.cover_text = read(self.cover) if self.cover.exists() else ""
        self.bib = manuscript / "references.bib"
        self.bib_text = read(self.bib) if self.bib.exists() else ""

    def add(self, rule_id: str, description: str, ok: bool | None, detail=None):
        self.rules.append({"id": rule_id, "description": description,
                           "status": "SKIP" if ok is None else ("PASS" if ok else "FAIL"), "detail": detail})

    # --- front matter ------------------------------------------------------
    def rule_abstract(self):
        abstract = command_arg(self.main_text, "abstract")
        n = word_count(abstract) if abstract else 0
        self.add("abstract_word_count", "abstract has 150-250 words (R8)", 150 <= n <= 250, {"words": n})
        caps = sorted({w for w in re.findall(r"\b[A-Z]{2,}\b", detex(abstract or "")) if w != "EQN"})
        self.add("abstract_no_undefined_abbreviations", "abstract contains no all-capital abbreviations (R8)", not caps, {"abbreviations": caps})

    def rule_keywords(self):
        kw = command_arg(self.main_text, "keywords")
        items = [k.strip() for k in (kw or "").split(",") if k.strip()]
        self.add("keyword_count", "4-6 keywords (R9)", 4 <= len(items) <= 6, {"keywords": items})

    def rule_msc(self):
        msc = command_arg(self.main_text, "pacs", optional=True)
        codes = re.findall(r"\b\d\d[A-Z]\d\d\b", msc or "")
        self.add("msc_codes_present", "MSC 2020 codes present in \\pacs[MSC Classification] (R11)", bool(codes), {"codes": codes})

    def rule_title_and_authors(self):
        title = command_arg(self.main_text, "title", optional=True)
        au = authors(self.main_text)
        ok = bool(title) and len(au) >= 1 and sum(a["corresponding"] for a in au) == 1 and all("@" in a["email"] for a in au)
        affils = set(re.findall(r"\\affil\*?\[(\d+)\]", self.main_text))
        used = {x for a in au for x in a["affils"]}
        ok = ok and used == affils and affils == {str(i) for i in range(1, len(affils) + 1)}
        self.add("title_page_complete", "title, authors with e-mails, exactly one corresponding author, affiliation numbering consecutive and used (R6)", ok,
                 {"title": title, "authors": au, "affiliations": sorted(affils)})

    def rule_authors_sync(self):
        au = authors(self.main_text)
        names = [f"{a['given']} {a['family']}" for a in au]
        problems = []
        esm_norm = re.sub(r"\s+", " ", detex(self.esm_text))
        for n in names:
            if n not in esm_norm and n.replace("G. Blake", "G.~Blake") not in self.esm_text:
                problems.append(f"{n} missing from ESM_1.tex")
        corr = [a for a in au if a["corresponding"]]
        if corr and corr[0]["email"] not in self.esm_text:
            problems.append("corresponding e-mail missing from ESM_1.tex")
        if self.code and (self.code / "CITATION.cff").exists():
            cff = read(self.code / "CITATION.cff")
            for a in au:
                if a["family"] not in cff:
                    problems.append(f"{a['family']} missing from CITATION.cff")
        if self.sub and (self.sub / "BIT_Submission_Metadata.md").exists():
            meta = read(self.sub / "BIT_Submission_Metadata.md")
            for a in au:
                if a["email"] not in meta or a["family"] not in meta:
                    problems.append(f"{a['family']}/{a['email']} missing from BIT_Submission_Metadata.md")
        self.add("authors_synchronized", "author names/e-mails identical in article, Online Resource, CITATION.cff and metadata sheet", not problems, {"authors": names, "problems": problems})

    # --- body --------------------------------------------------------------
    def rule_refs(self):
        lab, rf = labels(self.article_text), refs(self.article_text)
        missing = sorted(rf - lab)
        self.add("article_refs_resolve", "every \\ref in the article has a \\label", not missing, {"missing": missing, "labels": len(lab)})
        if self.esm_text:
            lab2, rf2 = labels(self.esm_text), refs(self.esm_text)
            for extra in ("article_numbers", "nine_candidate_certificate", "extrapolation_table", "tolerance_sweep_table", "guard_telemetry_table"):
                p = self.m / f"{extra}.tex"
                if p.exists():
                    lab2 |= labels(read(p))
                    rf2 |= refs(read(p))
            missing2 = sorted(rf2 - lab2)
            self.add("esm_refs_resolve", "every \\ref in the Online Resource has a \\label", not missing2, {"missing": missing2})

    def rule_citations(self):
        keys = cite_keys(self.article_text)
        bk = bib_keys(self.bib_text)
        missing = sorted(keys - bk)
        unused = sorted(bk - keys)
        self.add("citation_keys_in_bibliography", "every cited key exists in references.bib", not missing, {"cited": len(keys), "missing": missing, "uncited_entries": unused})
        self.add("bibliography_has_no_uncited_entries", "references.bib contains only cited works (R17)", not unused, {"uncited_entries": unused})
        dois = re.findall(r"doi\s*=\s*\{([^}]+)\}", self.bib_text)
        bad = [d for d in dois if d.lower().startswith("http")]
        self.add("doi_fields_bare", "DOI fields are bare identifiers so the style renders https://doi.org links (R17)", not bad, {"dois": len(dois), "bad": bad})

    def rule_headings(self):
        deep = re.findall(r"\\subparagraph\{", self.article_text)
        levels = {"section": len(re.findall(r"\\section\{", self.article_text)), "subsection": len(re.findall(r"\\subsection\{", self.article_text)),
                  "subsubsection": len(re.findall(r"\\subsubsection\{", self.article_text)), "paragraph_runin": len(re.findall(r"\\paragraph\{", self.article_text))}
        self.add("heading_levels", "numbered decimal headings use at most three levels (R13); run-in \\paragraph titles are unnumbered", not deep, levels)

    def rule_captions(self):
        bad = []
        for m in re.finditer(r"\\caption\s*\{", self.article_text):
            cap, _ = braced_arg(self.article_text, m.end() - 1)
            text = detex(cap).strip()
            if text.endswith("."):
                bad.append(text[-60:])
        self.add("caption_style", "captions do not end with a period (Springer figure/table caption style)", not bad, {"offending": bad})

    def rule_figures(self):
        inc = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", self.article_text)
        names = [Path(i).name for i in inc]
        expected = [f"Fig{i}.eps" for i in range(1, len(names) + 1)]
        exist = [n for n in names if (self.m / n).exists()]
        problems = []
        if names != expected:
            problems.append(f"included figures {names} are not the sequence {expected}")
        if len(exist) != len(names):
            problems.append("missing figure files")
        fig_labels = re.findall(r"\\begin\{figure\}.*?\\label\{([^}]+)\}.*?\\end\{figure\}", self.article_text, re.S)
        rf = refs(self.article_text)
        uncited = [l for l in fig_labels if l not in rf]
        if uncited:
            problems.append(f"figures not referenced in the text: {uncited}")
        if len(fig_labels) != len(names):
            problems.append("every figure environment must contain one \\includegraphics and one \\label")
        self.add("figures_sequential_and_cited", "figures are Fig1.eps...FigN.eps, present, and cited in order (R19)", not problems, {"figures": names, "problems": problems})
        if self.esm_text:
            inc2 = [Path(i).name for i in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", self.esm_text)]
            exp2 = [f"FigS{i}.eps" for i in range(1, len(inc2) + 1)]
            ok2 = inc2 == exp2 and all((self.m / n).exists() for n in inc2)
            self.add("esm_figures_sequential", "Online Resource figures are FigS1.eps...FigSM.eps and present", ok2, {"figures": inc2})

    def rule_declarations(self):
        heads = re.findall(r"\\bmhead\{([^}]+)\}", self.article_text)
        missing = [h for h in REQUIRED_DECLARATIONS if h not in heads]
        self.add("required_declarations_present", "Declarations section contains every required heading (R10, R22, R24)", not missing, {"headings": heads, "missing": missing})

        def section_words(head: str) -> int:
            m = re.search(r"\\bmhead\{" + re.escape(head) + r"\}(.*?)(?=\\bmhead|\\section|\\bibliography|\Z)", self.article_text, re.S)
            return word_count(m.group(1)) if m else 0

        self.add("ai_disclosure_present", "generative-AI statement present with substantive content (R7)", section_words("Use of generative artificial intelligence") >= 40,
                 {"words": section_words("Use of generative artificial intelligence")})
        self.add("data_availability_present", "Data Availability Statement present (R24)", section_words("Data availability") >= 15, {"words": section_words("Data availability")})

    def rule_online_resource(self):
        ok = "Online Resource~1" in self.article_text or "Online Resource 1" in self.article_text
        ok = ok and ("ESM\\_1.pdf" in self.article_text or "ESM_1.pdf" in self.article_text) and self.esm.exists()
        title = command_arg(self.main_text, "title", optional=True) or ""
        norm = lambda s: re.sub(r"\s+", " ", detex(s).replace("-", "")).strip().lower()
        ok = ok and "Online Resource 1" in detex(self.esm_text) and norm(title)[:40] in norm(self.esm_text)
        self.add("online_resource_synchronized", "article cites Online Resource 1 as ESM_1.pdf; ESM_1.tex exists and carries the article title (R20)", ok, {"esm_exists": self.esm.exists()})

    def rule_cover_letter(self):
        m = re.search(r"\\begin\{letterbody\}(.*?)\\end\{letterbody\}", strip_comments(self.cover_text), re.S)
        if not m:
            self.add("cover_letter_word_count", "cover letter body at most 150 words (R1)", False, {"error": "letterbody environment not found"})
            return
        body = m.group(1)
        markers = [mk for mk in DELIBERATE_MARKERS if mk in body]
        for mk in markers:
            body = body.replace(mk, "")
        n = word_count(body)
        estimate = n + (7 if markers else 0)  # a seven-word statement replaces the marker at final build
        self.add("cover_letter_word_count", "cover letter body at most 150 words (R1)", estimate <= 150,
                 {"words_without_markers": n, "markers": markers, "estimated_final_words": estimate})

    # --- residues and markers -----------------------------------------------
    def text_files(self) -> list[Path]:
        files = [p for p in self.m.iterdir() if p.is_file() and p.suffix.lower() in TEXT_EXT]
        if self.sub and self.sub.exists():
            files += [p for p in self.sub.iterdir() if p.is_file() and p.suffix.lower() in TEXT_EXT]
        return sorted(files)

    def rule_markers(self):
        deliberate, unknown, pending = [], [], []
        for p in self.text_files():
            if p.name in ("sn-jnl.cls", "sn-mathphys-num.bst"):
                continue
            # normalise line breaks and Markdown block-quote prefixes so a marker wrapped over two lines is
            # still recognised as the deliberate marker (and cannot hide from a FINAL check)
            t = re.sub(r"\s+", " ", re.sub(r"\n>\s?", "\n", read(p)))
            pending.extend(f"{p.name}: {m.group(0)}" for m in PENDING_FINAL_STATUS.finditer(t))
            for mk in DELIBERATE_MARKERS:
                if mk in t:
                    deliberate.append(f"{p.name}: {mk}")
            t2 = t
            for mk in DELIBERATE_MARKERS:
                t2 = t2.replace(mk, "")
            for m in MARKER_WORDS.finditer(t2):
                unknown.append(f"{p.name}: {m.group(0)}")
            for m in BRACKET_MARKER.finditer(t2):
                unknown.append(f"{p.name}: {m.group(0)}")
        if self.status == "FINAL":
            ok = not deliberate and not unknown and not pending
        else:
            ok = not unknown
        self.add("unresolved_markers", "no TODO/TBD/PLACEHOLDER/INSERT markers; deliberate pre-submission markers only while PRE_SUBMISSION", ok,
                 {"deliberate_markers": deliberate, "unknown_markers": unknown, "pending_final_status": pending})
        self.add("package_status_consistent", "status FINAL requires zero markers; PRE_SUBMISSION whenever a human gate remains", ok if self.status == "FINAL" else True,
                 {"status": self.status, "open_markers": len(deliberate), "pending_final_status": pending})

    def rule_residue(self):
        hits = []
        for p in self.text_files():
            if p.name in ("sn-jnl.cls", "sn-mathphys-num.bst", "references.bib") or p.suffix.lower() in (".cls", ".bst"):
                continue
            t = read(p)
            for m in JCP_RESIDUE.finditer(t):
                hits.append(f"{p.name}: {m.group(0)}")
            for term in STALE_TERMS:
                if term in t:
                    hits.append(f"{p.name}: {term}")
        if self.bib.exists():
            bt = self.bib_text
            for m in re.finditer(r"elsarticle|graphical abstract|\bHighlights\b", bt, re.I):
                hits.append(f"references.bib: {m.group(0)}")
        self.add("residue_scan", "no Elsevier/JCP boilerplate, former title, or retired jargon in manuscript and submission text files (plan 5.3)", not hits, {"hits": hits})

    def rule_jcp_files(self):
        bad = sorted(p.name for p in self.m.iterdir() if JCP_ONLY_FILES.match(p.name) or p.is_dir())
        self.add("no_jcp_boilerplate_files", "no JCP/Elsevier-only files or subdirectories in the manuscript directory (R15, plan 5.3)", not bad, {"offending": bad})

    # --- build products ------------------------------------------------------
    def rule_build(self):
        if not self.build or not (self.build / "main.log").exists():
            self.add("build_log_clean", "latexmk log has no undefined references/citations, missing characters, or overfull boxes", None, {"reason": "no build directory"})
            self.add("esm_article_numbers_synchronized", "article_numbers.tex matches the compiled article's main.aux", None, {"reason": "no build directory"})
            return
        log = read(self.build / "main.log")
        problems = []
        for pat in (r"Citation `[^']*' on page \d+ undefined", r"Reference `[^']*' on page \d+ undefined", r"There were undefined", r"Missing character", r"^! ", r"Overfull \\hbox"):
            for m in re.finditer(pat, log, re.M):
                problems.append(m.group(0)[:80])
        self.add("build_log_clean", "latexmk log has no undefined references/citations, missing characters, errors, or overfull boxes", not problems, {"problems": problems[:20]})
        aux = self.build / "main.aux"
        an = self.m / "article_numbers.tex"
        if aux.exists() and an.exists():
            sys.path.insert(0, str(Path(__file__).resolve().parent))
            import sync_article_numbers as san

            expected = san.render(read(aux))
            self.add("esm_article_numbers_synchronized", "article_numbers.tex matches the compiled article's main.aux", read(an) == expected, {})
        else:
            self.add("esm_article_numbers_synchronized", "article_numbers.tex matches the compiled article's main.aux", None, {"reason": "aux or article_numbers.tex missing"})
        bbl = self.build / "main.bbl"
        if bbl.exists():
            bt = read(bbl)
            dois = re.findall(r"doi\s*=\s*\{([^}]+)\}", self.bib_text)
            missing = [d for d in dois if not any(f"{macro}{{{d}}}" in bt for macro in (r"\doiurl", r"\doi", r"\url{https://doi.org/"))]
            txt = self.build / "main.txt"
            rendered = len(re.findall(r"https://doi\.", read(txt))) if txt.exists() else None
            self.add("doi_links_rendered", "every DOI is passed to the style's DOI macro, which renders it as a full https://doi.org link (R17)", not missing,
                     {"missing": missing, "dois": len(dois), "doi_links_in_pdf_text": rendered})

    # --- submission package --------------------------------------------------
    def rule_submission(self):
        if not self.sub or not self.sub.exists():
            self.add("submission_files_present", "all required submission artifacts exist", None, {"reason": "no submission directory"})
            return
        present = {n: (self.sub / n).exists() for n in SUBMISSION_FILES}
        outer = sorted(self.sub.glob("BIT_Submission_Package_*.zip"))
        self.add("submission_files_present", "all required submission artifacts exist (plan Section 3)", all(present.values()) and len(outer) == 1 and (self.sub / (outer[0].name + ".sha256")).exists() if outer else False,
                 {"present": present, "outer_packages": [p.name for p in outer]})
        if outer:
            ok = outer[0].name.endswith(f"_{self.status}.zip")
            readme = read(self.sub / "BIT_README_FIRST.md") if (self.sub / "BIT_README_FIRST.md").exists() else ""
            self.add("package_status_labelled", "outer package name and README carry the package status", ok and self.status in readme, {"outer": outer[0].name})
        # PDFs
        try:
            from pypdf import PdfReader
        except ImportError:  # pragma: no cover
            PdfReader = None
        pdf_info = {}
        pdf_ok = True
        for name in ("Manuscript_BIT_review.pdf", "Manuscript_BIT_clean.pdf", "Cover_Letter_BIT.pdf", "ESM_1.pdf"):
            p = self.sub / name
            if not p.exists() or PdfReader is None:
                pdf_ok = False
                continue
            r = PdfReader(str(p))
            text = "".join((pg.extract_text() or "") for pg in r.pages[:3])
            fonts_ok, fonts = fonts_embedded(p)
            info = {"pages": len(r.pages), "encrypted": bool(r.is_encrypted), "searchable": len(text) > 200, "fonts_embedded": fonts_ok, "fonts": fonts,
                    "title": (r.metadata or {}).get("/Title")}
            pdf_info[name] = info
            if info["encrypted"] or not info["searchable"] or fonts_ok is not True:
                pdf_ok = False
            if name == "Cover_Letter_BIT.pdf" and info["pages"] != 1:
                pdf_ok = False
        self.add("pdfs_readable_searchable_fonts_embedded", "every PDF is unencrypted, searchable, has embedded outline fonts (no bitmap Type 3 text fonts); cover letter is one page (plan 5.5)", pdf_ok, pdf_info)
        if "Manuscript_BIT_review.pdf" in pdf_info and "Manuscript_BIT_clean.pdf" in pdf_info:
            same_pages = pdf_info["Manuscript_BIT_review.pdf"]["pages"] == pdf_info["Manuscript_BIT_clean.pdf"]["pages"]
            self.add("review_and_clean_from_same_source", "review and clean PDFs have the same page count (line numbers only)", same_pages,
                     {"review": pdf_info["Manuscript_BIT_review.pdf"]["pages"], "clean": pdf_info["Manuscript_BIT_clean.pdf"]["pages"]})
        # ZIPs
        for name, need_main in (("BIT_Manuscript_Source.zip", True), ("BIT_ESM_1_Source.zip", True), ("BIT_Code_and_Data.zip", False)):
            p = self.sub / name
            if not p.exists():
                continue
            self.add(f"zip_safe_{name}", f"{name}: CRC ok, safe paths, {'flat and exactly one TeX entry point' if need_main else 'no environments/caches/temporaries'}", *zip_rule(p, need_main))
        # checksums
        man = self.sub / "BIT_SHA256.txt"
        if man.exists():
            import hashlib

            bad, listed = [], set()
            for line in read(man).splitlines():
                if not line.strip():
                    continue
                h, _, fname = line.partition("  ")
                fname = fname.strip().lstrip("*")
                listed.add(fname)
                q = self.sub / fname
                if not q.exists() or hashlib.sha256(q.read_bytes()).hexdigest() != h.strip():
                    bad.append(fname)
            expected = {n for n in SUBMISSION_FILES if n != "BIT_SHA256.txt" and (self.sub / n).exists()}
            self.add("checksums_match", "BIT_SHA256.txt covers every payload file and matches the bytes on disk (plan 5.6)", not bad and listed == expected,
                     {"mismatched": bad, "listed": sorted(listed), "not_listed": sorted(expected - listed), "extra": sorted(listed - expected)})

    def run(self) -> dict:
        if not self.main.exists():
            self.add("main_tex_present", "manuscript/main.tex exists", False, {"path": str(self.main)})
        else:
            self.add("main_tex_present", "manuscript/main.tex exists", True, {"inputs": [p.name for p in self.article_files]})
            self.rule_abstract()
            self.rule_keywords()
            self.rule_msc()
            self.rule_title_and_authors()
            self.rule_authors_sync()
            self.rule_refs()
            self.rule_citations()
            self.rule_headings()
            self.rule_captions()
            self.rule_figures()
            self.rule_declarations()
            self.rule_online_resource()
            self.rule_cover_letter()
            self.rule_markers()
            self.rule_residue()
            self.rule_jcp_files()
            self.rule_build()
            self.rule_submission()
        counts = {s: sum(1 for r in self.rules if r["status"] == s) for s in ("PASS", "FAIL", "SKIP")}
        return {"manuscript": str(self.m), "status": self.status, "counts": counts, "rules": self.rules}


def parse_pdffonts(rows: list[str]):
    """Interpret `pdffonts` rows: every font embedded, and no bitmap Type 3 text font (pdfTeX names those F<n>)."""
    cols = [l.split() for l in rows if l.strip()]
    if not cols:
        return None, "no font rows"
    # columns: name type encoding emb sub uni object ID (the type may contain a space, so count from the end)
    emb = [c[-5] if len(c) >= 7 else "?" for c in cols]
    types = [" ".join(c[1:-6]) if len(c) >= 7 else "?" for c in cols]
    bitmap_text = [c[0] for c, ty in zip(cols, types) if ty == "Type 3" and re.fullmatch(r"F\d+", c[0])]
    ok = all(e == "yes" for e in emb) and not bitmap_text
    return ok, {"fonts": len(cols), "not_embedded": sum(1 for e in emb if e != "yes"), "type3_fonts": sum(1 for ty in types if ty == "Type 3"),
                "type3_text_fonts": bitmap_text}


def fonts_embedded(pdf: Path):
    exe = shutil.which("pdffonts")
    if not exe:
        return None, "pdffonts not available"
    out = subprocess.run([exe, str(pdf)], capture_output=True, text=True, errors="replace").stdout.splitlines()
    return parse_pdffonts(out[2:])


def zip_rule(p: Path, need_main: bool):
    problems = []
    # the reproducibility archive legitimately carries the historical reproduction logs (*.log) of the
    # sealed v3.5 release; LaTeX temporaries are still rejected everywhere
    junk = JUNK_EXT if need_main else JUNK_EXT - {".log", ".out"}
    with zipfile.ZipFile(p) as zf:
        bad = zf.testzip()
        if bad:
            problems.append(f"CRC error in {bad}")
        names = zf.namelist()
        for n in names:
            if n.startswith("/") or ".." in n.split("/") or re.match(r"^[A-Za-z]:", n):
                problems.append(f"unsafe path {n}")
            if any(part in ("__pycache__", ".venv", ".pytest_cache", ".git", ".matplotlib-cache") for part in n.split("/")):
                problems.append(f"environment/cache member {n}")
            if n.lower().endswith(tuple(junk)) and not n.endswith("main.bbl"):
                problems.append(f"temporary file {n}")
            if need_main and "/" in n:
                problems.append(f"nested member {n}")
        if need_main:
            entry = [n for n in names if n.endswith(".tex") and b"\\documentclass" in zf.read(n)]
            if len(entry) != 1:
                problems.append(f"expected exactly one TeX entry point, found {entry}")
    return (not problems, {"members": len(names), "problems": problems[:20]})


def run_audit(manuscript: Path, status: str = "PRE_SUBMISSION", submission: Path | None = None, build_dir: Path | None = None, code_dir: Path | None = None) -> dict:
    return Audit(Path(manuscript), status, Path(submission) if submission else None, Path(build_dir) if build_dir else None, Path(code_dir) if code_dir else None).run()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manuscript", required=True, type=Path)
    ap.add_argument("--submission", type=Path, default=None)
    ap.add_argument("--build-dir", type=Path, default=None, help="directory holding main.log/main.aux/main.bbl of the compiled article")
    ap.add_argument("--code-dir", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--status", choices=["PRE_SUBMISSION", "FINAL"], required=True)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args(argv)
    result = run_audit(args.manuscript, args.status, args.submission, args.build_dir, args.code_dir)
    for r in result["rules"]:
        print(f"{r['status']:4} {r['id']:45} {r['description']}")
        if r["status"] == "FAIL":
            print("     " + json.dumps(r["detail"], default=str)[:600])
    c = result["counts"]
    print(f"audit: {c['PASS']} passed, {c['FAIL']} failed, {c['SKIP']} skipped (status {args.status})")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=1, default=str) + "\n", encoding="utf-8")
        print(f"wrote {args.json}")
    return 0 if c["FAIL"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
