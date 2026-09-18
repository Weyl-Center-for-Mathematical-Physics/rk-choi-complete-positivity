#!/usr/bin/env python3
"""Build the BIT submission candidate package from allowlisted sources (Task 16/19).

    python scripts/build_bit_package.py --revision-root .. --status PRE_SUBMISSION [--date 2026-09-13]

Steps (each compile runs in a fresh empty directory under <root>/build/package/ so that the flattened
sources are proven self-contained):
  1. flatten manuscript/main.tex (every \\input inlined) and compile the clean article; compile the review
     copy from the same flattened source with the class option `lineno`;
  2. regenerate article_numbers.tex from the clean build's main.aux and require it to equal the working copy;
  3. compile Online Resource 1 from its own flat source set; compile the cover letter; export its plain text;
  4. write BIT_Manuscript_Source.zip, BIT_ESM_1_Source.zip, BIT_Code_and_Data.zip (deterministic, flat where
     required, allowlisted);
  5. write BIT_README_FIRST.md, BIT_SHA256.txt (verified from scratch), the outer package ZIP and its
     sibling .sha256; re-open the outer ZIP, test CRC, compare every member byte-for-byte with the staged
     files, and re-verify the internal manifest;
  6. write build/audit/bit_package.json.
Nothing is uploaded or sent anywhere. A FINAL local build additionally requires --gates-resolved for local
finalization conditions (not later portal review/submission authority) and a marker-free
source audit.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
CODE_ROOT = HERE.parent
ZIP_DATE = (2026, 8, 16, 0, 0, 0)  # fixed member timestamp (SOURCE_DATE_EPOCH 1787184000)

MAIN_SOURCE_EXTRA = ["references.bib", "sn-jnl.cls", "sn-mathphys-num.bst"]
ESM_SOURCE = ["ESM_1.tex", "article_numbers.tex", "nine_candidate_certificate.tex", "extrapolation_table.tex",
              "tolerance_sweep_table.tex", "guard_telemetry_table.tex"]
CODE_EXCLUDE_DIRS = {".venv", "__pycache__", ".pytest_cache", ".matplotlib-cache", ".git", "logs"}
CODE_EXCLUDE_FILES = {"bit_audit.json"}
CODE_EXCLUDE_PREFIX = ("figures_jcp/Graphical_Abstract",)
CODE_EXCLUDE_SUFFIX = (".pyc", ".aux", ".synctex.gz", ".fdb_latexmk", ".fls")
PACKAGE_DIRS = {"main_clean", "main_review", "esm", "cover", "extract_BIT_Manuscript_Source",
                "extract_BIT_ESM_1_Source", "outer_check"}

HUMAN_GATES = [
    ("Exclusive consideration", "the manuscript is neither under consideration nor accepted for publication elsewhere"),
    ("Scientific approval", "all authors approve the final theorem statements, title, abstract, figures, Online Resource, limitations, and comparisons with prior work"),
    ("Changed metadata approval", "approve revised article fields with the manuscript; retain approved identities, order, affiliations (including institution-only ODU), contact details, ORCIDs and roles unless authors report a factual change"),
    ("Funding/acknowledgment approval", "confirm the Wahab/ODU acknowledgment and whether the 'no specific grant' statement is accurate alongside NSF-supported computing infrastructure"),
    ("AI disclosure approval", "approve the generative-AI statement in the Declarations"),
    ("Prior-release statement", "the cited Weyl Center note and historical GitHub/Zenodo release are distinguished from the present article and code release"),
    ("Data/code repository", "under the project's chosen release requirement, authorize deposit of the verified archive and supply its persistent DOI/URL; retain the historical regression fixture; BIT's mandatory Data Availability Statement is distinct from this release requirement"),
    ("Reviewer/editor metadata", "optional editor suggestions require conflict checks if used; they do not block a final local package"),
    ("Final visual approval", "an author inspects the exact PDFs of this package, including every figure at page scale"),
    ("Portal-generated PDF", "an author inspects the journal portal's compiled PDF and confirms file order, rendering, declarations, and metadata"),
    ("Submission authorization", "an author explicitly instructs the operator to submit; package preparation is not authorization"),
]


# ---------------------------------------------------------------------------
# helpers (unit-tested)
# ---------------------------------------------------------------------------
def flatten_tex(main: Path) -> str:
    def expand(p: Path) -> str:
        out = []
        for line in p.read_text(encoding="utf-8").splitlines(keepends=True):
            stripped = line.split("%", 1)[0] if not line.lstrip().startswith("%") else ""
            m = re.search(r"\\input\{([^}]+)\}", stripped)
            if m:
                name = m.group(1)
                q = p.parent / (name if name.endswith(".tex") else name + ".tex")
                out.append(f"%% ---- begin {q.name} ----\n")
                out.append(expand(q))
                if not out[-1].endswith("\n"):
                    out.append("\n")
                out.append(f"%% ---- end {q.name} ----\n")
            else:
                out.append(line)
        return "".join(out)

    return expand(main)


def review_variant(flat: str) -> str:
    new, n = re.subn(r"\\documentclass\[pdflatex,sn-mathphys-num\]\{sn-jnl\}", r"\\documentclass[pdflatex,sn-mathphys-num,lineno]{sn-jnl}", flat, count=1)
    if n != 1:
        raise ValueError("documentclass line not found")
    return new


def write_flat_zip(target: Path, members: dict[str, Path]) -> None:
    for name in members:
        if "/" in name or "\\" in name or name.startswith("..") or Path(name).is_absolute():
            raise ValueError(f"flat archive member name not allowed: {name}")
    _write_zip(target, {name: members[name] for name in sorted(members)})


def write_tree_zip(target: Path, members: dict[str, Path]) -> None:
    for name in members:
        parts = name.replace("\\", "/").split("/")
        if ".." in parts or name.startswith("/") or Path(name).is_absolute():
            raise ValueError(f"unsafe member name: {name}")
    _write_zip(target, {name.replace("\\", "/"): members[name] for name in sorted(members)})


def _write_zip(target: Path, members: dict[str, Path]) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, src in members.items():
            info = zipfile.ZipInfo(name, date_time=ZIP_DATE)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, src.read_bytes())


def verify_flat_zip(target: Path, members: dict[str, Path]) -> dict:
    with zipfile.ZipFile(target) as zf:
        crc_ok = zf.testzip() is None
        names = zf.namelist()
        flat = all("/" not in n for n in names)
        safe = all(not n.startswith("/") and ".." not in n.split("/") for n in names)
        bytes_match = set(names) == set(members) and all(zf.read(n) == members[n].read_bytes() for n in names if n in members)
    return {"crc_ok": crc_ok, "flat": flat, "safe": safe, "bytes_match": bytes_match, "members": len(names)}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write_lf(path: Path, text: str) -> None:
    """Write UTF-8 text with LF line endings on every platform (sha256sum -c and POSIX tools expect LF)."""
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


def sha256_manifest(files: list[Path], base: Path) -> str:
    return "".join(f"{sha256(p)}  {p.relative_to(base).as_posix()}\n" for p in sorted(files, key=lambda q: q.name))


def verify_sha256_manifest(text: str, base: Path) -> list[str]:
    bad = []
    for line in text.splitlines():
        if not line.strip():
            continue
        h, _, name = line.partition("  ")
        p = base / name.strip()
        if not p.exists() or sha256(p) != h.strip():
            bad.append(name.strip())
    return bad


def code_archive_files(code_root: Path) -> list[Path]:
    inventory = code_root / "CODE_ARCHIVE_ALLOWLIST.txt"
    if not inventory.is_file():
        raise ValueError("explicit CODE_ARCHIVE_ALLOWLIST.txt is required")
    out = []
    for line in inventory.read_text(encoding="utf-8").splitlines():
        name = line.strip()
        if not name or name.startswith("#"):
            continue
        rel = Path(name)
        if rel.is_absolute() or ".." in rel.parts or "\\" in name or ":" in name:
            raise ValueError(f"unsafe allowlist path: {name}")
        parts = rel.parts
        if any(part in CODE_EXCLUDE_DIRS or part.endswith(".egg-info") for part in parts[:-1]):
            raise ValueError(f"environment/cache in allowlist: {name}")
        if (rel.name in CODE_EXCLUDE_FILES or rel.name.lower().endswith(CODE_EXCLUDE_SUFFIX)
                or rel.as_posix().startswith(CODE_EXCLUDE_PREFIX) or rel.name.startswith(".")):
            raise ValueError(f"excluded file in allowlist: {name}")
        p = code_root / rel
        reject_reparse_path(p)
        if not p.is_file() or not p.resolve().is_relative_to(code_root.resolve()):
            raise ValueError(f"missing or escaping allowlist file: {name}")
        if rel in out:
            raise ValueError(f"duplicate allowlist path: {name}")
        out.append(rel)
    return sorted(out)


def validate_reproduction_manifest(code_root: Path, manifest: dict, manuscript: Path | None = None) -> None:
    """Reject partial runs, missing outputs, changed frozen data and stale proof code."""
    sys.path.insert(0, str(CODE_ROOT))
    import reproduce_bit as rb
    manuscript = code_root.parent / "manuscript" if manuscript is None else manuscript

    records = manifest.get("stages", [])
    if ([r.get("stage") for r in records] != rb.STAGES
            or manifest.get("stage_counts") != {"PASS": 5, "FAIL": 0, "SKIPPED": 0}
            or any(r.get("status") != "PASS" or not r.get("commands")
                   or any(c.get("exit") != 0 for c in r["commands"]) for r in records)):
        raise ValueError("reproduction requires all five successful stages")
    collected, passed = manifest.get("tests_collected"), manifest.get("tests_passed")
    actual_passed = sum(c.get("tests_passed", 0) for r in records for c in r["commands"])
    actual_collected = [c["tests_collected"] for r in records for c in r["commands"] if "tests_collected" in c]
    if (type(collected) is not int or collected <= 0 or collected != passed or passed != actual_passed
            or actual_collected != [collected] or manifest.get("tests_skipped", 0) != 0
            or any(c.get("tests_skipped", 0) for r in records for c in r["commands"])):
        raise ValueError("reproduction test totals are incomplete or inconsistent")
    expected = {"artifacts": rb.EXPECTED_ARTIFACTS, "frozen_inputs_sha256": set(rb.FROZEN_INPUTS),
                "source_inputs_sha256": set(rb.source_hashes(code_root)),
                "manuscript_inputs_sha256": set(rb.manuscript_hashes(manuscript))}
    for field, names in expected.items():
        entries = manifest.get(field, {})
        if not names or set(entries) != names:
            raise ValueError(f"reproduction {field} inventory is incomplete or stale")
        for rel, entry in entries.items():
            p = (manuscript if field == "manuscript_inputs_sha256" else code_root) / rel
            reject_reparse_path(p)
            expected_hash = entry.get("sha256") if field == "artifacts" else entry
            if not p.is_file() or sha256(p) != expected_hash:
                raise ValueError(f"reproduction hash mismatch: {rel}")
            if field == "artifacts" and p.stat().st_size != entry.get("bytes"):
                raise ValueError(f"reproduction size mismatch: {rel}")


def reject_reparse_path(path: Path) -> None:
    for part in (path.absolute(), *path.absolute().parents):
        if part.is_symlink() or getattr(part, "is_junction", lambda: False)():
            raise ValueError(f"symlink/junction path is not allowed: {part}")


def preserve_previous_packages(submission: Path, history: Path) -> None:
    """Keep every previous outer ZIP and sibling checksum in a unique recoverable directory."""
    reject_reparse_path(submission)
    reject_reparse_path(history)
    if history.name != "history" or history.parent.name != "build" or history.parent.parent.resolve() != submission.parent.resolve():
        raise ValueError("package history must be under the revision root's build/history")
    old = sorted(submission.glob("BIT_Submission_Package_*.zip"))
    files = old + [p.with_name(p.name + ".sha256") for p in old if p.with_name(p.name + ".sha256").exists()]
    for p in files:
        reject_reparse_path(p)
        if not p.is_file() or p.resolve().parent != submission.resolve():
            raise ValueError(f"unsafe previous package: {p}")
    if files:
        history.mkdir(parents=True, exist_ok=True)
        dest = Path(tempfile.mkdtemp(prefix="outer_packages_", dir=history))
        for p in files:
            shutil.move(str(p), str(dest / p.name))


# ---------------------------------------------------------------------------
# compile helpers
# ---------------------------------------------------------------------------
def latexmk(workdir: Path, tex: str, log_name: str) -> Path:
    env = dict(os.environ)
    env["MIKTEX_ENABLEINSTALLER"] = "1"
    env["SOURCE_DATE_EPOCH"] = "1787184000"
    cmd = ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", tex]
    print("  $", " ".join(cmd), "  (cwd", workdir, ")")
    proc = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True, errors="replace")
    (workdir / log_name).write_text(proc.stdout + "\n[stderr]\n" + proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise SystemExit(f"latexmk failed for {tex} in {workdir} (exit {proc.returncode}); see {workdir / log_name}")
    pdf = workdir / (Path(tex).stem + ".pdf")
    if not pdf.exists():
        raise SystemExit(f"no PDF produced for {tex}")
    return pdf


def fresh(d: Path) -> Path:
    reject_reparse_path(d)
    if d.name not in PACKAGE_DIRS or d.parent.name != "package" or d.parent.parent.name != "build":
        raise ValueError(f"cleanup is confined to dedicated build/package children: {d}")
    if d.resolve().parent != d.parent.resolve():
        raise ValueError(f"cleanup target escapes build/package: {d}")
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def pdf_pages(p: Path) -> int:
    from pypdf import PdfReader

    return len(PdfReader(str(p)).pages)


def pdf_text(p: Path) -> str:
    exe = shutil.which("pdftotext")
    if exe:
        proc = subprocess.run([exe, "-enc", "UTF-8", "-layout", str(p), "-"],
                              capture_output=True, text=True, encoding="utf-8", errors="strict")
        if proc.returncode == 0:
            return proc.stdout
    from pypdf import PdfReader

    return "\n".join((pg.extract_text() or "") for pg in PdfReader(str(p)).pages)


# ---------------------------------------------------------------------------
# main build
# ---------------------------------------------------------------------------
def build(root: Path, status: str, date: str, gates_resolved: bool) -> dict:
    manuscript = root / "manuscript"
    submission = root / "submission"
    build_dir = root / "build" / "package"
    audit_dir = root / "build" / "audit"
    code_root = root / "code_and_data"
    if status == "FINAL" and not gates_resolved:
        raise SystemExit("FINAL requires --gates-resolved (local finalization conditions 1-7 and 9 closed by the authors)")
    submission.mkdir(exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    report: dict = {"status": status, "date": date, "steps": []}

    # 0. source audit before building
    sys.path.insert(0, str(HERE))
    import audit_bit_source as abs_

    pre = abs_.run_audit(manuscript, status=status, code_dir=code_root)
    if pre["counts"]["FAIL"]:
        failed = [r["id"] for r in pre["rules"] if r["status"] == "FAIL"]
        raise SystemExit(f"source audit failed before packaging: {failed}")
    report["steps"].append({"step": "pre-build source audit", "counts": pre["counts"]})

    # 0b. the reproduction manifest must describe the artifacts about to be packaged (review finding B1)
    manifest_path = code_root / "results" / "bit_revision" / "manifest.json"
    if not manifest_path.exists():
        raise SystemExit("results/bit_revision/manifest.json is missing: run `reproduce_bit.py all` before packaging")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    try:
        validate_reproduction_manifest(code_root, manifest)
    except ValueError as exc:
        raise SystemExit(f"{exc}: run `reproduce_bit.py all` before packaging") from exc
    report["steps"].append({"step": "reproduction manifest verified against the files to be packaged", "artifacts": len(manifest["artifacts"]), "tests_passed": manifest.get("tests_passed")})

    # 1. flattened article: clean and review copies
    flat = flatten_tex(manuscript / "main.tex")
    # the flattened article is also shipped inside the code archive so that the article-table contract tests can run
    # without the revision workspace (review finding B9)
    ms_dir = code_root / "manuscript_source"
    ms_dir.mkdir(exist_ok=True)
    write_lf(ms_dir / "main.tex", flat)
    write_lf(ms_dir / "README.txt", "Flattened LaTeX source of the article as packaged in BIT_Manuscript_Source.zip (generated by scripts/build_bit_package.py). "
             "It is read by tests/test_bit_table_contract.py to check the typeset Tables 1 and 2 against the frozen records; it is not a build input of the archive.\n")
    figs = sorted(p.name for p in manuscript.glob("Fig[0-9].eps")) + sorted(p.name for p in manuscript.glob("Fig[0-9][0-9].eps"))
    clean_dir = fresh(build_dir / "main_clean")
    (clean_dir / "main.tex").write_text(flat, encoding="utf-8")
    for name in MAIN_SOURCE_EXTRA + figs:
        shutil.copy2(manuscript / name, clean_dir / name)
    clean_pdf = latexmk(clean_dir, "main.tex", "latexmk_main_clean.txt")
    review_dir = fresh(build_dir / "main_review")
    (review_dir / "main.tex").write_text(review_variant(flat), encoding="utf-8")
    for name in MAIN_SOURCE_EXTRA + figs:
        shutil.copy2(manuscript / name, review_dir / name)
    review_pdf = latexmk(review_dir, "main.tex", "latexmk_main_review.txt")
    report["steps"].append({"step": "article compiled from flattened source in fresh directories", "clean_pages": pdf_pages(clean_pdf), "review_pages": pdf_pages(review_pdf), "figures": figs})

    # 2. article numbers for the Online Resource
    import sync_article_numbers as san

    expected_numbers = san.render((clean_dir / "main.aux").read_text(encoding="utf-8"))
    current = (manuscript / "article_numbers.tex").read_text(encoding="utf-8") if (manuscript / "article_numbers.tex").exists() else ""
    if current != expected_numbers:
        write_lf(manuscript / "article_numbers.tex", expected_numbers)
        report["steps"].append({"step": "article_numbers.tex regenerated (was stale)"})
    # post-build audit with logs
    post = abs_.run_audit(manuscript, status=status, build_dir=clean_dir, code_dir=code_root)
    if post["counts"]["FAIL"]:
        failed = [r["id"] for r in post["rules"] if r["status"] == "FAIL"]
        raise SystemExit(f"post-build source audit failed: {failed}")

    # 3. Online Resource and cover letter
    esm_dir = fresh(build_dir / "esm")
    esm_figs = sorted(p.name for p in manuscript.glob("FigS[0-9].eps"))
    for name in ESM_SOURCE + esm_figs:
        shutil.copy2(manuscript / name, esm_dir / name)
    esm_pdf = latexmk(esm_dir, "ESM_1.tex", "latexmk_esm.txt")
    cover_dir = fresh(build_dir / "cover")
    shutil.copy2(manuscript / "cover_letter.tex", cover_dir / "cover_letter.tex")
    cover_pdf = latexmk(cover_dir, "cover_letter.tex", "latexmk_cover.txt")
    cover_txt = pdf_text(cover_pdf).replace("­", "–")  # pdftotext maps the en dash glyph to a soft hyphen
    report["steps"].append({"step": "Online Resource and cover letter compiled", "esm_pages": pdf_pages(esm_pdf), "cover_pages": pdf_pages(cover_pdf)})

    # 4. stage payload
    preserve_previous_packages(submission, root / "build" / "history")
    staged: dict[str, Path] = {}
    shutil.copy2(review_pdf, submission / "Manuscript_BIT_review.pdf")
    shutil.copy2(clean_pdf, submission / "Manuscript_BIT_clean.pdf")
    shutil.copy2(esm_pdf, submission / "ESM_1.pdf")
    shutil.copy2(cover_pdf, submission / "Cover_Letter_BIT.pdf")
    write_lf(submission / "Cover_Letter_BIT.txt", cover_txt)
    main_members = {"main.tex": clean_dir / "main.tex", "main.bbl": clean_dir / "main.bbl"}
    for name in MAIN_SOURCE_EXTRA + figs:
        main_members[name] = manuscript / name
    write_flat_zip(submission / "BIT_Manuscript_Source.zip", main_members)
    esm_members = {name: manuscript / name for name in ESM_SOURCE + esm_figs}
    write_flat_zip(submission / "BIT_ESM_1_Source.zip", esm_members)
    code_files = code_archive_files(code_root)
    write_tree_zip(submission / "BIT_Code_and_Data.zip", {f"BIT_Code_and_Data/{rel.as_posix()}": code_root / rel for rel in code_files})
    zip_checks = {
        "BIT_Manuscript_Source.zip": verify_flat_zip(submission / "BIT_Manuscript_Source.zip", main_members),
        "BIT_ESM_1_Source.zip": verify_flat_zip(submission / "BIT_ESM_1_Source.zip", esm_members),
    }
    with zipfile.ZipFile(submission / "BIT_Code_and_Data.zip") as zf:
        zip_checks["BIT_Code_and_Data.zip"] = {"crc_ok": zf.testzip() is None, "members": len(zf.namelist())}
    for name, chk in zip_checks.items():
        if not chk["crc_ok"] or chk.get("bytes_match") is False or chk.get("flat") is False:
            raise SystemExit(f"archive verification failed for {name}: {chk}")
    report["steps"].append({"step": "archives written and verified", "checks": zip_checks, "code_files": len(code_files)})

    # 4b. clean-extraction compile tests of both source ZIPs
    for zname, tex in (("BIT_Manuscript_Source.zip", "main.tex"), ("BIT_ESM_1_Source.zip", "ESM_1.tex")):
        d = fresh(build_dir / ("extract_" + zname.replace(".zip", "")))
        with zipfile.ZipFile(submission / zname) as zf:
            zf.extractall(d)
        nested = [p for p in d.rglob("*") if p.is_file() and p.parent != d]
        if nested:
            raise SystemExit(f"nested files in {zname}: {nested}")
        pdf = latexmk(d, tex, "latexmk_extracted.txt")
        report["steps"].append({"step": f"{zname} compiled independently after clean extraction", "pages": pdf_pages(pdf)})

    for name in ("BIT_Submission_Metadata.md", "BIT_Declarations.md"):
        if not (submission / name).exists():
            raise SystemExit(f"missing authored file submission/{name}")

    # 5. README, manifest, outer package
    readme = render_readme(status, date, report, code_files, gates_resolved)
    write_lf(submission / "BIT_README_FIRST.md", readme)
    payload = [submission / n for n in abs_.SUBMISSION_FILES if n != "BIT_SHA256.txt"]
    manifest_text = sha256_manifest(payload, submission)
    write_lf(submission / "BIT_SHA256.txt", manifest_text)
    bad = verify_sha256_manifest((submission / "BIT_SHA256.txt").read_text(encoding="utf-8"), submission)
    if bad:
        raise SystemExit(f"checksum verification failed: {bad}")
    outer = submission / f"BIT_Submission_Package_{date}_{status}.zip"
    members = {p.name: p for p in payload + [submission / "BIT_SHA256.txt"]}
    write_flat_zip(outer, members)
    outer_hash = sha256(outer)
    write_lf(submission / (outer.name + ".sha256"), f"{outer_hash}  {outer.name}\n")
    chk = verify_flat_zip(outer, members)
    with zipfile.ZipFile(outer) as zf:
        inner_manifest = zf.read("BIT_SHA256.txt").decode("utf-8")
        tmp = fresh(build_dir / "outer_check")
        zf.extractall(tmp)
    inner_bad = verify_sha256_manifest(inner_manifest, tmp)
    if not (chk["crc_ok"] and chk["bytes_match"] and chk["flat"] and not inner_bad):
        raise SystemExit(f"outer package verification failed: {chk} {inner_bad}")
    report["steps"].append({"step": "outer package written and verified", "outer": outer.name, "sha256": outer_hash, "members": chk["members"]})

    # 6. final audit of the staged package
    final = abs_.run_audit(manuscript, status=status, submission=submission, build_dir=clean_dir, code_dir=code_root)
    report["final_audit_counts"] = final["counts"]
    report["final_audit_failures"] = [r["id"] for r in final["rules"] if r["status"] == "FAIL"]
    report["payload_sha256"] = {p.name: sha256(p) for p in payload + [submission / "BIT_SHA256.txt"]}
    (audit_dir / "bit_package.json").write_text(json.dumps(report, indent=1, default=str) + "\n", encoding="utf-8")
    (audit_dir / "bit_audit.json").write_text(json.dumps(final, indent=1, default=str) + "\n", encoding="utf-8")
    if report["final_audit_failures"]:
        raise SystemExit(f"final package audit failed: {report['final_audit_failures']}")
    return report


def render_readme(status: str, date: str, report: dict, code_files, gates_resolved: bool) -> str:
    pages = {s.get("step"): s for s in report["steps"]}
    art = pages.get("article compiled from flattened source in fresh directories", {})
    esm = pages.get("Online Resource and cover letter compiled", {})
    def gate_status(i):
        return {8: "OPTIONAL", 10: "LATER", 11: "NOT AUTHORIZED"}.get(i, "RESOLVED" if gates_resolved else "OPEN")
    gates = "\n".join(f"| {i} | {name} | {gate_status(i)} | {action} |" for i, (name, action) in enumerate(HUMAN_GATES, 1))
    return f"""# READ ME FIRST — BIT Numerical Mathematics submission candidate

**Package status: `{status}`.** Built {date} by `code_and_data/scripts/build_bit_package.py` from the revision
workspace. {'One or more local finalization conditions remain open; this package must not be uploaded.' if status == 'PRE_SUBMISSION' else 'The authors declared the local finalization conditions (1-7 and 9) resolved; later portal inspection and submission authorization are not implied.'}
**Completing this local package is not a journal submission.** No journal files have been uploaded, and no portal
declaration has been made. Submission requires an author's explicit authorization after the gates below.

## Files

| File | What it is |
|---|---|
| `Manuscript_BIT_review.pdf` | article compiled with line numbers (class option `lineno`), {art.get('review_pages', '?')} pages |
| `Manuscript_BIT_clean.pdf` | the same source without line numbers, {art.get('clean_pages', '?')} pages |
| `ESM_1.pdf` | Online Resource 1 (Electronic Supplementary Material), {esm.get('esm_pages', '?')} pages |
| `Cover_Letter_BIT.pdf`, `Cover_Letter_BIT.txt` | cover letter (one page; body at most 150 words) and its plain text |
| `BIT_Manuscript_Source.zip` | flat LaTeX source of the article: one `main.tex`, `main.bbl`, `references.bib`, `sn-jnl.cls`, `sn-mathphys-num.bst`, {', '.join(art.get('figures', []))} |
| `BIT_ESM_1_Source.zip` | flat LaTeX source of Online Resource 1 (separate entry point `ESM_1.tex`); kept for the archive — upload `ESM_1.pdf` as the Online Resource and do **not** upload this ZIP as manuscript LaTeX source |
| `BIT_Code_and_Data.zip` | reproducibility archive ({len(code_files)} files, directory `BIT_Code_and_Data/`); its `manuscript_contract/` folder is the historical v3.5 source-contract fixture (earlier-journal files) required by the historical test suite and `reproduce.sh source-audit`, not part of the BIT submission |
| `BIT_Submission_Metadata.md` | portal field values with approval status per field |
| `BIT_Declarations.md` | declarations as typeset, with approval status |
| `BIT_SHA256.txt` | SHA-256 of every file above |
| `BIT_Submission_Package_{date}_{status}.zip` (+ `.sha256`) | the outer package containing all of the above |

## Reproduction

```powershell
Expand-Archive BIT_Code_and_Data.zip; Set-Location BIT_Code_and_Data\\BIT_Code_and_Data
python -m venv .venv; .\\.venv\\Scripts\\python.exe -m pip install -e ".[test]"
.\\.venv\\Scripts\\python.exe .\\reproduce_bit.py all
```

(POSIX: `python3 -m venv .venv && ./.venv/bin/python -m pip install -e ".[test]" && ./.venv/bin/python reproduce_bit.py all`.)

## Data-link status

{'This staging build does not attest that repository publication is complete; verify the versioned URL in the Code availability statement before finalizing.' if status == 'PRE_SUBMISSION' else 'The Code availability statement cites the verified public versioned release; the earlier Zenodo DOI remains a separate historical reference.'}

## Human gates

Only conditions 1-7 and 9 govern local FINAL preparation. Optional editor choices, later portal-proof
review, and explicit submission authorization remain separate. Unchanged approved author fields are retained.

| # | Condition | Status | Requirement |
|---|---|---|---|
{gates}

A `{status}` package name is truthful only while the statuses above are current; rebuild the package after
any change to the sources or the gate statuses.
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--revision-root", type=Path, default=CODE_ROOT.parent)
    ap.add_argument("--status", choices=["PRE_SUBMISSION", "FINAL"], required=True)
    ap.add_argument("--date", default=_dt.date.today().isoformat())
    ap.add_argument("--gates-resolved", action="store_true", help="the authors have closed local finalization conditions 1-7 and 9 (FINAL only; not submission authority)")
    args = ap.parse_args(argv)
    root = args.revision_root.resolve()
    report = build(root, args.status, args.date, args.gates_resolved)
    print(json.dumps({k: v for k, v in report.items() if k != "payload_sha256"}, indent=1, default=str))
    print(f"package built: {root / 'submission'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
