#!/usr/bin/env python3
"""Cross-platform reproduction orchestrator for the BIT Numerical Mathematics revision.

Stages (run in this order by ``all``):
  certificates        independent first-principles verification of the article's claims and of the
                      extrapolation corollary (SymPy exact arithmetic; writes results/bit_revision/*.json|md)
  tests               the complete pytest suite, one isolated process per module (historical v3.5 modules
                      plus the BIT contract modules); a collection gate records the number of tests
  audits              the retained v3.x audit scripts whose checks the article still relies on
  figures             deterministic regeneration of the BIT figure set and the Online Resource tables
  manuscript-contract the manuscript source audit (needs the revision-root ``manuscript`` directory;
                      skipped by ``all`` when the archive is used stand-alone)

Every stage prints the exact command lines it runs, stops at the first failure, returns that failure's
exit code, and appends a timestamped log under results/bit_revision/logs/.  ``all`` finishes by writing
results/bit_revision/manifest.json (stage outcomes, pass/fail/skip counts, test totals, SHA-256 of every
generated artifact).  No file is written outside this archive except the optional manuscript-audit JSON,
whose location is given explicitly on the command line.

Examples (PowerShell):   .\\.venv\\Scripts\\python.exe .\\reproduce_bit.py all
          (POSIX):        ./.venv/bin/python reproduce_bit.py tests
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "bit_revision"
DEFAULT_LOG_DIR = RESULTS / "logs"
MANIFEST = RESULTS / "manifest.json"
PY = sys.executable

TEST_MODULES = [
    # historical v3.5 modules (order of reproduce.sh) ...
    "tests/test_closeout.py",
    "tests/test_core.py",
    "tests/test_topology_operational.py",
    "tests/test_v26.py",
    "tests/test_v27_certification.py",
    "tests/test_v27_noncommuting_adaptive.py",
    "tests/test_v28_candidates.py",
    "tests/test_v29_executed_candidate.py",
    "tests/test_v29_theory.py",
    "tests/test_v30_execution_contract.py",
    "tests/test_v31_closeout.py",
    "tests/test_v32_external_review.py",
    "tests/test_v33_endpoint_membership.py",
    "tests/test_v34_remediation.py",
    "tests/test_v35_editorial_layout.py",
    "tests/test_v35_submission_revision.py",
    # ... followed by the BIT revision contract modules
    "tests/test_bit_scientific_contract.py",
    "tests/test_bit_figure_contract.py",
    "tests/test_bit_table_contract.py",
    "tests/test_bit_reproduce_contract.py",
    "tests/test_bit_submission_contract.py",
]

# Artifacts hashed into the manifest after ``all`` (relative to ROOT).
ARTIFACT_GLOBS = [
    "results/bit_revision/claim_verification.json",
    "results/bit_revision/claim_verification.md",
    "results/bit_revision/extrapolation_corollary.json",
    "results/bit_revision/extrapolation_corollary.md",
    "results/bit_revision/figures/manifest.json",
    "results/bit_revision/tables/manifest.json",
    "results/bit_revision/tables/*.tex",
    "figures_bit/*.eps",
    "figures_bit/*.pdf",
    "figures_bit/*.png",
]

# Exact release contract: a missing output must not disappear silently from a glob.
EXPECTED_ARTIFACTS = {
    f"results/bit_revision/{stem}.{ext}"
    for stem in ("claim_verification", "extrapolation_corollary") for ext in ("json", "md")
} | {"results/bit_revision/figures/manifest.json", "results/bit_revision/tables/manifest.json"} | {
    f"results/bit_revision/tables/{stem}.tex"
    for stem in ("extrapolation_table", "guard_telemetry_table", "tolerance_sweep_table")
} | {f"figures_bit/{stem}.{ext}" for stem in ("Fig1", "Fig2", "Fig3", "FigS1", "FigS2", "FigS3", "FigS4")
     for ext in ("eps", "pdf", "png")}

# Frozen records consumed by the figure and table generators (scientific inputs; hashed for provenance).
FROZEN_INPUTS = [
    "results/jcp_figures_v34/figure1_exact_branches.csv",
    "results/jcp_figures_v34/figure1_conditioning.csv",
    "results/jcp_figures_v34/figure2_candidate_regions.csv",
    "results/jcp_figures_v34/figure2_richardson_margin.csv",
    "results/jcp_figures_v34/figure3_noncommuting_endpoints.csv",
    "results/jcp_figures_v34/figure4_conditioning.csv",
    "results/jcp_figures_v34/figure5_tolerance_sweep.csv",
    "results/jcp_figures_v34/figure5_action_composition.csv",
    "results/jcp_figures_v34/figureS2_operational.csv",
    "results/v27/v27_verification.json",
    "results/v34/v34_verification.json",
    "results/closeout/root_sign_certificate_v2_7.json",
]


def _env() -> dict[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + str(ROOT / "scripts") + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    env["SOURCE_DATE_EPOCH"] = "1787184000"
    env["MPLBACKEND"] = "Agg"
    env["PYTHONIOENCODING"] = "utf-8"
    return env


def stage_commands(stage: str, manuscript: Path | None, submission: Path | None, audit_json: Path | None) -> list[tuple[str, list[str]]]:
    if stage == "certificates":
        return [
            ("Independent first-principles verification of the article's claims", [PY, "-u", "scripts/verify_bit_claims.py"]),
            ("Verification of the extrapolation corollary for the archived methods", [PY, "-u", "scripts/verify_bit_extension.py"]),
        ]
    if stage == "tests":
        cmds = [("Collect the complete pytest suite", [PY, "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider", "tests"])]
        for mod in TEST_MODULES:
            if (ROOT / mod).exists():
                cmds.append((f"Test {mod}", [PY, "-m", "pytest", "-q", "-rs", "-p", "no:cacheprovider", mod]))
        return cmds
    if stage == "audits":
        return [
            ("Verify theorem identities (v2.x archive audit)", [PY, "-u", "scripts/verify_theorems.py"]),
            ("Verify phase-covariant identities", [PY, "-u", "scripts/verify_phase_covariant.py"]),
            ("Verify method and candidate audit", [PY, "-u", "scripts/verify_section4.py"]),
            ("Independent 15,000-sample Choi audit", [PY, "-u", "independent_validation/independent_audit_v31.py"]),
            ("Executed-candidate contract audit", [PY, "-u", "independent_validation/executed_candidate_contract_v31.py"]),
            ("Endpoint-domain adversarial audit", [PY, "-u", "independent_validation/endpoint_membership_v33.py"]),
        ]
    if stage == "figures":
        return [
            ("Generate the BIT figure set with overlap and colour audits", [PY, "-u", "scripts/generate_bit_figures.py"]),
            ("Generate the Online Resource tables from frozen records", [PY, "-u", "scripts/generate_bit_tables.py"]),
        ]
    if stage == "manuscript-contract":
        if manuscript is None:
            return []
        cmd = [PY, "-u", "scripts/audit_bit_source.py", "--manuscript", str(manuscript), "--status", "PRE_SUBMISSION"]
        if submission is not None:
            cmd += ["--submission", str(submission)]
        if audit_json is not None:
            cmd += ["--json", str(audit_json)]
        return [("Audit the manuscript source against the archive contract", cmd)]
    raise ValueError(stage)


STAGES = ["certificates", "tests", "audits", "figures", "manuscript-contract"]


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _inside(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def run_stage(stage: str, log_dir: Path, manuscript: Path | None, submission: Path | None, audit_json: Path | None) -> dict:
    cmds = stage_commands(stage, manuscript, submission, audit_json)
    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f"{stamp}_{stage}.log"
    record = {"stage": stage, "log": str(log.relative_to(ROOT)) if _inside(log, ROOT) else str(log), "commands": [], "status": "PASS"}
    if not cmds:
        record["status"] = "SKIPPED"
        record["reason"] = "manuscript directory not available; run with --manuscript to audit the article source"
        log.write_text(f"[{stamp}] stage {stage}: SKIPPED ({record['reason']})\n", encoding="utf-8")
        print(f"==> {stage}: SKIPPED ({record['reason']})")
        return record
    with log.open("w", encoding="utf-8") as fh:
        fh.write(f"[{stamp}] stage {stage} in {ROOT}\n")
        for title, argv in cmds:
            line = " ".join(argv)
            print(f"==> {title}\n    $ {line}")
            fh.write(f"\n==> {title}\n$ {line}\n")
            t0 = time.time()
            proc = subprocess.run(argv, cwd=ROOT, env=_env(), capture_output=True, text=True, encoding="utf-8", errors="replace")
            dt = time.time() - t0
            fh.write(proc.stdout)
            if proc.stderr:
                fh.write("\n[stderr]\n" + proc.stderr)
            fh.write(f"\n<== exit {proc.returncode} after {dt:.1f} s\n")
            # the manifest records portable command lines: the interpreter as "python" and workspace paths relative to the archive
            shown = [("python" if a == PY else a.replace(str(ROOT.parent), "<revision-root>").replace(str(ROOT), "<archive>")) for a in argv]
            entry = {"title": title, "argv": shown, "exit": proc.returncode, "seconds": round(dt, 1)}
            m = re.search(r"(\d+) passed", proc.stdout)
            if m and "pytest" in line:
                entry["tests_passed"] = int(m.group(1))
            m = re.search(r"(\d+) skipped", proc.stdout)
            if m and "pytest" in line:
                entry["tests_skipped"] = int(m.group(1))
            m = re.search(r"(\d+) tests? collected", proc.stdout)
            if m and "--collect-only" in line:
                entry["tests_collected"] = int(m.group(1))
            record["commands"].append(entry)
            tail = proc.stdout.strip().splitlines()[-1:] or [""]
            print(f"<== exit {proc.returncode} ({dt:.1f} s) {tail[0][:120]}")
            if proc.returncode != 0:
                record["status"] = "FAIL"
                record["failed_exit"] = proc.returncode
                print(f"<!! {title} failed (exit {proc.returncode}); see {log}", file=sys.stderr)
                break
    return record


def source_hashes(root: Path | None = None) -> dict[str, str]:
    """Bind the computation, tests, dependency pins and packaging contract, not transient outputs."""
    root = ROOT if root is None else root
    files = {root / n for n in ("reproduce_bit.py", "reproduce.sh", "pyproject.toml", "CODE_ARCHIVE_ALLOWLIST.txt")
             if (root / n).is_file()}
    for directory in ("src", "scripts", "tests", "independent_validation", "manuscript_contract"):
        files.update(p for p in (root / directory).rglob("*.py") if p.is_file())
    return {p.relative_to(root).as_posix(): sha256(p) for p in sorted(files)}


def frozen_hashes() -> dict[str, str]:
    return {rel: sha256(ROOT / rel) for rel in FROZEN_INPUTS if (ROOT / rel).is_file()}


def manuscript_hashes(manuscript: Path | None) -> dict[str, str]:
    if manuscript is None or not manuscript.is_dir():
        return {}
    # article_numbers.tex is derived by the package builder and checked against main.aux there.
    return {p.name: sha256(p) for p in sorted(manuscript.iterdir())
            if p.is_file() and p.suffix in {".tex", ".bib", ".eps", ".cls", ".bst"}
            and p.name != "article_numbers.tex"}


def write_manifest(records: list[dict], path: Path = MANIFEST, source_inputs: dict | None = None,
                   frozen_inputs: dict | None = None, manuscript_inputs: dict | None = None) -> Path:
    artifacts = {}
    for pattern in ARTIFACT_GLOBS:
        for p in sorted(ROOT.glob(pattern)):
            if p.is_file():
                artifacts[p.relative_to(ROOT).as_posix()] = {"bytes": p.stat().st_size, "sha256": sha256(p)}
    frozen = frozen_hashes() if frozen_inputs is None else frozen_inputs
    tests_passed = sum(c.get("tests_passed", 0) for r in records for c in r["commands"])
    collected = next((c["tests_collected"] for r in records for c in r["commands"] if "tests_collected" in c), None)
    counts = {s: sum(1 for r in records if r["status"] == s) for s in ("PASS", "FAIL", "SKIPPED")}
    manifest = {
        "generator": "reproduce_bit.py",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "stage_counts": counts,
        "stages": records,
        "tests_collected": collected,
        "tests_passed": tests_passed,
        "tests_skipped": sum(c.get("tests_skipped", 0) for r in records for c in r["commands"]),
        "frozen_inputs_sha256": frozen,
        "source_inputs_sha256": source_hashes() if source_inputs is None else source_inputs,
        "manuscript_inputs_sha256": manuscript_hashes(ROOT.parent / "manuscript") if manuscript_inputs is None else manuscript_inputs,
        "artifacts": artifacts,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=STAGES + ["all"])
    ap.add_argument("--manuscript", type=Path, default=None, help="revision-root manuscript directory for the manuscript-contract stage")
    ap.add_argument("--submission", type=Path, default=None, help="submission directory (optional, forwarded to the source audit)")
    ap.add_argument("--audit-json", type=Path, default=None, help="where the source audit writes its JSON (default: results/bit_revision/bit_audit.json)")
    ap.add_argument("--log-dir", type=Path, default=DEFAULT_LOG_DIR, help="log directory; must lie inside this archive")
    ap.add_argument("--list", action="store_true", help="print the commands of the selected stage(s) and exit")
    args = ap.parse_args(argv)

    log_dir = args.log_dir if args.log_dir.is_absolute() else ROOT / args.log_dir
    if not _inside(log_dir, ROOT):
        print(f"error: --log-dir must be inside the archive root {ROOT}", file=sys.stderr)
        return 2
    manuscript = args.manuscript
    if manuscript is None and (ROOT.parent / "manuscript" / "main.tex").exists():
        manuscript = ROOT.parent / "manuscript"
    if manuscript is not None and not (manuscript / "main.tex").exists():
        print(f"error: {manuscript} does not contain main.tex", file=sys.stderr)
        return 2
    submission = args.submission
    audit_json = args.audit_json or (RESULTS / "bit_audit.json")

    stages = STAGES if args.stage == "all" else [args.stage]
    if args.list:
        for st in stages:
            print(f"[{st}]")
            for title, cmd in stage_commands(st, manuscript, submission, audit_json):
                print("  " + " ".join(cmd) + f"    # {title}")
        return 0
    if args.stage == "manuscript-contract" and manuscript is None:
        print("error: the manuscript-contract stage needs --manuscript <dir containing main.tex>", file=sys.stderr)
        return 2

    initial_sources = source_hashes() if args.stage == "all" else None
    initial_frozen = frozen_hashes() if args.stage == "all" else None
    initial_manuscript = manuscript_hashes(manuscript) if args.stage == "all" else None
    records = []
    rc = 0
    for st in stages:
        rec = run_stage(st, log_dir, manuscript, submission, audit_json)
        if st == "tests" and rec["status"] == "PASS":
            collected = next((c.get("tests_collected") for c in rec["commands"] if "tests_collected" in c), None)
            passed = sum(c.get("tests_passed", 0) for c in rec["commands"])
            skipped = sum(c.get("tests_skipped", 0) for c in rec["commands"])
            if not collected or collected != passed + skipped or (skipped and manuscript is not None):
                rec.update(status="FAIL", failed_exit=1, reason=f"test coverage mismatch: collected {collected}, passed {passed}, skipped {skipped}")
                print(rec["reason"], file=sys.stderr)
        records.append(rec)
        if rec["status"] == "FAIL":
            rc = rec["failed_exit"] or 1
            break
    if args.stage == "all":
        changes = [label for label, before, after in (
            ("computation sources", initial_sources, source_hashes()),
            ("frozen inputs", initial_frozen, frozen_hashes()),
            ("tested manuscript", initial_manuscript, manuscript_hashes(manuscript))) if before != after]
        if changes:
            records.append({"stage": "input-consistency", "status": "FAIL", "commands": [],
                            "reason": "inputs changed during reproduction: " + ", ".join(changes)})
            rc = 1
        path = write_manifest(records, source_inputs=initial_sources, frozen_inputs=initial_frozen,
                              manuscript_inputs=initial_manuscript)
        counts = {s: sum(1 for r in records if r["status"] == s) for s in ("PASS", "FAIL", "SKIPPED")}
        print(f"==> manifest written to {path}")
        print(f"==> stages: {counts['PASS']} passed, {counts['FAIL']} failed, {counts['SKIPPED']} skipped; "
              f"tests passed: {sum(c.get('tests_passed', 0) for r in records for c in r['commands'])}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
