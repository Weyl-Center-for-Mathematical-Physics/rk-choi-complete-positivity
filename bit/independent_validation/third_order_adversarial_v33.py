#!/usr/bin/env python3
"""Current-version third-order adversarial release gate for JCP v3.3.

This gate verifies that the v3.1 hostile-audit defects remain closed in the
current v3.3 source/results: theorem-domain validation, current-result-bound
benchmark prose, honest fallback labeling, and strict public API validation.
It exits nonzero on any finding.
"""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "manuscript_contract"
OUT = ROOT / "independent_validation" / "third_order_adversarial_v33_report.json"
sys.path.insert(0, str(ROOT / "src"))

from rk_choi_margin import (  # noqa: E402
    build_exact_gksl_execution,
    build_strang_execution,
    ExactGKSLExponentialRecipe,
    ExactSubflowStrangRecipe,
)
from rk_choi_margin.candidates import equal_substep_candidate  # noqa: E402
from rk_choi_margin.controllers import run_adaptive_channel_benchmark  # noqa: E402

findings: list[dict] = []
raw: dict = {}

# 1. Structural theorem domain must reject negative time and accept controls.
neg = []
for h in (-1e-12, -1e-13, -1e-16, math.nextafter(0.0, -math.inf)):
    for name, fn, extra in (
        ("exact-gksl", build_exact_gksl_execution, {}),
        ("strang", build_strang_execution, {"substeps": 1}),
    ):
        try:
            fn(total_x=h, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0, dps=80, **extra)
        except (TypeError, ValueError) as exc:
            neg.append({"builder": name, "h": h, "rejected": True, "exception": type(exc).__name__})
        else:
            neg.append({"builder": name, "h": h, "rejected": False})
            findings.append({"id": "NEGATIVE-TIME-STRUCTURAL-ACCEPT", "builder": name, "h": h})
raw["negative_time"] = neg
for h in (0.0, 1e-16, 0.1):
    if build_exact_gksl_execution(total_x=h, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0).cp_status.value not in {"PASS", "STRUCTURAL_PASS"}:
        findings.append({"id": "POSITIVE-CONTROL-FAILED", "h": h})

# 2. Current benchmark and TeX macros must be v3.3-bound and labels honest.
csv_path = ROOT / "results" / "v33" / "adaptive_benchmark_v33.csv"
macro_path = SOURCE / "generated" / "benchmark_macros_v33.tex"
section_path = SOURCE / "sections" / "07_step_control.tex"
article_path = SOURCE / "jcp_article.tex"
if not (csv_path.exists() and macro_path.exists() and section_path.exists() and article_path.exists()):
    findings.append({"id": "MISSING-V33-BENCHMARK-BINDING"})
else:
    rows = list(csv.DictReader(csv_path.open(newline="", encoding="utf-8")))
    tex = section_path.read_text(encoding="utf-8")
    article = article_path.read_text(encoding="utf-8")
    macros = macro_path.read_text(encoding="utf-8")
    boundary_guard = [r for r in rows if r["scenario"] == "boundary_tolerance_sweep" and r["policy"] == "guard_rotating_fallback"]
    if not boundary_guard:
        findings.append({"id": "MISSING-BOUNDARY-GUARD-ROWS"})
    else:
        accepted = sum(int(r["accepted_steps"]) for r in boundary_guard)
        rotating = sum(int(r["rotating_accepts"]) for r in boundary_guard)
        direct = sum(int(r["direct_accepts"]) for r in boundary_guard)
        projected = sum(int(r["projected_accepts"]) for r in boundary_guard)
        if rotating != accepted or direct != 0 or projected != 0:
            findings.append({"id": "BOUNDARY-ACTION-COMPOSITION-UNEXPECTED", "accepted": accepted, "rotating": rotating, "direct": direct, "projected": projected})
        if not any(label in tex for label in ("CP guard + rotating fallback", "CP guard with rotating-frame fallback", "guard-plus-fallback")):
            findings.append({"id": "FALLBACK-POLICY-MISLABELED"})
    if "benchmark_macros_v33" not in article or "BenchmarkRotatingError" not in macros:
        findings.append({"id": "V33-MACROS-NOT-BOUND"})
    if re.search(r"3\.7\\times10\^\{-5\}.*96 stage units", tex, flags=re.S):
        findings.append({"id": "STALE-V28-BENCHMARK-LITERAL"})
    useful = {r["scenario"]: r for r in rows if r["scenario"] in {"buffered_projection", "buffered_direct", "near_saddle_projection"}}
    if int(useful.get("buffered_projection", {}).get("projected_accepts", 0)) < 1:
        findings.append({"id": "NO-USEFUL-PROJECTION"})
    if int(useful.get("buffered_direct", {}).get("direct_accepts", 0)) < 1:
        findings.append({"id": "NO-DIRECT-RETENTION"})
    if int(useful.get("near_saddle_projection", {}).get("component_searches", 0)) < 1:
        findings.append({"id": "NO-NEAR-SADDLE-SEARCH"})
    raw["benchmark_rows"] = len(rows)

# 3. Public exact-integer and controller-domain validation must fail early.
api = []
for value in (True, 2.0, 2.9, 0, -1):
    try:
        equal_substep_candidate(total_x=1, method="rk4", substeps=value, theta=0, kappa=0, varpi=0)
    except (TypeError, ValueError) as exc:
        api.append({"api": "equal_substep_candidate", "input": repr(value), "rejected": True, "exception": type(exc).__name__})
    else:
        api.append({"api": "equal_substep_candidate", "input": repr(value), "rejected": False})
        findings.append({"id": "SUBSTEPS-INVALID-ACCEPTED", "value": repr(value)})
for value in (True, 1.5, 0, -5):
    for cls in (ExactGKSLExponentialRecipe, ExactSubflowStrangRecipe):
        try:
            cls.from_values(dps=value)
        except (TypeError, ValueError) as exc:
            api.append({"api": cls.__name__, "dps": repr(value), "rejected": True, "exception": type(exc).__name__})
        else:
            api.append({"api": cls.__name__, "dps": repr(value), "rejected": False})
            findings.append({"id": "DPS-INVALID-ACCEPTED", "class": cls.__name__, "value": repr(value)})
for case in (
    {"final_time": -1.0, "initial_h": 0.1, "tolerance": 1e-2, "min_h": 1e-12},
    {"final_time": 1.0, "initial_h": -0.1, "tolerance": 1e-2, "min_h": 1e-12},
    {"final_time": 1.0, "initial_h": 0.1, "tolerance": -1e-2, "min_h": 1e-12},
    {"final_time": 1.0, "initial_h": 0.1, "tolerance": 1e-2, "min_h": -1e-12},
):
    try:
        run_adaptive_channel_benchmark(theta=0, kappa=0, varpi=2, policy="candidate_guard", fallback="exact", max_steps=3, **case)
    except (TypeError, ValueError) as exc:
        api.append({"api": "controller", "input": case, "rejected": True, "exception": type(exc).__name__})
    else:
        api.append({"api": "controller", "input": case, "rejected": False})
        findings.append({"id": "CONTROLLER-INVALID-ACCEPTED", "case": case})
raw["api_validation"] = api

report = {
    "script": Path(__file__).name,
    "code_root": str(ROOT),
    "source_root": str(SOURCE),
    "finding_count": len(findings),
    "findings": findings,
    "raw_evidence": raw,
    "overall": "PASS" if not findings else "FAIL",
}
OUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps({"overall": report["overall"], "finding_count": len(findings), "json_out": str(OUT)}, indent=2))
raise SystemExit(0 if not findings else 1)
