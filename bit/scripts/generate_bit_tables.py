#!/usr/bin/env python3
"""Generate the data tables of Online Resource 1 (BIT revision) from frozen/verified records.

Outputs (LaTeX fragments, deterministic):
  results/bit_revision/tables/extrapolation_table.tex   from results/bit_revision/extrapolation_corollary.json
  results/bit_revision/tables/tolerance_sweep_table.tex from results/jcp_figures_v34/figure5_tolerance_sweep.csv
  results/bit_revision/tables/guard_telemetry_table.tex from results/v34/v34_verification.json
  results/bit_revision/tables/manifest.json             hashes of inputs and outputs
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "bit_revision" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

INPUTS = {
    "corollary": ROOT / "results" / "bit_revision" / "extrapolation_corollary.json",
    "sweep": ROOT / "results" / "jcp_figures_v34" / "figure5_tolerance_sweep.csv",
    "v34": ROOT / "results" / "v34" / "v34_verification.json",
}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def sci(v: float, digits: int = 3) -> str:
    if v == 0:
        return "0"
    m, e = f"{v:.{digits - 1}e}".split("e")
    return rf"${m.rstrip('0').rstrip('.')}\times10^{{{int(e)}}}$"


def tex_frac(s: str) -> str:
    s = s.strip()
    if "/" in s:
        n, d = s.split("/")
        sign = "-" if n.startswith("-") else ""
        return rf"${sign}\tfrac{{{n.lstrip('-')}}}{{{d}}}$"
    return f"${s}$"


def extrapolation_table() -> str:
    data = json.loads(INPUTS["corollary"].read_text(encoding="utf-8"))
    rows = data["method_table"]
    by_method: dict[str, dict] = {}
    for r in rows:
        by_method.setdefault(r["method"], {"p": r["p"], "eta_p1": r["eta_p1"], "eta_p2": r["eta_p2"], "degenerate": r["degenerate"], "n": {}})
        by_method[r["method"]]["n"][r["n"]] = r
    lines = [r"\begin{tabular}{@{}lcccccccc@{}}", r"\toprule",
             r"Method & $p$ & $\eta_{p+1}$ & $\eta_{p+2}$ & $\eta^E_{p+2}$, $n=2$ & $n=3$ & $n=4$ & defect order & at $\varpi=0$ \\", r"\midrule"]
    display = {"Dormand-Prince 5 principal formula": "DP5(4), order-5 formula", "Dormand-Prince embedded 4 formula": "DP5(4), embedded order-4"}
    for name, m in by_method.items():
        name = display.get(name, name)
        cells = []
        for n in (2, 3, 4):
            r = m["n"][n]
            cells.append("0" if r["degenerate"] else tex_frac(r["defect_coefficient"]))
        r2 = m["n"][2]
        order = f"{r2['defect_index']}" + (" (degenerate)" if m["degenerate"] else "")
        orient = "degenerate" if m["degenerate"] else ("CPTP" if r2["locally_CP_at_varpi0"] else "non-CPTP")
        lines.append(f"{name} & {m['p']} & {tex_frac(m['eta_p1'])} & {tex_frac(m['eta_p2'])} & {cells[0]} & {cells[1]} & {cells[2]} & {order} & {orient} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def tolerance_sweep_table() -> str:
    with INPUTS["sweep"].open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    names = {"error_only": "error only", "guard_rotating_fallback": "guard, rotating fallback", "rotating_frame": "rotating frame", "strang": "Strang (exact subflows)"}
    lines = [r"\begin{tabular}{@{}lcccccc@{}}", r"\toprule",
             r"Policy & Tolerance & Choi error & Min.\ eigenvalue & Accepted & Fallback & RK stages / exp.\ actions \\", r"\midrule"]
    for r in sorted(rows, key=lambda r: (list(names).index(r["policy"]), -float(r["tolerance"]))):
        eig = float(r["minimum_step_choi_eigenvalue"])
        eig_s = r"$\ge-10^{-15}$" if eig > -5e-16 else sci(eig, 2)
        lines.append(f"{names[r['policy']]} & {sci(float(r['tolerance']), 2)} & {sci(float(r['global_normalized_choi_error']))} & {eig_s} & "
                     f"{int(float(r['accepted_steps']))} & {int(float(r['fallback_events']))} & {int(float(r['rhs_stage_evaluations']))} / {int(float(r['exponential_actions']))} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def guard_telemetry_table() -> str:
    V = json.loads(INPUTS["v34"].read_text(encoding="utf-8"))
    rows = [r for r in V["benchmark_rows"] if r["controller_policy"] == "candidate_guard" and r["scenario"] != "projection_adversarial_path_test"]
    case = {"boundary_tolerance_sweep": "boundary", "buffered_projection": "buffered (1)", "buffered_direct": "buffered (2)", "near_saddle_projection": "near saddle"}
    lines = [r"\begin{tabular}{@{}lcccccccc@{}}", r"\toprule",
             r"Case & Tol. & Locator calls & Proj. & Empty & Unc. & Isolations & Certif. & Recert. \\", r"\midrule"]
    for r in rows:
        lines.append(f"{case[r['scenario']]} & {sci(float(r['tolerance']), 2)} & {int(r['locator_calls'])} & {int(r['successful_projections'])} & {int(r['empty_component_outcomes'])} & "
                     f"{int(r['uncertain_outcomes'])} & {int(r['root_isolation_ops'])} & {int(r['certification_ops'])} & {int(r['binary_recertifications'])} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    outputs = {"extrapolation_table.tex": extrapolation_table(), "tolerance_sweep_table.tex": tolerance_sweep_table(), "guard_telemetry_table.tex": guard_telemetry_table()}
    header = "% Generated by scripts/generate_bit_tables.py from frozen records; do not edit by hand.\n"
    for name, body in outputs.items():
        (OUT / name).write_text(header + body, encoding="utf-8")
    manifest = {"inputs": {k: sha256(v) for k, v in INPUTS.items()}, "outputs": {n: sha256(OUT / n) for n in outputs}}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("wrote", ", ".join(outputs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
