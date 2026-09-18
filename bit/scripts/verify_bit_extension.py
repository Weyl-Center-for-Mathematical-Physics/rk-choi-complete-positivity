#!/usr/bin/env python3
"""Symbolic evidence for the extrapolation corollary (BIT revision extension E1).

Writes results/bit_revision/extrapolation_corollary.json and .md.  Exit code 0 only if every check passes.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rk_choi_margin.extrapolation import (  # noqa: E402
    frequency_factor,
    local_orientation_coefficient,
    predicted_extrapolation_defect,
    richardson_extrapolate,
)
from rk_choi_margin.methods import X, all_methods  # noqa: E402
from rk_choi_margin.symbolic import first_exponential_defect  # noqa: E402

OUT = ROOT / "results" / "bit_revision"
OUT.mkdir(parents=True, exist_ok=True)

checks: list[dict] = []
ok_all = True

# 1. Generic symbolic identity: R = truncated exponential of degree p with symbolic defects eta1, eta2.
eta1, eta2 = sp.symbols("eta1 eta2")
for p in range(1, 7):
    Rg = sum(X**j / sp.factorial(j) for j in range(p + 1)) + (1 / sp.factorial(p + 1) + eta1) * X ** (p + 1) + (1 / sp.factorial(p + 2) + eta2) * X ** (p + 2)
    for n in (2, 3, 4):
        E = richardson_extrapolate(Rg, p, n)
        ser = sp.expand(sp.series(E - sp.exp(X), X, 0, p + 3).removeO())
        low = all(sp.simplify(ser.coeff(X, k)) == 0 for k in range(p + 2))
        lead = sp.simplify(ser.coeff(X, p + 2) - predicted_extrapolation_defect(eta1, eta2, p, n))
        good = low and lead == 0
        ok_all &= good
        checks.append({"kind": "generic", "p": p, "n": n, "pass": good})

# 2. Archived methods: computed defect versus corollary, and local orientation at varpi = 0.
table = []
for key, method in all_methods().items():
    R = method.stability_function()
    p = method.order
    ser = sp.series(R - sp.exp(X), X, 0, p + 3).removeO()
    e1 = sp.simplify(ser.coeff(X, p + 1))
    e2 = sp.simplify(ser.coeff(X, p + 2))
    for n in (2, 3, 4):
        E = richardson_extrapolate(R, p, n)
        d = first_exponential_defect(E)
        pred = predicted_extrapolation_defect(e1, e2, p, n)
        if e1 == e2:
            good = d.index >= p + 3
            orient = None
        else:
            good = d.index == p + 2 and sp.simplify(d.coefficient - pred) == 0
            orient = local_orientation_coefficient(d.index, d.coefficient, 0)
        ok_all &= good
        table.append({
            "method": method.name, "key": key, "p": p, "n": n,
            "eta_p1": str(e1), "eta_p2": str(e2), "degenerate": bool(e1 == e2),
            "defect_index": d.index, "defect_coefficient": str(d.coefficient), "predicted": str(pred),
            "orientation_coefficient_varpi0": None if orient is None else str(orient),
            "locally_CP_at_varpi0": None if orient is None else bool(orient > 0),
            "pass": good,
        })

# 3. RK4 two-fold special case and frequency bands (must agree with the archive certificate).
R4 = all_methods()["rk4"].stability_function()
E4 = sp.expand(richardson_extrapolate(R4, 4, 2))
d4 = first_exponential_defect(E4)
v = sp.Symbol("varpi", real=True)
q = 1 + 4 * v**2
bands_ok = (
    (d4.index, d4.coefficient) == (6, -sp.Rational(1, 4320))
    and sp.simplify(frequency_factor(6, v) - q * (q**2 - 18 * q + 48) / 32) == 0
    and local_orientation_coefficient(6, d4.coefficient, 0) < 0
    and local_orientation_coefficient(6, d4.coefficient, 1) > 0
    and local_orientation_coefficient(6, d4.coefficient, 5) < 0
)
ok_all &= bands_ok
checks.append({"kind": "rk4_special_case", "pass": bands_ok, "eta6": str(d4.coefficient)})

# 4. Truncated-exponential closed form and parity rule.
parity_ok = True
for key, p in (("euler", 1), ("rk2", 2), ("rk3", 3), ("rk4", 4)):
    for n in (2, 3):
        closed = -sp.Rational((n - 1) * (p + 1), n * (n**p - 1)) / sp.factorial(p + 2)
        E = richardson_extrapolate(all_methods()[key].stability_function(), p, n)
        d = first_exponential_defect(E)
        coeff = local_orientation_coefficient(d.index, d.coefficient, 0)
        parity_ok &= (d.index == p + 2) and sp.simplify(d.coefficient - closed) == 0 and ((coeff < 0) == (p % 2 == 0))
ok_all &= parity_ok
checks.append({"kind": "truncated_exponential_parity", "pass": parity_ok})

def _plain(obj):
    """Convert SymPy booleans/expressions to JSON-serializable values."""
    if isinstance(obj, dict):
        return {k: _plain(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_plain(v) for v in obj]
    if obj is sp.true or obj is sp.false:
        return bool(obj)
    if isinstance(obj, sp.Basic):
        return str(obj)
    return obj


summary = _plain({"script": "scripts/verify_bit_extension.py", "sympy": sp.__version__, "all_pass": bool(ok_all), "generic_and_special_checks": checks, "method_table": table})
(OUT / "extrapolation_corollary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
lines = ["# Extrapolation corollary evidence", "", f"All checks pass: {ok_all}", "",
         "| Method | p | n | eta_{p+1} | eta_{p+2} | defect order | eta^E | predicted | orientation at varpi=0 |", "|---|---|---|---|---|---|---|---|---|"]
for r in table:
    lines.append(f"| {r['method']} | {r['p']} | {r['n']} | {r['eta_p1']} | {r['eta_p2']} | {r['defect_index']} | {r['defect_coefficient']} | {r['predicted']} | {'degenerate' if r['degenerate'] else ('CP' if r['locally_CP_at_varpi0'] else 'non-CP')} |")
(OUT / "extrapolation_corollary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))
sys.exit(0 if ok_all else 1)
