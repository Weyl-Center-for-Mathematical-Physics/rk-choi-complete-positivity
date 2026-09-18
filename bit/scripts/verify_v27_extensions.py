#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from rk_choi_margin.certification import CPStatus, certify_cptp_step, locate_cptp_components
from rk_choi_margin.controllers import (
    halving_threshold_certificate,
    reduction_skip_intervals,
    rk4_detached_interval_exact,
    run_adaptive_channel_benchmark,
)
from rk_choi_margin.liouvillian import (
    certify_rk4_gksl_step_exact,
    hermitian_eigenvalue_brackets_exact,
    locate_rk4_gksl_cptp_components_exact,
)
from rk_choi_margin.topology import rk4_general_margin

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v27"
OUT.mkdir(parents=True, exist_ok=True)


def fstr(x: sp.Expr, n: int = 18) -> str:
    return str(sp.N(x, n))


def exact_margin_varpi2(x: sp.Rational) -> sp.Rational:
    q = sp.Integer(17)
    y = q * x
    cubic = y**3 - 16*y**2 - 32*(q-6)*y + 384*(q-4)
    return sp.Rational(sp.factor(-x**5*q*cubic/sp.Integer(147456)))


def direct_double_margin(x: float) -> float:
    z = complex(-x / 2, 2 * x)
    r = lambda w: 1 + w + w**2 / 2 + w**3 / 6 + w**4 / 24
    return float(r(-x).real - abs(r(z)) ** 2)


def main() -> None:
    payload: dict[str, object] = {"schema_version": "2.7"}

    bug_rows = []
    for text in ("0.004", "0.001", "0.0001", "0.00001"):
        x = sp.Rational(text)
        exact = exact_margin_varpi2(x)
        legacy_accept = direct_double_margin(float(x)) >= -1e-12
        status = certify_cptp_step("rk4", x=x, theta=0, kappa=0, varpi=2).status
        bug_rows.append(
            {
                "x": text,
                "exact_margin": str(exact),
                "exact_margin_decimal": fstr(exact, 30),
                "binary64_direct_margin": repr(direct_double_margin(float(x))),
                "legacy_tol_accept": legacy_accept,
                "certified_status": status.value,
            }
        )
    payload["blocking_false_acceptance_reproduction"] = bug_rows

    # Conditioning table: exact defect, direct binary64 subtraction, and scaled
    # sign polynomial.  The scaled sign remains order one when M~x^5.
    conditioning = []
    q = sp.Integer(17)
    y = sp.Symbol("y")
    C = y**3 - 16 * y**2 - 32 * (q - 6) * y + 384 * (q - 4)
    for exp in range(1, 13):
        x = sp.Rational(1, 10**exp)
        exact = exact_margin_varpi2(x)
        scaled = -sp.factor(C.subs(y, q * x))
        conditioning.append(
            {
                "x": str(x),
                "exact_margin": fstr(exact, 28),
                "binary64_margin": repr(direct_double_margin(float(x))),
                "scaled_sign_polynomial_minus_C": fstr(scaled, 20),
            }
        )
    payload["floating_point_conditioning"] = conditioning

    # Companion-matrix failure near saddle-node.
    qplus = (sp.Integer(123) + 11 * sp.sqrt(33)) / 16
    vp = float(sp.sqrt((qplus - 1) / 4))
    root_conditioning = []
    for ulps in range(1, 25):
        v = vp
        for _ in range(ulps):
            v = math.nextafter(v, math.inf)
        qq = 1 + 4 * v**2
        roots = np.roots([1.0, -16.0, -32.0 * (qq - 6.0), 384.0 * (qq - 4.0)])
        eligible = sorted(float(r.real / qq) for r in roots if abs(float(r.imag)) < 1e-9 and float(r.real) > 1e-12)
        root_conditioning.append({"ulps_above": ulps, "varpi": repr(v), "eligible_count": len(eligible), "eligible": eligible})
    payload["companion_matrix_saddle_node"] = root_conditioning

    half = halving_threshold_certificate(50)
    payload["halving_threshold"] = half

    # Verify resultant identity behind the threshold polynomial.
    qq, yy = sp.symbols("q y")
    cubic = yy**3 - 16 * yy**2 - 32 * (qq - 6) * yy + 384 * (qq - 4)
    resultant = sp.factor(sp.resultant(cubic, cubic.subs(yy, 2 * yy), yy))
    payload["halving_resultant"] = str(resultant)

    # Reduction theorem examples at varpi=2.
    br = rk4_detached_interval_exact(2, digits=40)
    assert br is not None
    xm = (br[0] + br[1]) / 2
    xp = (br[2] + br[3]) / 2
    payload["varpi2_detached_interval"] = {
        "x_minus_bracket": [str(br[0]), str(br[1])],
        "x_plus_bracket": [str(br[2]), str(br[3])],
        "ratio_mid": fstr(xm / xp, 25),
        "halving_skip_intervals_first_four": reduction_skip_intervals(float(xm), float(xp), 0.5, 4),
    }

    # Buffered attached-endpoint checks and exact components.
    buffered_cases = []
    for theta, kappa in ((0, "0.001"), ("0.001", 0), ("0.001", "0.001"), (0, "0.01"), ("0.01", 0)):
        loc = locate_cptp_components("rk4", theta=theta, kappa=kappa, varpi=2, upper=1.5, root_digits=30)
        comps = [(float(c.lower_value), float(c.upper_value)) for c in loc.components]
        buffered_cases.append({"theta": str(theta), "kappa": str(kappa), "components": comps})
    payload["buffered_components"] = buffered_cases

    # Noncommuting transverse-field persistence with exact component isolation
    # and exact eigenvalue brackets at four strict samples.
    noncommuting = []
    for omega_x in ("0", "0.05", "0.1", "0.2"):
        loc = locate_rk4_gksl_cptp_components_exact(
            upper="1.5", theta="0.001", kappa="0.001", omega_z="2", omega_x=omega_x, root_digits=22
        )
        comps = [(fstr(lo, 16), fstr(hi, 16)) for lo, hi, _sample in loc["components"]]
        sample_rows = []
        for h in ("0.1", "0.5", "0.9", "1.4"):
            cert = certify_rk4_gksl_step_exact(h=h, theta="0.001", kappa="0.001", omega_z="2", omega_x=omega_x)
            eig_brackets = hermitian_eigenvalue_brackets_exact(cert.choi, digits=24)
            sample_rows.append(
                {
                    "h": h,
                    "status": cert.status.value,
                    "lambda_min_bracket": [str(eig_brackets[0][0]), str(eig_brackets[0][1])],
                    "lambda_min_mid": fstr((eig_brackets[0][0] + eig_brackets[0][1]) / 2, 18),
                }
            )
        noncommuting.append({"omega_x": omega_x, "components": comps, "strict_samples": sample_rows})
    payload["noncommuting_transverse_field"] = noncommuting

    # End-to-end generator-level PI benchmark.
    adaptive = []
    scenarios = [
        ("boundary", 0.0, 0.0, 2.0),
        ("buffered", 0.001, 0.001, 2.0),
        ("near_saddle", 0.0, 0.0, 1.64),
    ]
    for scenario, theta, kappa, varpi in scenarios:
        for policy in ("error_only", "certified_guard", "rotating_frame"):
            result = run_adaptive_channel_benchmark(
                final_time=2.0,
                initial_h=1.4,
                tolerance=0.002,
                theta=theta,
                kappa=kappa,
                varpi=varpi,
                policy=policy,
                fallback="rotating_frame",
            )
            adaptive.append(
                {
                    "scenario": scenario,
                    "theta": theta,
                    "kappa": kappa,
                    "varpi": varpi,
                    "policy": policy,
                    "global_normalized_choi_error": result.global_normalized_choi_error,
                    "minimum_step_choi_eigenvalue": result.min_step_choi_eigenvalue,
                    "stats": result.stats,
                    "steps": [
                        {
                            "t": row.t,
                            "proposed_h": row.proposed_h,
                            "used_h": row.used_h,
                            "error_estimate": row.error_estimate,
                            "cp_status": row.cp_status,
                            "minimum_choi_eigenvalue": row.min_choi_eigenvalue,
                            "action": row.action,
                        }
                        for row in result.records
                    ],
                }
            )
    payload["adaptive_benchmarks"] = adaptive

    (OUT / "v27_verification.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    # Flat CSVs for plotting/manuscript tables.
    with (OUT / "conditioning.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=conditioning[0].keys())
        writer.writeheader(); writer.writerows(conditioning)
    with (OUT / "adaptive_benchmarks.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["scenario", "theta", "kappa", "varpi", "policy", "global_normalized_choi_error", "minimum_step_choi_eigenvalue", "accepted_steps", "rejected_error", "rejected_cp", "fallback_steps", "rhs_stage_evaluations"]
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for row in adaptive:
            writer.writerow({
                "scenario": row["scenario"], "theta": row["theta"], "kappa": row["kappa"], "varpi": row["varpi"], "policy": row["policy"],
                "global_normalized_choi_error": row["global_normalized_choi_error"], "minimum_step_choi_eigenvalue": row["minimum_step_choi_eigenvalue"],
                "accepted_steps": row["stats"]["accepted_steps"], "rejected_error": row["stats"]["rejected_error"], "rejected_cp": row["stats"]["rejected_cp"], "fallback_steps": row["stats"]["fallback_steps"], "rhs_stage_evaluations": row["stats"]["rhs_stage_evaluations"],
            })

    md = ["# Version 2.7 certified verification", "", "## Blocking guard defect", ""]
    for row in bug_rows:
        md.append(f"- x={row['x']}: exact M={row['exact_margin_decimal']}; legacy tolerance accepts={row['legacy_tol_accept']}; certified={row['certified_status']}.")
    md += ["", "## Halving threshold", "", f"- resultant: `{resultant}`", f"- varpi_1/2: {half['varpi_mid']}", "", "## Noncommuting persistence", ""]
    for row in noncommuting:
        md.append(f"- Omega_x/Gamma={row['omega_x']}: components {row['components']}")
    md += ["", "## Adaptive benchmark summary", ""]
    for row in adaptive:
        md.append(f"- {row['scenario']} / {row['policy']}: accepted={row['stats']['accepted_steps']}, min lambda={row['minimum_step_choi_eigenvalue']:.6g}, global Bell/Choi error={row['global_normalized_choi_error']:.6g}.")
    (OUT / "v27_verification.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"Wrote {OUT / 'v27_verification.json'}")
    print(f"Wrote {OUT / 'v27_verification.md'}")


if __name__ == "__main__":
    main()
