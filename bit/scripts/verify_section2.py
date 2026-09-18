from __future__ import annotations

from pathlib import Path

import sympy as sp


def main() -> None:
    gamma, gamma_phi = sp.symbols(
        "gamma gamma_phi", nonnegative=True, real=True
    )
    a, c = sp.symbols("a c", real=True)

    # Vectorization order: (E00, E01, E10, E11).
    L = sp.Matrix(
        [
            [0, 0, 0, gamma],
            [0, -(gamma / 2 + gamma_phi), 0, 0],
            [0, 0, -(gamma / 2 + gamma_phi), 0],
            [0, 0, 0, -gamma],
        ]
    )

    e00 = sp.Matrix([1, 0, 0, 0])
    e01 = sp.Matrix([0, 1, 0, 0])
    e10 = sp.Matrix([0, 0, 1, 0])
    e11_minus_e00 = sp.Matrix([-1, 0, 0, 1])

    checks = {
        "L(E00)=0": L * e00 == sp.zeros(4, 1),
        "L(E11-E00)=-gamma(E11-E00)": (
            L * e11_minus_e00 == -gamma * e11_minus_e00
        ),
        "L(E01)=-(gamma/2+gamma_phi)E01": (
            L * e01 == -(gamma / 2 + gamma_phi) * e01
        ),
        "L(E10)=-(gamma/2+gamma_phi)E10": (
            L * e10 == -(gamma / 2 + gamma_phi) * e10
        ),
    }

    T = sp.Matrix(
        [
            [1, 0, 0, 1 - a],
            [0, c, 0, 0],
            [0, 0, c, 0],
            [0, 0, 0, a],
        ]
    )
    trace_row = sp.Matrix([[1, 0, 0, 1]])
    checks["trace_row*T=trace_row"] = trace_row * T == trace_row
    checks["T(E00)=E00"] = T * e00 == e00
    checks["T(E11)=(1-a)E00+aE11"] = T * sp.Matrix([0, 0, 0, 1]) == sp.Matrix(
        [1 - a, 0, 0, a]
    )

    failed = [name for name, ok in checks.items() if not ok]
    if failed:
        raise AssertionError("Section 2 verification failed: " + "; ".join(failed))

    report = ["# Section 2 symbolic verification", ""]
    report.extend(f"- PASS: {name}" for name in checks)
    report.extend(
        [
            "",
            "The exact numerical superoperator in vectorization order "
            "(E00,E01,E10,E11) is",
            "",
            "```text",
            str(T),
            "```",
        ]
    )

    root = Path(__file__).resolve().parents[1]
    out = root / "results" / "section2" / "symbolic_verification.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
