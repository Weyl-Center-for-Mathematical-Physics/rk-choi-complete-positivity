#!/usr/bin/env python3
"""Generate the three figures used in the benchmark submission.

The script uses only exact endpoint values already certified by the symbolic audit
and high-precision repeated-map formulas. It writes PDF, EPS, PNG, and the source
CSV tables used by the plots.
"""
from __future__ import annotations

import csv
import datetime as dt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
# Embed TrueType outlines rather than Type 3 glyphs in PDF/EPS outputs.
# This improves portability, searchability, and journal-production compatibility.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures_benchmark"
DATA = ROOT / "results" / "benchmark_figures"
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

PDF_META = {
    "Title": "benchmark submission figure",
    "Author": "G. Blake Pierpoint, Olivier Bernard, and Yichen Liu",
    "Creator": "generate_benchmark_figures.py",
    "CreationDate": dt.datetime(2026, 8, 9, 0, 0, 0, tzinfo=dt.timezone.utc),
    "ModDate": dt.datetime(2026, 8, 9, 0, 0, 0, tzinfo=dt.timezone.utc),
}
PNG_META = {
    "Title": "benchmark submission figure",
    "Author": "G. Blake Pierpoint, Olivier Bernard, and Yichen Liu",
    "Software": "generate_benchmark_figures.py",
}


def save(fig: plt.Figure, stem: str) -> None:
    pdf_path = OUT / f"{stem}.pdf"
    eps_path = OUT / f"{stem}.eps"
    png_path = OUT / f"{stem}.png"
    fig.savefig(pdf_path, bbox_inches="tight", metadata=PDF_META)
    fig.savefig(eps_path, bbox_inches="tight", metadata={"Creator": "generate_benchmark_figures.py"})
    fig.savefig(png_path, dpi=600, bbox_inches="tight", metadata=PNG_META)
    plt.close(fig)

    # Matplotlib inserts the wall-clock time into EPS headers even when all
    # numerical and graphical content is deterministic.  Normalize it so a
    # clean reproduction run is byte-for-byte checkable.
    eps_text = eps_path.read_text(encoding="latin-1")
    eps_lines = [
        "%%CreationDate: Sun Aug  9 00:00:00 2026" if line.startswith("%%CreationDate:") else line
        for line in eps_text.splitlines()
    ]
    eps_path.write_text("\n".join(eps_lines) + "\n", encoding="latin-1", newline="\n")


def geometry() -> None:
    # Panel (a): determinant boundary of the autonomous phase-covariant Choi block.
    theta_values = [0.0, 0.25, 0.5]
    a_grid = np.linspace(-1.0, 1.0, 1200)
    rows: list[list[float | str | int]] = []

    v_rk4 = np.sqrt(3.0) / 2.0
    v_ssp3 = np.sqrt(7.0) / 2.0
    x_max = 2.0

    with (DATA / "figure1_channel_geometry.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["panel", "series", "kind", "x_left", "x_right", "y_value", "orientation", "admissible"])
        for theta in theta_values:
            lower = -min(theta / (1 - theta), (1 - theta) / theta) if 0.0 < theta < 1.0 else 0.0
            boundary = a_grid + theta * (1 - theta) * (1 - a_grid) ** 2
            for av, bv in zip(a_grid, boundary):
                admissible = int(av >= lower - 1e-15 and av <= 1.0 and bv >= -1e-15)
                rows.append(["a", f"theta={theta:g}", "curve", av, av, bv, "", admissible])
        # Piecewise orientation data, including the analytically resolved equality cases.
        rows.extend([
            ["b", "SSPRK(3,3)", "interval", 0.0, v_ssp3, 1.0, "outward", 1],
            ["b", "SSPRK(3,3)", "interval", v_ssp3, x_max, 1.0, "inward", 1],
            ["b", "SSPRK(3,3)", "threshold", v_ssp3, v_ssp3, 1.0, "inward", 1],
            ["b", "classical RK4", "interval", 0.0, v_rk4, 0.0, "inward", 1],
            ["b", "classical RK4", "interval", v_rk4, x_max, 0.0, "outward", 1],
            ["b", "classical RK4", "threshold", v_rk4, v_rk4, 0.0, "outward", 1],
        ])
        w.writerows(rows)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.05, 2.75))

    # Generalized-amplitude-damping Choi boundary.
    for theta, linestyle in zip(theta_values, ["-", "--", "-."]):
        lower = -min(theta / (1 - theta), (1 - theta) / theta) if 0.0 < theta < 1.0 else 0.0
        mask = a_grid >= lower
        boundary = a_grid + theta * (1 - theta) * (1 - a_grid) ** 2
        label = (r"$\theta=0$ / exact $u=0$ trajectory" if theta == 0.0
                 else rf"$\theta={theta:g}$")
        ax1.plot(a_grid[mask], boundary[mask], linestyle=linestyle, linewidth=1.55,
                 label=label)
    ax1.axvline(0.0, color="0.55", linewidth=0.7, linestyle=":")
    ax1.set_xlim(-1.0, 1.0)
    ax1.set_ylim(0.0, 1.05)
    ax1.set_xlabel(r"population-mode multiplier $a$")
    ax1.set_ylabel(r"maximum allowed $|c|^2$")
    ax1.set_title("(a) Phase-covariant Choi boundary", fontsize=8.5)
    ax1.tick_params(labelsize=7)
    ax1.grid(alpha=0.16, linewidth=0.55)
    ax1.legend(loc="lower right", fontsize=5.8, frameon=True, borderpad=0.3, handlelength=2.1)

    # Frequency-dependent orientation as a categorical phase diagram.
    method_specs = [
        (1.0, "SSPRK(3,3)", v_ssp3, "outward", "inward", "#2C7FB8"),
        (0.0, "classical RK4", v_rk4, "inward", "outward", "#D95F0E"),
    ]
    for y, label, threshold, left_orientation, right_orientation, color in method_specs:
        left_style = "-" if left_orientation == "inward" else (0, (5, 3))
        right_style = "-" if right_orientation == "inward" else (0, (5, 3))
        ax2.hlines(y, 0.0, threshold, color=color, linewidth=3.0, linestyles=left_style)
        ax2.hlines(y, threshold, x_max, color=color, linewidth=3.0, linestyles=right_style)
        ax2.plot(threshold, y, "o", markersize=5.0, markerfacecolor=color, markeredgecolor="white",
                 markeredgewidth=0.7, zorder=5)
        ax2.text(threshold / 2, y + 0.14, left_orientation, ha="center", va="bottom", fontsize=6.2)
        ax2.text((threshold + x_max) / 2, y + 0.14, right_orientation, ha="center", va="bottom", fontsize=6.2)
    ax2.axvline(v_rk4, color="0.65", linewidth=0.65, linestyle=":", zorder=0)
    ax2.axvline(v_ssp3, color="0.65", linewidth=0.65, linestyle=":", zorder=0)
    ax2.text(v_rk4, -0.34, r"$\sqrt{3}/2$", ha="center", va="top", fontsize=6.5)
    ax2.text(v_ssp3, 0.66, r"$\sqrt{7}/2$", ha="center", va="top", fontsize=6.5)
    ax2.set_xlim(0.0, x_max)
    ax2.set_ylim(-0.48, 1.48)
    ax2.set_yticks([1.0, 0.0])
    ax2.set_yticklabels([r"SSPRK(3,3)", r"classical RK4"], fontsize=7)
    ax2.set_xlabel(r"frequency ratio $|\varpi|=|\omega|/\Gamma$")
    ax2.set_title("(b) Small-step orientation", fontsize=8.5)
    ax2.tick_params(axis="x", labelsize=7)
    ax2.grid(axis="x", alpha=0.16, linewidth=0.55)
    ax2.spines["left"].set_visible(False)
    ax2.tick_params(axis="y", length=0)
    # Legend encodes orientation independently of color.
    from matplotlib.lines import Line2D
    ax2.legend(handles=[
        Line2D([0], [0], color="0.25", linewidth=2.5, linestyle="-", label="inward"),
        Line2D([0], [0], color="0.25", linewidth=2.5, linestyle=(0, (5, 3)), label="outward"),
    ], loc="lower right", fontsize=6.0, frameon=True, borderpad=0.3, handlelength=2.2)

    fig.tight_layout(pad=0.45, w_pad=1.25)
    save(fig, "Fig1_channel_geometry")


def admissible_sets() -> None:
    alpha4 = 2.785293563405282
    alpha5 = 3.306567892634947
    beta_minus = 3.068856548050381
    beta_plus = 4.384986320801944
    rows = [
        ("Forward Euler", "point", 0.0, 0.0),
        ("Heun RK2", "interval", 0.0, 2.0),
        ("SSPRK(3,3)", "point", 0.0, 0.0),
        ("Classical RK4", "interval", 0.0, alpha4),
        ("Dormand--Prince 5 (principal)", "interval", 0.0, alpha5),
        ("Dormand--Prince 4 (embedded)", "point", 0.0, 0.0),
        ("Dormand--Prince 4 (embedded)", "isolated_interval", beta_minus, beta_plus),
        ("Backward Euler", "unbounded", 0.0, float("inf")),
        ("Implicit midpoint", "point", 0.0, 0.0),
    ]
    with (DATA / "figure2_admissible_sets.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["method", "component_type", "left", "right"])
        for method, kind, left, right in rows:
            w.writerow([method, kind, left, "infinity" if np.isinf(right) else right])

    labels = ["Euler", "Heun RK2", "SSPRK3", "RK4", "DP5", "DP4 emb.", "Backward Euler", "Implicit midpoint"]
    y = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(3.42, 3.12))
    blue = "#2C7FB8"
    orange = "#E66101"
    dark = "0.20"

    for idx in [0, 2, 5, 7]:
        ax.plot(0, y[idx], "o", markersize=4.2, color=dark, zorder=4)
    for idx, endpoint in [(1, 2.0), (3, alpha4), (4, alpha5)]:
        ax.hlines(y[idx], 0, endpoint, color=blue, linewidth=2.4)
        ax.vlines([0, endpoint], y[idx] - 0.12, y[idx] + 0.12, color=blue, linewidth=1.0)
    ax.hlines(y[5], beta_minus, beta_plus, color=orange, linewidth=2.4, linestyles=(0, (5, 3)))
    ax.vlines([beta_minus, beta_plus], y[5] - 0.12, y[5] + 0.12, color=orange, linewidth=1.0)
    ax.annotate("", xy=(5.16, y[6]), xytext=(0, y[6]),
                arrowprops=dict(arrowstyle="->", color=dark, linewidth=2.4))

    ax.text(2.0, y[1] - 0.26, r"$2$", ha="center", va="center", fontsize=6.5)
    ax.text(alpha4, y[3] - 0.26, r"$\alpha_4$", ha="center", va="center", fontsize=6.5)
    ax.text(alpha5, y[4] + 0.25, r"$\alpha_5$", ha="center", va="center", fontsize=6.5)
    ax.text(beta_minus, y[5] + 0.25, r"$\beta_-$", ha="center", va="center", fontsize=6.5)
    ax.text(beta_plus, y[5] + 0.25, r"$\beta_+$", ha="center", va="center", fontsize=6.5)
    ax.text((beta_minus + beta_plus) / 2, y[5] - 0.28, "isolated", color=orange,
            ha="center", va="center", fontsize=6.2)
    ax.text(4.95, y[6] + 0.25, "unbounded", color=dark, ha="right", va="center", fontsize=6.2)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.7)
    ax.set_xlabel(r"relaxation step $x$", fontsize=7.5)
    ax.set_xlim(-0.10, 5.22)
    ax.set_ylim(-0.55, len(labels) - 0.45)
    ax.grid(axis="x", alpha=0.18, linewidth=0.55)
    ax.tick_params(axis="x", labelsize=6.7)
    ax.tick_params(axis="y", length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout(pad=0.25)
    save(fig, "Fig2_admissible_sets")


def r_fe(z: mp.mpf) -> mp.mpf:
    return 1 + z


def r_ssp3(z: mp.mpf) -> mp.mpf:
    return 1 + z + z**2 / 2 + z**3 / 6


def lambda_minus(a: mp.mpf, c: mp.mpf) -> mp.mpf:
    return (1 + a - mp.sqrt((1 - a) ** 2 + 4 * c**2)) / 2


def fixed_horizon() -> None:
    mp.mp.dps = 90
    X = mp.mpf("1")
    Ns = np.unique(np.logspace(np.log10(8), 4, 15).astype(int))
    configs = [
        ("Forward Euler", r_fe, 2, -mp.mpf(1) / 4, "o"),
        ("SSPRK(3,3)", r_ssp3, 4, -mp.mpf(7) / 192, "s"),
    ]
    records: list[list[str | int]] = []

    fig, ax = plt.subplots(figsize=(3.42, 2.78))
    colors = ["#2C7FB8", "#D95F0E"]
    for (label, R, m, K, marker), color in zip(configs, colors):
        hs: list[float] = []
        exact: list[float] = []
        asym: list[float] = []
        for Nint in Ns:
            N = mp.mpf(int(Nint))
            h = X / N
            a = R(-h) ** int(Nint)
            c = R(-h / 2) ** int(Nint)
            lam = lambda_minus(a, c)
            approx = (K * X * mp.e**(-X) / (1 + mp.e**(-X))) * h ** (m - 1)
            hs.append(float(h))
            exact.append(float(-lam))
            asym.append(float(-approx))
            records.append([label, int(Nint), mp.nstr(h, 30), mp.nstr(lam, 35), mp.nstr(approx, 35)])
        ax.loglog(hs, exact, marker=marker, markersize=3.3, linewidth=1.25, color=color,
                  label=f"{label}: exact")
        ax.loglog(hs, asym, linestyle="--", linewidth=1.15, color=color,
                  label=f"{label}: asymptotic")

    with (DATA / "figure3_fixed_horizon.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["method", "N", "x_N", "lambda_minus_exact", "lambda_minus_asymptotic"])
        w.writerows(records)

    ax.set_xlabel(r"step $x_N=X/N$ at fixed $X=1$")
    ax.set_ylabel(r"defect magnitude $-\lambda_{-,N}$")
    ax.tick_params(labelsize=7)
    ax.grid(which="both", alpha=0.18, linewidth=0.6)
    ax.legend(loc="upper left", fontsize=5.8, frameon=True, borderpad=0.3, handlelength=2.0)
    fig.tight_layout(pad=0.25)
    save(fig, "Fig3_fixed_horizon")


if __name__ == "__main__":
    geometry()
    admissible_sets()
    fixed_horizon()
    print(f"Wrote benchmark figures to {OUT}")
