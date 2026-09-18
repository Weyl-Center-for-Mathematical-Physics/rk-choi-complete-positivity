#!/usr/bin/env python3
"""Deterministic BIT Numerical Mathematics figure generator (Task 11 of the BIT revision).

Reads only the frozen v3.4 plotting data under results/jcp_figures_v34/, results/v27/, and results/v34/
(no scientific quantity is recomputed except dense plotting grids of exact closed-form polynomials, which
are recorded next to the outputs).  Produces Fig1--Fig3 (article) and FigS1--FigS4 (Online Resource 1) at
the final width of 119 mm with 8--10 pt lettering, colour-blind-safe colours with redundant styles, and
embedded fonts; writes a manifest, an overlap audit, and grayscale / deuteranopia proof renders.

Usage (from code_and_data):  python scripts/generate_bit_figures.py [--out DIR]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("SOURCE_DATE_EPOCH", "1787184000")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib-cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.text  # noqa: E402
import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import Collection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402

DATA = ROOT / "results" / "jcp_figures_v34"
OUT = ROOT / "figures_bit"
RESULTS = ROOT / "results" / "bit_revision" / "figures"

SOURCE_FILES = [
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
]

WIDTH_IN = 119.0 / 25.4
FIXED_DATE = dt.datetime(2026, 9, 13, tzinfo=dt.timezone.utc)
PDF_META = {"Title": "Complete-positivity regions of Runge-Kutta discretizations: article figures",
            "Author": "G. Blake Pierpoint, Olivier Bernard, Yichen Liu",
            "Creator": "BIT deterministic Matplotlib figure generator", "CreationDate": FIXED_DATE, "ModDate": FIXED_DATE}
PNG_META = {"Software": "BIT deterministic Matplotlib figure generator", "Creation Time": "2026-09-13T00:00:00Z"}

# Okabe--Ito palette
BLUE, ORANGE, GREEN, PURPLE, SKY, YELLOW, VERMIL = "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#D55E00"
INK, NEUTRAL, GRID = "#222222", "#666666", "#E3E3E3"
LIGHT_BLUE, LIGHT_ORANGE = "#D9EDF7", "#FBE3BD"

EXPORT_RC = {
    "pdf.fonttype": 42,
    # Type 3 glyph programs survive EPS->PDF conversion in common pipelines; every glyph is embedded as vector outlines.
    "ps.fonttype": 3,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8.5,
    "mathtext.fontset": "stixsans",
    "axes.labelsize": 8.5,
    "axes.titlesize": 9.5,
    "axes.linewidth": 0.7,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "legend.fontsize": 8.0,
    "legend.handlelength": 2.0,
    "legend.handletextpad": 0.5,
    "legend.columnspacing": 0.9,
    "legend.framealpha": 1.0,
    "legend.borderpad": 0.3,
    "lines.linewidth": 1.1,
    "lines.markersize": 4.5,
    "hatch.linewidth": 0.6,
    "figure.dpi": 100,
    "savefig.dpi": 300,
}
ALPHA4 = 2.785293563405282
VARPI_M = 0.8273622514511260
VARPI_R = np.sqrt(3) / 2
VARPI_P = 1.6307120238943857
VARPI_C = 0.7950166505008560


# ---------------------------------------------------------------------------
# data access
# ---------------------------------------------------------------------------
def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(name: str) -> list[dict]:
    with (DATA / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def col(rows, key, cast=float):
    return np.array([cast(r[key]) if r[key] not in ("", "nan") else np.nan for r in rows])


def load_json(rel: str):
    return json.loads((ROOT / rel).read_text(encoding="utf-8"))


def source_hashes() -> dict[str, str]:
    return {rel: sha256(ROOT / rel) for rel in SOURCE_FILES}


def _cubic_roots(varpi: float) -> list[float]:
    """Companion roots of the exact RK4 cubic, dense rendering only (endpoints are certified elsewhere)."""
    q = 1 + 4 * varpi**2
    roots = np.roots([1.0, -16.0, -32.0 * (q - 6.0), 384.0 * (q - 4.0)])
    return sorted(float(r.real / q) for r in roots if abs(float(r.imag)) < 1e-7 and r.real > 1e-10)


def _branches(vgrid):
    """Attached upper / detached lower / detached upper edges of A_4 on a frequency grid (plot-only)."""
    qm = (123 - 11 * np.sqrt(33)) / 16
    qp = (123 + 11 * np.sqrt(33)) / 16
    au = np.full_like(vgrid, np.nan)
    dl = np.full_like(vgrid, np.nan)
    du = np.full_like(vgrid, np.nan)
    for i, v in enumerate(vgrid):
        q = 1 + 4 * v**2
        roots = _cubic_roots(float(v))
        eps = 2e-7
        if q < qm - eps and roots:
            au[i] = min(ALPHA4, roots[-1])
        elif abs(q - qm) <= eps:
            au[i] = min(ALPHA4, (2 * np.sqrt(33) - 2) / qm)
        elif q < 4 - eps and len(roots) == 3:
            au[i] = roots[0]
            dl[i], du[i] = roots[1], min(ALPHA4, roots[2])
        elif q > qp + eps and len(roots) == 2:
            dl[i], du[i] = roots[0], min(ALPHA4, roots[1])
    return au, dl, du


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def panel(ax, label):
    ax.set_title(label, loc="left", fontsize=9.5, fontweight="bold", pad=4)


def fill_components(ax, vgrid, au, dl, du, labels=True):
    ax.fill_between(vgrid, 0, au, where=np.isfinite(au), facecolor=LIGHT_BLUE, edgecolor=BLUE, linewidth=0.5,
                    hatch="///", label="attached component" if labels else None)
    ax.fill_between(vgrid, dl, du, where=np.isfinite(dl), facecolor=LIGHT_ORANGE, edgecolor=ORANGE, linewidth=0.5,
                    hatch="\\\\\\", label="detached component" if labels else None)
    ax.plot(vgrid, au, color=BLUE, linestyle="-", linewidth=1.1)
    ax.plot(vgrid, dl, color=ORANGE, linestyle="--", linewidth=1.1)
    ax.plot(vgrid, du, color=ORANGE, linestyle="--", linewidth=1.1)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def fig1():
    rows = read_csv("figure1_exact_branches.csv")
    vgrid, au, dl, du = col(rows, "varpi"), col(rows, "attached_upper"), col(rows, "detached_lower"), col(rows, "detached_upper")
    fig, ax = plt.subplots(figsize=(WIDTH_IN, 3.1), layout="constrained")
    fill_components(ax, vgrid, au, dl, du)
    ax.axhline(ALPHA4, color=INK, linestyle="-.", linewidth=0.9, label=r"population ceiling $\alpha_4$")
    for v in (VARPI_C, VARPI_M, VARPI_R, VARPI_P):
        ax.axvline(v, color=NEUTRAL, linestyle=":", linewidth=0.8)
    ax.plot([VARPI_R, VARPI_P], [2.0, (9 + np.sqrt(33)) / ((123 + 11 * np.sqrt(33)) / 16)], marker="D", markersize=5,
            markerfacecolor=YELLOW, markeredgecolor=INK, markeredgewidth=0.8, linestyle="None", label="isolated admissible step")
    ax.text(VARPI_P + 0.04, 0.12, r"$\varpi_+$", fontsize=8.5, ha="left", va="bottom")
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 3.4)
    ax.set_xlabel(r"frequency ratio $\varpi=|\omega|/\Gamma$")
    ax.set_ylabel(r"step $x=\Gamma h$")
    ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, borderaxespad=0.0)
    # zoom into the reentrant band
    inset = ax.inset_axes([0.60, 0.45, 0.36, 0.36])
    vz = np.linspace(0.78, 0.88, 1201)
    azu, zdl, zdu = _branches(vz)
    fill_components(inset, vz, azu, zdl, zdu, labels=False)
    for v, lab in ((VARPI_C, r"$\varpi_c$"), (VARPI_M, r"$\varpi_-$"), (VARPI_R, r"$\sqrt{3}/2$")):
        inset.axvline(v, color=NEUTRAL, linestyle=":", linewidth=0.7)
        inset.text(v + 0.002, 2.85, lab, fontsize=8.0, rotation=90, ha="left", va="top")
    inset.set_xlim(0.78, 0.88)
    inset.set_ylim(0, 2.95)
    inset.set_xticks([0.78, 0.83, 0.88])
    inset.set_yticks([0, 1, 2])
    inset.tick_params(labelsize=8.0)
    inset.set_facecolor("white")
    for spine in inset.spines.values():
        spine.set_linewidth(0.8)
    ax.indicate_inset([0.78, 0.0, 0.10, 2.95], inset_ax=None, edgecolor=NEUTRAL, linestyle=(0, (2, 2)), linewidth=0.8, alpha=1.0)
    with (RESULTS / "fig1_zoom_grid.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["varpi", "attached_upper", "detached_lower", "detached_upper"])
        w.writerows(zip(vz, azu, zdl, zdu))
    return fig


def conditioning_panel(bx):
    cond = read_csv("figure1_conditioning.csv")
    deltas, widths = col(cond, "delta_varpi"), col(cond, "window_width")
    ok = np.isfinite(widths)
    bx.loglog(deltas[ok], widths[ok], color=BLUE, marker="o", markevery=12, linestyle="-", label=r"$x_+-x_-$")
    ref = np.sqrt(deltas[ok]) * widths[ok][-1] / np.sqrt(deltas[ok][-1])
    bx.loglog(deltas[ok], ref, color=ORANGE, linestyle="--", label=r"$O((\varpi-\varpi_+)^{1/2})$")
    bx.set_xlabel(r"$\varpi-\varpi_+$")
    bx.set_ylabel("detached-window width")
    bx.grid(which="major", color=GRID, linewidth=0.5)
    bx.legend(frameon=False, loc="lower right")


def fig2():
    V = load_json("results/v34/v34_verification.json")
    direct = V["candidate_regions_varpi_2"]["full"][0]
    fine = V["candidate_regions_varpi_2"]["two_half"][0]
    marg = read_csv("figure2_richardson_margin.csv")
    xs, scaled = col(marg, "x"), col(marg, "M_ext_over_x6")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH_IN, 2.45), layout="constrained")
    panel(ax, "(a)")
    rows = [(1, direct["lower_decimal"], direct["upper_decimal"], BLUE, "o"), (0, fine["lower_decimal"], fine["upper_decimal"], ORANGE, "s")]
    for y, lo, hi, color, mk in rows:
        ax.hlines(y, lo, hi, linewidth=5, color=color, linestyles="-")
        ax.plot([lo, hi], [y, y], color=color, marker=mk, markerfacecolor="white" if mk == "s" else color, markeredgewidth=0.9, linestyle="None")
    ax.plot([0, 0], [1, 0], marker="o", markersize=3.5, color=INK, linestyle="None")
    for xv in (1.0, 2.0):
        ax.axvline(xv, color=NEUTRAL, linestyle=":", linewidth=0.9)
    ax.text(1.0, 1.55, r"$x=1$", ha="center", va="bottom", fontsize=8.5)
    ax.text(2.0, 1.55, r"$x=2$", ha="center", va="bottom", fontsize=8.5)
    ax.set_yticks([1, 0])
    ax.set_yticklabels(["one full step", "two half steps"])
    ax.set_xlim(-0.05, 2.75)
    ax.set_ylim(-0.6, 2.05)
    ax.set_xlabel(r"total step $x=\Gamma h$ at $\varpi=2$")
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    panel(bx, "(b)")
    bx.plot(xs, scaled, color=PURPLE, linestyle="-", label=r"$M_{\rm ext}(x)/x^6$")
    bx.axhline(0, color=INK, linewidth=0.8)
    bx.fill_between(xs, scaled, 0, where=scaled < 0, facecolor=LIGHT_ORANGE, edgecolor=ORANGE, linewidth=0.5, hatch="xx", label="extrapolate non-CPTP")
    bx.axvline(ALPHA4, color=NEUTRAL, linestyle=":", linewidth=0.9)
    bx.text(ALPHA4 - 0.05, -2.45e-4, r"$\alpha_4$", ha="right", va="bottom", fontsize=8.5)
    bx.set_xlim(0, 2.85)
    bx.set_ylim(-2.6e-4, 1.4e-4)
    bx.set_xlabel(r"nonrotating step $x$")
    bx.set_ylabel("scaled Choi margin")
    bx.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)
    bx.grid(color=GRID, linewidth=0.5)
    bx.legend(frameon=False, loc="upper right", borderaxespad=0.6, handlelength=1.4, handletextpad=0.4)
    return fig


def fig3():
    V = load_json("results/v27/v27_verification.json")
    buffered = V["buffered_components"]
    nc = read_csv("figure3_noncommuting_endpoints.csv")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH_IN, 2.55), layout="constrained", gridspec_kw={"width_ratios": [1.1, 1.0]})
    panel(ax, "(a)")
    ypos = np.arange(len(buffered))[::-1]
    labels = []
    for y, row in zip(ypos, buffered):
        th, ka = row["theta"], row["kappa"]
        labels.append(fr"$\theta={th}$, $\kappa={ka}$".replace("0.001", "10^{-3}").replace("0.01", "10^{-2}"))
        for j, (lo, hi) in enumerate(row["components"]):
            color, mk = (BLUE, "o") if j == 0 else (ORANGE, "s")
            ax.hlines(y, float(lo), float(hi), linewidth=3.2, color=color, linestyles="-")
            ax.plot([float(lo), float(hi)], [y, y], marker=mk, linestyle="None", markersize=4.0, color=color,
                    markerfacecolor="white" if j else color, markeredgewidth=0.8)
    ax.plot([], [], linewidth=3.2, color=BLUE, linestyle="-", marker="o", label="attached")
    ax.plot([], [], linewidth=3.2, color=ORANGE, linestyle="-", marker="s", markerfacecolor="white", label="detached")
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_xlim(-0.02, 1.36)
    ax.set_ylim(-0.7, 6.0)
    ax.set_xlabel(r"$x=\Gamma h$ at $\varpi=2$")
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    ax.legend(frameon=False, loc="upper center", ncol=2, handlelength=1.3, columnspacing=0.8, handletextpad=0.4)
    panel(bx, "(b)")
    ox = col(nc, "omega_x")
    bx.plot(ox, col(nc, "attached_upper"), color=BLUE, marker="o", linestyle="None", markersize=5, label="attached")
    bx.plot(ox, col(nc, "detached_lower"), color=ORANGE, marker="s", markerfacecolor="white", markeredgewidth=1.0, linestyle="None", markersize=5, label="detached, lower")
    bx.plot(ox, col(nc, "detached_upper"), color=GREEN, marker="^", markerfacecolor="white", markeredgewidth=1.0, linestyle="None", markersize=5.5, label="detached, upper")
    bx.set_xticks(ox)
    bx.set_xticklabels(["0", "0.05", "0.1", "0.2"])
    bx.set_xlim(-0.03, 0.235)
    bx.set_ylim(0.15, 1.40)
    bx.set_xlabel(r"certified coupling $\Omega_x/\Gamma$")
    bx.set_ylabel(r"certified endpoint $x$")
    bx.grid(color=GRID, linewidth=0.5)
    bx.legend(frameon=False, loc="center", bbox_to_anchor=(0.5, 0.28), handletextpad=0.3, handlelength=1.2)
    return fig


def fig4():
    sweep = read_csv("figure5_tolerance_sweep.csv")
    comp = read_csv("figure5_action_composition.csv")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH_IN, 2.5), layout="constrained", gridspec_kw={"width_ratios": [1.0, 1.2]})
    panel(ax, "(a)")
    styles = {"error_only": (BLUE, "o", "-", "error only"), "guard_rotating_fallback": (ORANGE, "s", "--", "guard + fallback"),
              "rotating_frame": (GREEN, "^", "-.", "rotating frame")}
    for policy, (color, mk, ls, label) in styles.items():
        data = sorted([r for r in sweep if r["policy"] == policy], key=lambda r: float(r["tolerance"]), reverse=True)
        work = np.array([float(r["rhs_stage_evaluations"]) + 3 * float(r["exponential_actions"]) for r in data])
        err = np.array([float(r["global_normalized_choi_error"]) for r in data])
        ax.loglog(work, err, color=color, marker=mk, markerfacecolor="white" if policy != "error_only" else color, markeredgewidth=0.9, linestyle=ls, label=label)
    ax.set_xlabel("propagator work")
    ax.set_ylabel("global normalized Choi error")
    ax.set_ylim(1e-6, 2e-1)
    ax.set_xlim(22, 260)
    ax.set_xticks([30, 60, 100, 200])
    ax.set_xticklabels(["30", "60", "100", "200"])
    ax.set_xticks([], minor=True)
    ax.grid(which="major", color=GRID, linewidth=0.5)
    ax.legend(frameon=False, loc="upper left")
    panel(bx, "(b)")
    order = ["buffered_projection", "buffered_direct", "near_saddle_projection", "boundary"]
    names = ["buffered\n(1)", "buffered\n(2)", "near\nsaddle", "boundary"]
    rows = []
    for key in order[:3]:
        rows.append(next(r for r in comp if r["scenario"] == key))
    rows.append(next(r for r in sweep if r["policy"] == "guard_rotating_fallback" and abs(float(r["tolerance"]) - 1e-3) < 1e-15))
    direct = np.array([float(r["direct_accepts"]) for r in rows])
    projected = np.array([float(r["projected_accepts"]) for r in rows])
    rotating = np.array([float(r["rotating_accepts"]) for r in rows])
    exact = np.array([float(r["exact_fallback_accepts"]) for r in rows])
    pos = np.arange(len(rows))
    for heights, bottoms, label, color, hatch in ((direct, 0 * direct, "direct", BLUE, "///"), (projected, direct, "projected", ORANGE, "\\\\\\"),
                                                  (rotating, direct + projected, "rotating fallback", GREEN, "xx"), (exact, direct + projected + rotating, "exact fallback", PURPLE, "..")):
        if np.any(heights):
            bx.bar(pos, heights, bottom=bottoms, label=label, color=color, edgecolor=INK, linewidth=0.5, hatch=hatch, width=0.62)
    bx.set_xticks(pos)
    bx.set_xticklabels(names)
    bx.set_ylabel("accepted candidates")
    bx.set_ylim(0, 8.2)
    bx.grid(axis="y", color=GRID, linewidth=0.5)
    bx.legend(frameon=False, loc="upper right", ncol=1)
    return fig


def figS1():
    fig, ax = plt.subplots(figsize=(WIDTH_IN, 1.95), layout="constrained")
    vmax, vlo, vr, vhi = 2.3, np.sqrt(2 - np.sqrt(33) / 4), VARPI_R, np.sqrt(2 + np.sqrt(33) / 4)
    rows = [("Richardson extrapolate of RK4", [(0, vlo, "out"), (vlo, vhi, "in"), (vhi, vmax, "out")]),
            ("Dormand–Prince 5", [(0, vlo, "in"), (vlo, vhi, "out"), (vhi, vmax, "in")]),
            ("classical RK4", [(0, vr, "in"), (vr, vmax, "out")])]
    for y, (name, segs) in enumerate(rows):
        for lo, hi, state in segs:
            ax.barh(y, hi - lo, left=lo, height=0.55, facecolor=LIGHT_BLUE if state == "in" else LIGHT_ORANGE,
                    edgecolor=BLUE if state == "in" else ORANGE, hatch="" if state == "in" else "///", linewidth=0.6)
            ax.text((lo + hi) / 2, y, "CPTP" if state == "in" else "non-CPTP", ha="center", va="center", fontsize=8.0, gid="inbar")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows])
    for v in (vlo, vr, vhi):
        ax.axvline(v, color=NEUTRAL, linestyle=":", linewidth=0.8)
    ax.text(vlo - 0.02, 2.62, "0.751", ha="right", va="bottom", fontsize=8.0)
    ax.text(vr + 0.02, 2.62, "0.866", ha="left", va="bottom", fontsize=8.0)
    ax.text(vhi, 2.62, "1.854", ha="center", va="bottom", fontsize=8.0)
    ax.set_xlim(0, vmax)
    ax.set_ylim(-0.5, 3.15)
    ax.set_xlabel(r"frequency ratio $\varpi=|\omega|/\Gamma$")
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    return fig


def figS2():
    V = load_json("results/v27/v27_verification.json")
    rows = read_csv("figure4_conditioning.csv")
    x = col(rows, "x")
    exact = col(rows, "exact_margin")
    dbl = col(rows, "binary64_margin")
    from fractions import Fraction

    lo = [float(Fraction(v)) for v in V["varpi2_detached_interval"]["x_minus_bracket"]]
    hi = [float(Fraction(v)) for v in V["varpi2_detached_interval"]["x_plus_bracket"]]
    xminus, xplus = sum(lo) / 2, sum(hi) / 2
    fig = plt.figure(figsize=(WIDTH_IN, 4.9), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05])
    cx = fig.add_subplot(gs[0, :])
    ax = fig.add_subplot(gs[1, 0])
    bx = fig.add_subplot(gs[1, 1])
    panel(cx, "(a)")
    conditioning_panel(cx)
    panel(ax, "(b)")
    xs = np.linspace(0.02, 1.45, 900)
    q = 17.0
    signed = -((q * xs) ** 3 - 16 * (q * xs) ** 2 - 32 * (q - 6) * (q * xs) + 384 * (q - 4))
    ax.plot(xs, signed, color=BLUE, linestyle="-", label=r"$-\mathcal{C}_{17}(17x)$")
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.axvspan(xminus, xplus, facecolor=LIGHT_ORANGE, edgecolor=ORANGE, linewidth=0.5, hatch="///", label="admissible")
    seq = np.array([1.4 / 2**j for j in range(6)])
    seq_y = -((q * seq) ** 3 - 16 * (q * seq) ** 2 - 32 * (q - 6) * (q * seq) + 384 * (q - 4))
    ax.plot(seq, seq_y, color=PURPLE, marker="o", markerfacecolor="white", markeredgewidth=0.9, linestyle="None", markersize=4.5, label="halving steps")
    ax.set_yscale("symlog", linthresh=100)
    ax.set_xlabel(r"$x=\Gamma h$ at $\varpi=2$")
    ax.set_ylabel("signed polynomial (symlog)")
    ax.grid(color=GRID, linewidth=0.5)
    ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=1, borderaxespad=0.4, handlelength=1.4, handletextpad=0.4)
    panel(bx, "(c)")
    direct = np.abs(dbl)
    direct[direct == 0] = np.nan
    bx.loglog(x, np.abs(exact), color=BLUE, marker="o", linestyle="-", label=r"exact $|M_4|$")
    bx.loglog(x, direct, color=ORANGE, marker="s", markerfacecolor="white", markeredgewidth=0.9, linestyle="--", label="binary64")
    bx.axhline(1e-12, color=NEUTRAL, linestyle=":", linewidth=0.9, label=r"$10^{-12}$")
    false = np.array([0.004, 0.001, 1e-4, 1e-5])
    false_y = [abs(float(r["exact_margin_decimal"])) for r in V["blocking_false_acceptance_reproduction"]]
    bx.scatter(false, false_y, color=PURPLE, marker="x", s=32, linewidths=1.1, label="false accept")
    bx.set_xlabel(r"$x$")
    bx.set_ylabel("margin magnitude")
    bx.grid(which="major", color=GRID, linewidth=0.5)
    bx.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=1, borderaxespad=0.4, handlelength=1.4, handletextpad=0.4)
    return fig


def figS3():
    rows = read_csv("figureS2_operational.csv")
    vscan, ratios, absdef = col(rows, "varpi"), col(rows, "ratio"), col(rows, "absolute_defect")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH_IN, 2.8), layout="constrained")
    panel(ax, "(a)")
    l1 = ax.plot(vscan, ratios, color=BLUE, linestyle="-", label=r"$\delta_{\rm CP}/E_J$")
    ax.set_xlabel(r"$\varpi$ at $x=10^{-2}$")
    ax.set_ylabel("fraction of Bell-input trace error")
    ax.set_ylim(0, 1)
    ax.grid(color=GRID, linewidth=0.5)
    ax2 = ax.twinx()
    l2 = ax2.semilogy(vscan, absdef, color=ORANGE, linestyle="--", label=r"absolute $\delta_{\rm CP}$")
    ax2.set_ylabel("absolute negative mass")
    ax.legend(l1 + l2, [l.get_label() for l in l1 + l2], frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=1, borderaxespad=0.4, handlelength=1.6, handletextpad=0.4)
    panel(bx, "(b)")
    mp.mp.dps = 60
    Ns = np.unique(np.logspace(1, 5, 60).astype(int))
    for name, fun, color, mk, ls in (("Forward Euler", lambda z: 1 + z, BLUE, "o", "-"), ("SSPRK(3,3)", lambda z: 1 + z + z**2 / 2 + z**3 / 6, ORANGE, "s", "--")):
        vals = []
        for N in Ns:
            h = mp.mpf(1) / int(N)
            A = fun(-h) ** int(N)
            C = fun(-h / 2) ** int(N)
            lam = (1 + A - mp.sqrt((1 - A) ** 2 + 4 * abs(C) ** 2)) / 4
            vals.append(float(-lam))
        bx.loglog(1 / Ns, vals, color=color, marker=mk, markevery=8, markersize=3.5, linestyle=ls, label=name)
    bx.set_xlabel(r"step $x_N=1/N$")
    bx.set_ylabel(r"$-p^{\rm RK}_{-,N}$")
    bx.grid(which="major", color=GRID, linewidth=0.5)
    bx.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=1, borderaxespad=0.4, handlelength=1.6, handletextpad=0.4)
    return fig


BUILDERS = {"Fig1": fig1, "Fig2": fig2, "Fig3": fig3, "FigS1": figS1, "FigS2": figS2, "FigS3": figS3, "FigS4": fig4}


def build(stem: str):
    RESULTS.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(EXPORT_RC):
        return BUILDERS[stem]()


def close(fig):
    plt.close(fig)


def all_text_artists(fig):
    return [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible()]


# ---------------------------------------------------------------------------
# overlap audit (rendered geometry)
# ---------------------------------------------------------------------------
def _data_points(ax):
    """Display-space points and rectangles of the data artists of an axes (interior points only, so that
    axvline/axhline endpoints on the frame do not count)."""
    pts = []
    rects = []
    for line in ax.lines:
        xy = np.column_stack([np.asarray(line.get_xdata(), dtype=float), np.asarray(line.get_ydata(), dtype=float)])
        xy = xy[np.isfinite(xy).all(axis=1)]
        if len(xy):
            pts.append(line.get_transform().transform(xy))
    for coll in ax.collections:
        if not isinstance(coll, Collection):
            continue
        tr = coll.get_transform()
        for path in coll.get_paths():
            v = np.asarray(path.vertices, dtype=float)
            v = v[np.isfinite(v).all(axis=1)]
            if len(v):
                pts.append(tr.transform(v))
    for patch in ax.patches:
        if isinstance(patch, Rectangle):
            rects.append(patch.get_window_extent())
    allpts = np.vstack(pts) if pts else np.zeros((0, 2))
    if len(allpts):
        ab = ax.get_window_extent()
        inside = (allpts[:, 0] > ab.x0 + 1.5) & (allpts[:, 0] < ab.x1 - 1.5) & (allpts[:, 1] > ab.y0 + 1.5) & (allpts[:, 1] < ab.y1 - 1.5)
        allpts = allpts[inside]
    return allpts, rects


def _bbox_hits(bbox: Bbox, pts, rects, pad=1.0):
    b = Bbox.from_extents(bbox.x0 - pad, bbox.y0 - pad, bbox.x1 + pad, bbox.y1 + pad)
    hit = False
    if len(pts):
        inside = (pts[:, 0] >= b.x0) & (pts[:, 0] <= b.x1) & (pts[:, 1] >= b.y0) & (pts[:, 1] <= b.y1)
        hit = bool(inside.any())
    for r in rects:
        if b.overlaps(r) and r.width > 0 and r.height > 0:
            hit = True
    return hit


def overlap_audit(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    report = {"legend_data_overlaps": [], "text_data_overlaps": [], "text_text_overlaps": [], "min_font_pt": None}
    boxes = []
    for ax in fig.get_axes():
        pts, rects = _data_points(ax)
        leg = ax.get_legend()
        if leg is not None:
            lb = leg.get_window_extent(renderer)
            boxes.append(("legend@" + repr(ax.get_position().bounds), lb))
            if _bbox_hits(lb, pts, rects):
                report["legend_data_overlaps"].append(f"legend of axes {ax.get_position().bounds}")
            ab = ax.get_window_extent(renderer)
            inside = lb.x0 >= ab.x0 - 0.5 and lb.x1 <= ab.x1 + 0.5 and lb.y0 >= ab.y0 - 0.5 and lb.y1 <= ab.y1 + 0.5
            outside = lb.x1 <= ab.x0 + 0.5 or lb.x0 >= ab.x1 - 0.5 or lb.y1 <= ab.y0 + 0.5 or lb.y0 >= ab.y1 - 0.5
            if not inside and not outside:
                report["legend_data_overlaps"].append(f"legend crosses the axes frame of axes {ax.get_position().bounds}")
        # panel labels are left-aligned titles; include them so legends placed above an axes cannot collide with them
        for t in (ax.title, ax._left_title, ax._right_title):
            if t.get_text().strip():
                boxes.append(("title:" + t.get_text(), t.get_window_extent(renderer)))
        for t in ax.texts:
            if not t.get_text().strip() or t.get_gid() == "inbar":
                continue
            tb = t.get_window_extent(renderer)
            boxes.append((t.get_text(), tb))
            if _bbox_hits(tb, pts, rects, pad=0.5):
                report["text_data_overlaps"].append(t.get_text())
        # child inset axes as opaque boxes against the parent's data
        for child in ax.child_axes:
            cb = child.get_window_extent(renderer)
            if _bbox_hits(cb, pts, rects, pad=0.0):
                report["text_data_overlaps"].append("inset axes")
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if boxes[i][1].overlaps(boxes[j][1]):
                report["text_text_overlaps"].append((boxes[i][0], boxes[j][0]))
    sizes = [t.get_fontsize() for t in all_text_artists(fig) if t.get_text().strip()]
    report["min_font_pt"] = min(sizes) if sizes else None
    report["max_font_pt"] = max(sizes) if sizes else None
    return report


def series_counts():
    out = {}
    f = build("Fig1")
    a = f.get_axes()[0]
    out["Fig1"] = {"a_fills": sum(1 for c in a.collections if c.get_label() and not c.get_label().startswith("_")), "inset": len(a.child_axes)}
    close(f)
    f = build("Fig2")
    a, b = f.get_axes()
    out["Fig2"] = {"a_intervals": len(a.collections), "b_lines": len([l for l in b.lines if len(l.get_xdata()) > 2])}
    close(f)
    f = build("Fig3")
    a, b = f.get_axes()
    out["Fig3"] = {"a_rows": len(a.get_yticks()), "b_marker_series": len([l for l in b.lines if l.get_linestyle() == "None"]),
                   "b_connecting_lines": len([l for l in b.lines if l.get_linestyle() != "None"])}
    close(f)
    f = build("FigS4")
    a, b = f.get_axes()
    out["FigS4"] = {"a_series": len([l for l in a.lines if l.get_label() and not l.get_label().startswith("_")]), "b_categories": len(b.get_xticks())}
    close(f)
    return out


# ---------------------------------------------------------------------------
# rendering, proofs, manifest
# ---------------------------------------------------------------------------
def save_all(fig, directory: Path, stem: str):
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / f"{stem}.pdf", metadata=PDF_META)
    fig.savefig(directory / f"{stem}.eps")
    fig.savefig(directory / f"{stem}.png", dpi=300, metadata=PNG_META)


def render_to(directory: Path, stem: str) -> Path:
    fig = build(stem)
    save_all(fig, Path(directory), stem)
    close(fig)
    return Path(directory)


DEUTERANOPIA = np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]])


def proofs(png: Path, audit_dir: Path):
    from matplotlib.image import imread, imsave

    img = imread(png)[..., :3].astype(float)
    lin = np.where(img <= 0.04045, img / 12.92, ((img + 0.055) / 1.055) ** 2.4)
    gray = 0.2126 * lin[..., 0] + 0.7152 * lin[..., 1] + 0.0722 * lin[..., 2]
    gray_s = np.where(gray <= 0.0031308, 12.92 * gray, 1.055 * gray ** (1 / 2.4) - 0.055)
    imsave(audit_dir / f"{png.stem}_grayscale.png", np.clip(gray_s, 0, 1), cmap="gray", metadata=PNG_META)
    deut = np.clip(lin @ DEUTERANOPIA.T, 0, 1)
    deut_s = np.where(deut <= 0.0031308, 12.92 * deut, 1.055 * deut ** (1 / 2.4) - 0.055)
    imsave(audit_dir / f"{png.stem}_deuteranopia.png", np.clip(deut_s, 0, 1), metadata=PNG_META)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(OUT))
    args = parser.parse_args(argv)
    out = Path(args.out)
    audit_dir = out / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    manifest = {"generator": "scripts/generate_bit_figures.py", "width_mm": 119.0, "sources": source_hashes(), "outputs": {}, "audit": {}}
    for stem in BUILDERS:
        fig = build(stem)
        w, h = fig.get_size_inches()
        rep = overlap_audit(fig)
        rep["size_mm"] = [round(w * 25.4, 2), round(h * 25.4, 2)]
        manifest["audit"][stem] = rep
        save_all(fig, out, stem)
        close(fig)
        proofs(out / f"{stem}.png", audit_dir)
        for ext in ("eps", "pdf", "png"):
            p = out / f"{stem}.{ext}"
            manifest["outputs"][p.name] = {"bytes": p.stat().st_size, "sha256": sha256(p)}
        print(f"{stem}: {rep['size_mm'][0]} x {rep['size_mm'][1]} mm, fonts {rep['min_font_pt']}-{rep['max_font_pt']} pt, "
              f"legend/data overlaps {len(rep['legend_data_overlaps'])}, text/data {len(rep['text_data_overlaps'])}, text/text {len(rep['text_text_overlaps'])}")
    (RESULTS / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("wrote", RESULTS / "manifest.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
