#!/usr/bin/env python3
"""Deterministic BIT Numerical Mathematics figure generator (Task 11 of the BIT revision).

Reads only the frozen v3.4 plotting data under results/jcp_figures_v34/, results/v27/, and results/v34/
(no scientific quantity is recomputed except dense plotting grids of exact closed-form polynomials, which
are recorded next to the outputs).  Produces Fig1--Fig3 (article) and FigS1--FigS4 (Online Resource 1) at
the final width of 119 mm with 8--10 pt Computer Modern lettering (the article is typeset with pdflatex in
Computer Modern), the Okabe--Ito palette with redundant line styles, markers and fill lightness, and embedded
fonts; writes a manifest, an overlap audit, and grayscale / deuteranopia proof renders.  No TeX installation
is needed: the lettering uses Matplotlib's bundled Computer Modern fonts through mathtext.

Usage (from code_and_data):  python scripts/generate_bit_figures.py [--out DIR]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import logging
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("SOURCE_DATE_EPOCH", "1787184000")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib-cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib._mathtext_data  # noqa: E402
import matplotlib.colors  # noqa: E402
import matplotlib.patheffects as patheffects  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.text  # noqa: E402
import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import Collection, PathCollection  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.transforms import Bbox, ScaledTranslation, blended_transform_factory  # noqa: E402

# fontTools complains about the (deliberately ancient) timestamps of the bundled Computer Modern TrueType files
logging.getLogger("fontTools").setLevel(logging.ERROR)

# Computer Modern \varpi.  Matplotlib's BaKoMa ("cm") mathtext table has no entry for \varpi, so the glyph would
# fall back to STIXGeneral-Italic, which looks bold next to the Computer Modern omega and Gamma.  The bundled
# cmmi10.ttf does contain TeX's varpi (glyph "pi1" at code 36, the OML slot of \varpi), so it is mapped here.
# The entry lives only in this process's copy of Matplotlib's table (no file is changed), it is added only when
# Matplotlib has none of its own, and the BaKoMa table is consulted only by the "cm" mathtext fontset, which no
# other archive script uses.
CM_VARPI = ("cmmi10", 36)


def install_cm_varpi() -> bool:
    table = matplotlib._mathtext_data.latex_to_bakoma
    if r"\varpi" not in table:
        table[r"\varpi"] = CM_VARPI
    return table[r"\varpi"] == CM_VARPI


install_cm_varpi()

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

# ---------------------------------------------------------------------------
# house style
# ---------------------------------------------------------------------------
# Okabe--Ito palette (hue pairs kept for colour-vision safety).  Fills are opaque blends with white rather than
# alpha, because PostScript has no transparency and the EPS must render exactly like the PDF.
BLUE, ORANGE, GREEN, PURPLE, SKY, YELLOW, VERMIL = "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#D55E00"
INK, GUIDE = "#1A1A1A", "#9C9C9C"


def tint(color: str, amount: float) -> str:
    """Opaque blend of `color` with white (amount = 1 gives the pure colour)."""
    r, g, b = matplotlib.colors.to_rgb(color)
    return matplotlib.colors.to_hex((1 - amount * (1 - r), 1 - amount * (1 - g), 1 - amount * (1 - b)))


def shade(color: str, amount: float) -> str:
    """Opaque blend of `color` with black (amount = 1 gives the pure colour)."""
    r, g, b = matplotlib.colors.to_rgb(color)
    return matplotlib.colors.to_hex((amount * r, amount * g, amount * b))


FILL_BLUE, FILL_ORANGE = tint(BLUE, 0.34), tint(ORANGE, 0.30)              # attached / CPTP vs detached / non-CPTP
# stacked bars: the same full blue and orange as the lines elsewhere, and a light green, so that the three shown tones
# stay separated in grayscale (dark 0.43 / mid 0.68 / light 0.82 relative luminance, sRGB-encoded)
BAR_BLUE, BAR_ORANGE, BAR_GREEN, BAR_PURPLE = BLUE, ORANGE, tint(GREEN, 0.35), tint(PURPLE, 0.85)
MARK_PURPLE = shade(PURPLE, 0.6)                                            # darker purple for isolated point marks
DASHED = (0, (3.6, 1.8))
DASHDOT = (0, (4.2, 1.6, 1.0, 1.6))
# The bundled Computer Modern roman has no U+2013 glyph; its cmap keeps the TeX OT1 layout, in which the en dash
# sits at code point U+007B (verified by rendering: U+007B is the short dash, U+007C the em dash).
ENDASH = "{"
PANEL_PAD_PT = 3.0

EXPORT_RC = {
    "pdf.fonttype": 42,
    # Type 3 glyph programs survive EPS->PDF conversion in common pipelines; every glyph is embedded as vector outlines.
    "ps.fonttype": 3,
    # Computer Modern text and mathematics, matching the pdflatex article; tick labels go through mathtext so that the
    # minus sign, the multiplication sign and exponents come from the Computer Modern math fonts.
    "font.family": "serif",
    "font.serif": ["cmr10"],
    "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True,
    "axes.unicode_minus": True,
    "font.size": 8.5,
    "text.color": INK,
    "axes.labelsize": 8.5,
    "axes.titlesize": 8.5,
    "axes.labelcolor": INK,
    "axes.labelpad": 3.0,
    "axes.edgecolor": INK,
    "axes.linewidth": 0.55,
    "axes.grid": False,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 2.8,
    "ytick.major.size": 2.8,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.major.pad": 2.6,
    "ytick.major.pad": 2.6,
    "legend.fontsize": 8.0,
    "legend.frameon": False,
    "legend.handlelength": 1.8,
    "legend.handleheight": 0.7,
    "legend.handletextpad": 0.5,
    "legend.labelspacing": 0.3,
    "legend.columnspacing": 1.2,
    "legend.borderpad": 0.2,
    "legend.borderaxespad": 0.4,
    "lines.linewidth": 1.0,
    "lines.markersize": 4.0,
    "lines.markeredgewidth": 0.8,
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
    """Springer-style bold lowercase panel letter outside the frame at the top-left corner of the panel.

    Never an axes title (titles live in the captions).  finish() aligns it with the left extent of the panel's
    tick labels and axis label once the constrained layout is known."""
    return ax.text(0.0, 1.0, label, transform=ax.transAxes, fontfamily="cmb10", fontweight="bold", fontsize=9.5,
                   ha="left", va="bottom", gid="panel", zorder=5)


def finish(fig):
    """Run the layout and align every panel letter with the left edge of its panel's decorations."""
    labels = [(ax, t) for ax in fig.get_axes() for t in ax.texts if t.get_gid() == "panel"]
    for _, t in labels:
        t.set_visible(False)
    for _ in range(2):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        to_fig = fig.transFigure.inverted()
        for ax, t in labels:
            t.set_visible(False)
            left = ax.get_tightbbox(renderer).x0
            t.set_transform(blended_transform_factory(fig.transFigure, ax.transAxes)
                            + ScaledTranslation(0, PANEL_PAD_PT / 72, fig.dpi_scale_trans))
            t.set_position((to_fig.transform((left, 0.0))[0], 1.0))
            t.set_visible(True)
    return fig


def guide(ax, x=None, y=None, **kw):
    """Thin neutral reference line (critical frequency, marked step, tolerance)."""
    style = dict(color=GUIDE, linewidth=0.55, linestyle="-", zorder=1.5)
    style.update(kw)
    return ax.axvline(x, **style) if x is not None else ax.axhline(y, **style)


def mark(ax, x, y, text, **kw):
    """Small upright label with a white backing, for labels that sit on a guide line."""
    style = dict(fontsize=8.0, ha="center", va="center", zorder=4,
                 bbox=dict(boxstyle="square,pad=0.12", facecolor="white", edgecolor="none"))
    style.update(kw)
    return ax.text(x, y, text, **style)


def categorical(ax, axis):
    """Ticks only on the labelled side of a categorical axis."""
    ax.tick_params(axis=axis, which="both", **({"right": False} if axis == "y" else {"top": False}))


def fill_components(ax, vgrid, au, dl, du, labels=True):
    ax.fill_between(vgrid, 0, au, where=np.isfinite(au), facecolor=FILL_BLUE, edgecolor="none", linewidth=0,
                    label="attached component" if labels else None, zorder=1)
    ax.fill_between(vgrid, dl, du, where=np.isfinite(dl), facecolor=FILL_ORANGE, edgecolor="none", linewidth=0,
                    label="detached component" if labels else None, zorder=1)
    ax.plot(vgrid, au, color=BLUE, linestyle="-", linewidth=1.0, zorder=3)
    ax.plot(vgrid, dl, color=ORANGE, linestyle=DASHED, linewidth=1.0, zorder=3)
    ax.plot(vgrid, du, color=ORANGE, linestyle=DASHED, linewidth=1.0, zorder=3)


def component_handles():
    return [Patch(facecolor=FILL_BLUE, edgecolor=BLUE, linewidth=0.8, label="attached component"),
            Patch(facecolor=FILL_ORANGE, edgecolor=ORANGE, linewidth=0.8, linestyle=DASHED, label="detached component")]


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def fig1():
    rows = read_csv("figure1_exact_branches.csv")
    vgrid, au, dl, du = col(rows, "varpi"), col(rows, "attached_upper"), col(rows, "detached_lower"), col(rows, "detached_upper")
    fig, ax = plt.subplots(figsize=(WIDTH_IN, 3.2), layout="constrained")
    fill_components(ax, vgrid, au, dl, du)
    ceiling = ax.axhline(ALPHA4, color=INK, linestyle=DASHDOT, linewidth=0.8, label=r"population ceiling $\alpha_4$", zorder=2.5)
    ZOOM = (0.78, 0.88, 0.0, 3.0)  # x0, x1, y0, y1 of the inset window 0.78 <= varpi <= 0.88, 0 <= x <= 3 (outlined in the main panel)
    YTOP = 3.6
    # only varpi_+ is marked in the main panel; the three frequencies of the reentrant band are marked and labelled in the
    # inset.  The varpi_+ line stops at x = 3 (the top of the zoom box), below the legend, so it never strikes the legend text.
    guide(ax, x=VARPI_P, ymax=ZOOM[3] / YTOP)
    isolated, = ax.plot([VARPI_R, VARPI_P], [2.0, (9 + np.sqrt(33)) / ((123 + 11 * np.sqrt(33)) / 16)], marker="D", markersize=4.6,
                        markerfacecolor=YELLOW, markeredgecolor=INK, markeredgewidth=0.7, linestyle="None",
                        label="isolated admissible step", zorder=4)
    # \! pulls the script "+" in to the tight placement of the other varpi subscripts (mathtext sets it loosely)
    ax.text(VARPI_P + 0.05, 0.10, r"$\varpi_{\!+}$", fontsize=8.5, ha="left", va="bottom")
    ax.set_xlim(0, 4)
    ax.set_ylim(0, YTOP)
    ax.set_xlabel(r"frequency ratio $\varpi=|\omega|/\Gamma$")
    ax.set_ylabel(r"step $x=\Gamma h$")
    ax.legend(handles=component_handles() + [ceiling, isolated], loc="upper right", ncol=2, borderaxespad=0.6, columnspacing=1.6)
    # zoom into the reentrant band
    # (placed so that the labels above the inset, whose mathtext radical box is tall, stay clear of the ceiling line)
    inset = ax.inset_axes([0.625, 0.38, 0.355, 0.30])
    vz = np.linspace(0.78, 0.88, 1201)
    azu, zdl, zdu = _branches(vz)
    fill_components(inset, vz, azu, zdl, zdu, labels=False)
    # Labels sit above the inset frame, each centred on its line.  The varpi_- label ends its subscript with an
    # invisible text-mode space: Matplotlib's mathtext mis-sizes the raster of an expression that ends with a
    # script-size Computer Modern minus (that glyph has negative depth) and would crop the minus away.
    above = blended_transform_factory(inset.transData, inset.transAxes) + ScaledTranslation(0, 2.0 / 72, fig.dpi_scale_trans)
    for v, lab in ((VARPI_C, r"$\varpi_c$"), (VARPI_M, r"$\varpi_{\!-\text{ }}$"), (VARPI_R, r"$\sqrt{3}/2$")):
        guide(inset, x=v)
        inset.text(v, 1.0, lab, transform=above, fontsize=8.0, ha="center", va="bottom")
    # The data window is exactly ZOOM (and matches the grey box); the x ticks sit inside it so that no tick label
    # overhangs the inset frame (into the main panel's right spine) or meets the y tick label "0" at the corner.
    inset.set_xlim(ZOOM[0], ZOOM[1])
    inset.set_ylim(ZOOM[2], ZOOM[3])
    inset.set_xticks([0.80, 0.83, 0.86])
    inset.set_yticks([0, 1, 2, 3])
    inset.tick_params(labelsize=8.0, length=2.2, pad=2.0)
    inset.set_facecolor("white")
    ax.indicate_inset([ZOOM[0], ZOOM[2], ZOOM[1] - ZOOM[0], ZOOM[3] - ZOOM[2]], inset_ax=None, edgecolor=GUIDE, linestyle="-",
                      linewidth=0.6, alpha=1.0, zorder=2)
    with (RESULTS / "fig1_zoom_grid.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["varpi", "attached_upper", "detached_lower", "detached_upper"])
        w.writerows(zip(vz, azu, zdl, zdu))
    return finish(fig)


def conditioning_panel(bx):
    cond = read_csv("figure1_conditioning.csv")
    deltas, widths = col(cond, "delta_varpi"), col(cond, "window_width")
    ok = np.isfinite(widths)
    bx.loglog(deltas[ok], widths[ok], color=BLUE, marker="o", markevery=12, linestyle="-", label=r"$x_{\!+}-x_{\!-}$", zorder=3)
    # square-root reference, offset upwards by a factor 3 from the last data point so that it does not hide under the data
    ref = 3.0 * np.sqrt(deltas[ok]) * widths[ok][-1] / np.sqrt(deltas[ok][-1])
    bx.loglog(deltas[ok], ref, color=ORANGE, linestyle=DASHED, label=r"$O((\varpi-\varpi_{\!+})^{1/2})$", zorder=4)
    bx.set_xlabel(r"$\varpi-\varpi_{\!+}$")
    bx.set_ylabel("detached-window width")
    bx.legend(loc="lower right")


def fig2():
    V = load_json("results/v34/v34_verification.json")
    direct = V["candidate_regions_varpi_2"]["full"][0]
    fine = V["candidate_regions_varpi_2"]["two_half"][0]
    marg = read_csv("figure2_richardson_margin.csv")
    xs, scaled = col(marg, "x"), col(marg, "M_ext_over_x6")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH_IN, 2.5), layout="constrained")
    panel(ax, "a")
    # At varpi = 2 both positive-step components are detached (A_4(2) = {0} u [x_-, x_+]), so both use the detached
    # style of Figs. 1 and 3 (orange, square markers); the two candidates differ by marker fill only (filled: one full
    # step, open: two half steps), and the rows are named on the axis, so no colour name is needed.
    rows = [(1, direct["lower_decimal"], direct["upper_decimal"], ORANGE, "s", ORANGE), (0, fine["lower_decimal"], fine["upper_decimal"], ORANGE, "s", "white")]
    for y, lo, hi, color, mk, face in rows:
        ax.hlines(y, lo, hi, linewidth=3.0, color=color, linestyles="-", zorder=2)
        # squares slightly larger than the 3 pt bar so that the filled ones still read as markers, not as bar ends
        ax.plot([lo, hi], [y, y], color=color, marker=mk, markerfacecolor=face, markeredgewidth=0.8, markersize=5.0, linestyle="None", zorder=3)
    ax.plot([0, 0], [1, 0], marker="o", markersize=3.0, color=INK, linestyle="None", zorder=3)
    for xv in (1.0, 2.0):
        guide(ax, x=xv)
        mark(ax, xv, 1.72, fr"$x={int(xv)}$", fontsize=8.5)
    ax.set_yticks([1, 0])
    ax.set_yticklabels(["one full step", "two half steps"])
    categorical(ax, "y")
    ax.set_xlim(-0.05, 2.75)
    ax.set_ylim(-0.6, 2.05)
    ax.set_xlabel(r"total step $x=\Gamma h$ at $\varpi=2$")
    panel(bx, "b")
    margin, = bx.plot(xs, scaled, color=PURPLE, linestyle="-", label=r"$M_{\rm ext}(x)/x^6$", zorder=3)
    bx.axhline(0, color=INK, linewidth=0.6, zorder=2)
    bx.fill_between(xs, scaled, 0, where=scaled < 0, facecolor=FILL_ORANGE, edgecolor="none", linewidth=0,
                    label="extrapolate non-CPTP", zorder=1)
    region = Patch(facecolor=FILL_ORANGE, edgecolor=ORANGE, linewidth=0.8, label="extrapolate non-CPTP")
    guide(bx, x=ALPHA4)
    bx.text(ALPHA4 - 0.06, -2.45e-4, r"$\alpha_4$", ha="right", va="bottom", fontsize=8.5)
    # clear space between the alpha_4 line and the right spine (a 1 mm gap read as a doubled frame line)
    bx.set_xlim(0, 3.0)
    bx.set_ylim(-2.6e-4, 1.4e-4)
    bx.set_xlabel(r"nonrotating step $x$")
    bx.set_ylabel("scaled Choi margin")
    bx.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)
    # upper left, in the empty band above zero and clear of the alpha_4 line
    bx.legend(handles=[margin, region], loc="upper left", borderaxespad=0.6, handlelength=1.5, handletextpad=0.5)
    return finish(fig)


def fig3():
    V = load_json("results/v27/v27_verification.json")
    buffered = V["buffered_components"]
    nc = read_csv("figure3_noncommuting_endpoints.csv")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH_IN, 2.6), layout="constrained", gridspec_kw={"width_ratios": [1.1, 1.0]})
    panel(ax, "a")
    ypos = np.arange(len(buffered))[::-1]
    labels = []
    for y, row in zip(ypos, buffered):
        th, ka = row["theta"], row["kappa"]
        labels.append(fr"$\theta={th}$, $\kappa={ka}$".replace("0.001", "10^{-3}").replace("0.01", "10^{-2}"))
        for j, (lo, hi) in enumerate(row["components"]):
            color, mk = (BLUE, "o") if j == 0 else (ORANGE, "s")
            ax.hlines(y, float(lo), float(hi), linewidth=2.6, color=color, linestyles="-", zorder=2)
            ax.plot([float(lo), float(hi)], [y, y], marker=mk, linestyle="None", markersize=4.0, color=color,
                    markerfacecolor="white" if j else color, markeredgewidth=0.8, zorder=3)
    ax.plot([], [], linewidth=2.6, color=BLUE, linestyle="-", marker="o", markersize=4.0, label="attached")
    ax.plot([], [], linewidth=2.6, color=ORANGE, linestyle="-", marker="s", markersize=4.0, markerfacecolor="white", label="detached")
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    categorical(ax, "y")
    ax.set_xlim(-0.02, 1.36)
    ax.set_ylim(-0.7, 6.0)
    ax.set_xlabel(r"$x=\Gamma h$ at $\varpi=2$")
    ax.legend(loc="upper center", ncol=2, handlelength=1.6, columnspacing=1.2, handletextpad=0.5)
    panel(bx, "b")
    ox = col(nc, "omega_x")
    bx.plot(ox, col(nc, "attached_upper"), color=BLUE, marker="o", linestyle="None", markersize=4.6, label="attached, upper")
    bx.plot(ox, col(nc, "detached_lower"), color=ORANGE, marker="s", markerfacecolor="white", markeredgewidth=0.9, linestyle="None", markersize=4.6, label="detached, lower")
    bx.plot(ox, col(nc, "detached_upper"), color=GREEN, marker="^", markerfacecolor="white", markeredgewidth=0.9, linestyle="None", markersize=5.2, label="detached, upper")
    bx.set_xticks(ox)
    bx.set_xticklabels(["$0$", "$0.05$", "$0.10$", "$0.20$"])
    bx.set_xlim(-0.03, 0.235)
    bx.set_ylim(0.15, 1.40)
    bx.set_xlabel(r"coupling $\Omega_x/\Gamma$")
    bx.set_ylabel(r"endpoint $x$")
    bx.legend(loc="center", bbox_to_anchor=(0.5, 0.28), handletextpad=0.4, handlelength=1.2)
    return finish(fig)


def fig4():
    sweep = read_csv("figure5_tolerance_sweep.csv")
    comp = read_csv("figure5_action_composition.csv")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH_IN, 2.55), layout="constrained", gridspec_kw={"width_ratios": [1.0, 1.2]})
    panel(ax, "a")
    styles = {"error_only": (BLUE, "o", "-", BLUE, "error only"), "guard_rotating_fallback": (ORANGE, "s", DASHED, "white", "guard + fallback"),
              "rotating_frame": (GREEN, "^", DASHDOT, "white", "rotating frame")}
    for policy, (color, mk, ls, face, label) in styles.items():
        data = sorted([r for r in sweep if r["policy"] == policy], key=lambda r: float(r["tolerance"]), reverse=True)
        work = np.array([float(r["rhs_stage_evaluations"]) + 3 * float(r["exponential_actions"]) for r in data])
        err = np.array([float(r["global_normalized_choi_error"]) for r in data])
        ax.loglog(work, err, color=color, marker=mk, markerfacecolor=face, markeredgewidth=0.8, linestyle=ls, label=label)
    ax.set_xlabel("propagator work")
    ax.set_ylabel("global normalized Choi error")
    ax.set_ylim(1e-6, 2e-1)
    ax.set_xlim(22, 260)
    ax.set_xticks([30, 60, 100, 200])
    ax.set_xticklabels(["30", "60", "100", "200"])
    ax.set_xticks([], minor=True)
    ax.legend(loc="upper left", handlelength=2.0)
    panel(bx, "b")
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
    for heights, bottoms, label, color in ((direct, 0 * direct, "direct", BAR_BLUE), (projected, direct, "projected", BAR_ORANGE),
                                           (rotating, direct + projected, "rotating fallback", BAR_GREEN),
                                           (exact, direct + projected + rotating, "exact fallback", BAR_PURPLE)):
        if np.any(heights):
            bx.bar(pos, heights, bottom=bottoms, label=label, color=color, edgecolor=INK, linewidth=0.5, width=0.6)
    bx.set_xticks(pos)
    bx.set_xticklabels(names)
    categorical(bx, "x")
    bx.set_ylabel("accepted candidates")
    bx.set_ylim(0, 8.2)
    bx.legend(loc="upper right", ncol=1, handlelength=1.4, handletextpad=0.5)
    return finish(fig)


def figS1():
    fig, ax = plt.subplots(figsize=(WIDTH_IN, 2.0), layout="constrained")
    vmax, vlo, vr, vhi = 2.3, np.sqrt(2 - np.sqrt(33) / 4), VARPI_R, np.sqrt(2 + np.sqrt(33) / 4)
    # method names as in the article's Table 1
    # the Table 1 name in full, on two lines so that the bars keep their width
    rows = [("Richardson extrapolate\n" r"of RK4 ($n=2$)", [(0, vlo, "out"), (vlo, vhi, "in"), (vhi, vmax, "out")]),
            (f"Dormand{ENDASH}Prince 5 (principal)", [(0, vlo, "in"), (vlo, vhi, "out"), (vhi, vmax, "in")]),
            ("Classical RK4", [(0, vr, "in"), (vr, vmax, "out")])]
    for y, (name, segs) in enumerate(rows):
        for lo, hi, state in segs:
            ax.barh(y, hi - lo, left=lo, height=0.56, facecolor=FILL_BLUE if state == "in" else FILL_ORANGE,
                    edgecolor=BLUE if state == "in" else ORANGE, linewidth=0.6, zorder=2)
            ax.text((lo + hi) / 2, y, "CPTP" if state == "in" else "non-CPTP", ha="center", va="center", fontsize=8.0, gid="inbar", zorder=4)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows])
    categorical(ax, "y")
    for v in (vlo, vr, vhi):
        guide(ax, x=v, zorder=3)
    # rounded (not truncated) three-decimal tick-style labels, each beside its line: sqrt(2 - sqrt33/4) = 0.7509...,
    # sqrt3/2 = 0.8660..., sqrt(2 + sqrt33/4) = 1.8536...
    ax.text(vlo - 0.025, 2.58, f"{vlo:.3f}", ha="right", va="bottom", fontsize=8.0)
    ax.text(vr + 0.025, 2.58, f"{vr:.3f}", ha="left", va="bottom", fontsize=8.0)
    ax.text(vhi + 0.025, 2.58, f"{vhi:.3f}", ha="left", va="bottom", fontsize=8.0)
    ax.set_xlim(0, vmax)
    ax.set_ylim(-0.5, 3.1)
    ax.set_xlabel(r"frequency ratio $\varpi=|\omega|/\Gamma$")
    return finish(fig)


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
    panel(cx, "a")
    conditioning_panel(cx)
    panel(ax, "b")
    xs = np.linspace(0.02, 1.45, 900)
    q = 17.0
    signed = -((q * xs) ** 3 - 16 * (q * xs) ** 2 - 32 * (q - 6) * (q * xs) + 384 * (q - 4))
    curve, = ax.plot(xs, signed, color=BLUE, linestyle="-", label=r"$-\mathcal{C}_{17}(17x)$", zorder=3)
    ax.axhline(0, color=INK, linewidth=0.6, zorder=2)
    ax.axvspan(xminus, xplus, facecolor=FILL_ORANGE, edgecolor="none", linewidth=0, zorder=1)
    for xv in (xminus, xplus):
        ax.axvline(xv, color=ORANGE, linewidth=0.8, linestyle=DASHED, zorder=2)
    window = Patch(facecolor=FILL_ORANGE, edgecolor=ORANGE, linewidth=0.8, linestyle=DASHED, label="admissible")
    seq = np.array([1.4 / 2**j for j in range(6)])
    seq_y = -((q * seq) ** 3 - 16 * (q * seq) ** 2 - 32 * (q - 6) * (q * seq) + 384 * (q - 4))
    halving, = ax.plot(seq, seq_y, color=PURPLE, marker="o", markerfacecolor="white", markeredgewidth=0.9, linestyle="None", markersize=4.6,
                       label="halving steps", zorder=4)
    ax.set_yscale("symlog", linthresh=100)
    ax.set_xlabel(r"$x=\Gamma h$ at $\varpi=2$")
    ax.set_ylabel("signed polynomial\n(symmetric log scale)")
    ax.legend(handles=[curve, window, halving], loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=1, borderaxespad=0.4, handlelength=1.6, handletextpad=0.5)
    panel(bx, "c")
    direct = np.abs(dbl)
    direct[direct == 0] = np.nan
    bx.loglog(x, np.abs(exact), color=BLUE, marker="o", linestyle="-", label=r"exact $|M_4|$", zorder=3)
    bx.loglog(x, direct, color=ORANGE, marker="s", markerfacecolor="white", markeredgewidth=0.8, linestyle=DASHED, label="double precision", zorder=3)
    guide(bx, y=1e-12, label=r"tolerance $10^{-12}$")
    false = np.array([0.004, 0.001, 1e-4, 1e-5])
    false_y = [abs(float(r["exact_margin_decimal"])) for r in V["blocking_false_acceptance_reproduction"]]
    # heavier, darker crosses with a thin white halo, so that they stay legible on the blue line in grayscale print
    bx.scatter(false, false_y, color=MARK_PURPLE, marker="x", s=36, linewidths=1.5, label="false accept", zorder=4,
               path_effects=[patheffects.withStroke(linewidth=2.9, foreground="white")])
    bx.set_xlabel(r"$x$")
    bx.set_ylabel("margin magnitude")
    bx.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=1, borderaxespad=0.4, handlelength=1.6, handletextpad=0.5)
    return finish(fig)


def figS3():
    rows = read_csv("figureS2_operational.csv")
    vscan, ratios, absdef = col(rows, "varpi"), col(rows, "ratio"), col(rows, "absolute_defect")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH_IN, 2.8), layout="constrained")
    panel(ax, "a")
    l1 = ax.plot(vscan, ratios, color=BLUE, linestyle="-", label=r"$\delta_{\rm CP}/E_J$")
    ax.set_xlabel(r"$\varpi$ at $x=10^{-2}$")
    ax.set_ylabel("fraction of Bell-input trace error")
    ax.set_ylim(0, 1)
    ax2 = ax.twinx()
    l2 = ax2.semilogy(vscan, absdef, color=ORANGE, linestyle=DASHED, label="absolute negative mass")  # same term as the right axis
    ax2.set_ylabel("absolute negative mass")
    ax.legend(l1 + l2, [l.get_label() for l in l1 + l2], loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=1, borderaxespad=0.4, handlelength=1.8, handletextpad=0.5)
    panel(bx, "b")
    mp.mp.dps = 60
    Ns = np.unique(np.logspace(1, 5, 60).astype(int))
    for name, fun, color, mk, ls, face in (("Forward Euler", lambda z: 1 + z, BLUE, "o", "-", BLUE),
                                           ("SSPRK(3,3)", lambda z: 1 + z + z**2 / 2 + z**3 / 6, ORANGE, "s", DASHED, "white")):
        vals = []
        for N in Ns:
            h = mp.mpf(1) / int(N)
            A = fun(-h) ** int(N)
            C = fun(-h / 2) ** int(N)
            lam = (1 + A - mp.sqrt((1 - A) ** 2 + 4 * abs(C) ** 2)) / 4
            vals.append(float(-lam))
        bx.loglog(1 / Ns, vals, color=color, marker=mk, markevery=8, markersize=3.6, markerfacecolor=face, markeredgewidth=0.8, linestyle=ls, label=name)
    bx.set_xlabel(r"step $x_N=1/N$")
    bx.set_ylabel(r"$-p^{\rm RK}_{-,N}$")
    bx.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=1, borderaxespad=0.4, handlelength=1.8, handletextpad=0.5)
    return finish(fig)


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
def _densify(xy, step=2.0):
    """Points every `step` display pixels along a display-space polyline, so that a legend or label crossed by a
    long straight segment (e.g. a full-height guide line, whose only vertices lie on the frame) is detected."""
    out = [xy[:1]]
    for a, b in zip(xy[:-1], xy[1:]):
        n = max(1, int(np.ceil(np.hypot(*(b - a)) / step)))
        out.append(a + np.linspace(0.0, 1.0, n + 1)[1:, None] * (b - a))
    return np.vstack(out)


def _runs(xy):
    """Split a vertex array at non-finite rows (line breaks) into finite runs."""
    ok = np.isfinite(xy).all(axis=1)
    runs, start = [], None
    for i, flag in enumerate(ok):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            runs.append(xy[start:i])
            start = None
    if start is not None:
        runs.append(xy[start:])
    return runs


def _data_points(ax, dense=False):
    """Display-space points and rectangles of the data artists of an axes (interior points only, so that
    axvline/axhline endpoints on the frame do not count).  With dense=True, drawn line segments and collection
    edges are sampled every 2 px instead of only at their vertices; marker-only lines and scatter markers
    contribute their centres."""
    pts = []
    rects = []
    for line in ax.lines:
        xy = np.column_stack([np.asarray(line.get_xdata(), dtype=float), np.asarray(line.get_ydata(), dtype=float)])
        drawn = line.get_linestyle() not in ("None", "none", "", " ")
        for run in _runs(xy):
            d = line.get_transform().transform(run)
            pts.append(_densify(d) if dense and drawn and len(d) > 1 else d)
    for coll in ax.collections:
        if not isinstance(coll, Collection):
            continue
        if isinstance(coll, PathCollection):          # scatter: marker centres are the data positions
            off = np.asarray(coll.get_offsets(), dtype=float)
            off = off[np.isfinite(off).all(axis=1)] if len(off) else off
            if len(off):
                pts.append(coll.get_offset_transform().transform(off))
            continue
        tr = coll.get_transform()
        for path in coll.get_paths():
            polys = path.to_polygons(closed_only=False) if dense else [np.asarray(path.vertices, dtype=float)]
            for poly in polys:
                for run in _runs(np.asarray(poly, dtype=float)):
                    d = tr.transform(run)
                    pts.append(_densify(d) if dense and len(d) > 1 else d)
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
    with plt.rc_context(EXPORT_RC):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        report = {"legend_data_overlaps": [], "text_data_overlaps": [], "text_text_overlaps": [], "min_font_pt": None}
        boxes = []
        for ax in fig.get_axes():
            pts, rects = _data_points(ax)
            dense, _ = _data_points(ax, dense=True)
            leg = ax.get_legend()
            if leg is not None:
                lb = leg.get_window_extent(renderer)
                boxes.append(("legend@" + repr(ax.get_position().bounds), lb))
                # frameless legends must not be crossed by any drawn line, guide lines included
                if _bbox_hits(lb, dense, rects):
                    report["legend_data_overlaps"].append(f"legend of axes {ax.get_position().bounds}")
                ab = ax.get_window_extent(renderer)
                inside = lb.x0 >= ab.x0 - 0.5 and lb.x1 <= ab.x1 + 0.5 and lb.y0 >= ab.y0 - 0.5 and lb.y1 <= ab.y1 + 0.5
                outside = lb.x1 <= ab.x0 + 0.5 or lb.x0 >= ab.x1 - 0.5 or lb.y1 <= ab.y0 + 0.5 or lb.y0 >= ab.y1 - 0.5
                if not inside and not outside:
                    report["legend_data_overlaps"].append(f"legend crosses the axes frame of axes {ax.get_position().bounds}")
            # axes titles are never used (panel letters are text artists), but keep them in the text-text test
            for t in (ax.title, ax._left_title, ax._right_title):
                if t.get_text().strip():
                    boxes.append(("title:" + t.get_text(), t.get_window_extent(renderer)))
            # axis labels and offset texts (e.g. a long y label running past the axes into the panel letter)
            for t in (ax.xaxis.label, ax.yaxis.label, ax.xaxis.offsetText, ax.yaxis.offsetText):
                if t.get_visible() and t.get_text().strip():
                    boxes.append(("axis:" + t.get_text(), t.get_window_extent(renderer)))
            # every text artist, including the panel letters, must clear the data and each other
            for t in ax.texts:
                if not t.get_text().strip() or t.get_gid() == "inbar":
                    continue
                tb = t.get_window_extent(renderer)
                boxes.append((t.get_text(), tb))
                # a label with a white backing (mark()) is meant to sit on a thin guide line, so only the data
                # vertices count for it; every other label must also clear the drawn line segments
                probe = pts if t.get_bbox_patch() is not None else dense
                if _bbox_hits(tb, probe, rects, pad=0.5):
                    report["text_data_overlaps"].append(t.get_text())
            # child inset axes as opaque boxes against the parent's data; their labels against both panels
            for child in ax.child_axes:
                cb = child.get_window_extent(renderer)
                if _bbox_hits(cb, dense, rects, pad=0.0):
                    report["text_data_overlaps"].append("inset axes")
                cpts, crects = _data_points(child)
                for t in child.texts:
                    if not t.get_text().strip():
                        continue
                    tb = t.get_window_extent(renderer)
                    boxes.append(("inset:" + t.get_text(), tb))
                    if _bbox_hits(tb, dense, rects, pad=0.5) or _bbox_hits(tb, cpts, crects, pad=0.5):
                        report["text_data_overlaps"].append("inset:" + t.get_text())
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
    with plt.rc_context(EXPORT_RC):
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
