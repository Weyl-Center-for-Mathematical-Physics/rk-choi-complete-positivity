#!/usr/bin/env python3
from __future__ import annotations

import csv
import datetime as dt
import json
import os
from pathlib import Path

os.environ.setdefault("SOURCE_DATE_EPOCH", "1786320000")

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import sympy as sp

from rk_choi_margin.topology import normalized_choi_defect_and_error

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures_jcp"
DATA = ROOT / "results" / "jcp_figures_v27"
VERIFY = ROOT / "results" / "v27" / "v27_verification.json"
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)
V = json.loads(VERIFY.read_text(encoding="utf-8"))

FIXED_DATE = dt.datetime(2026, 8, 10, tzinfo=dt.timezone.utc)
PDF_META = {
    "Title": "Certified Runge-Kutta complete-positivity figures, v2.7",
    "Author": "G. Blake Pierpoint, Olivier Bernard, Yichen Liu",
    "Creator": "generate_jcp_figures_v27.py",
    "CreationDate": FIXED_DATE,
    "ModDate": FIXED_DATE,
}
PNG_META = {"Software": "generate_jcp_figures_v27.py", "Creation Time": "2026-08-10T00:00:00Z"}


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight", metadata=PDF_META)
    fig.savefig(FIG / f"{stem}.png", dpi=500, bbox_inches="tight", metadata=PNG_META)
    fig.savefig(FIG / f"{stem}.eps", bbox_inches="tight")
    plt.close(fig)


def _plotting_roots(varpi: float) -> list[float]:
    """Companion roots for dense rendering only; no scientific decision uses them."""
    q = 1 + 4 * varpi**2
    roots = np.roots([1.0, -16.0, -32.0 * (q - 6.0), 384.0 * (q - 4.0)])
    return sorted(float(r.real / q) for r in roots if abs(float(r.imag)) < 1e-7 and r.real > 1e-10)


def rk4_intervals_plot(varpi: float, alpha4: float = 2.785293563405282):
    qm = (123 - 11 * np.sqrt(33)) / 16
    qp = (123 + 11 * np.sqrt(33)) / 16
    q = 1 + 4 * varpi**2
    roots = _plotting_roots(varpi)
    eps = 2e-7
    if q < qm - eps:
        return ([(0.0, min(alpha4, roots[-1]))] if roots else []), []
    if abs(q - qm) <= eps:
        upper = min(alpha4, (2 * np.sqrt(33) - 2) / qm)
        return [(0.0, upper)], [(9 - np.sqrt(33)) / qm]
    if q < 4 - eps:
        return ([(0.0, roots[0]), (roots[1], min(alpha4, roots[2]))] if len(roots) == 3 else []), []
    if abs(q - 4) <= eps:
        return [], [0.0, 2.0]
    if q < qp - eps:
        return [], [0.0]
    if abs(q - qp) <= eps:
        return [], [0.0, (9 + np.sqrt(33)) / qp]
    return ([(roots[0], min(alpha4, roots[1]))] if len(roots) == 2 else []), [0.0]


def exact_phase_fill(ax, vgrid):
    au = np.full_like(vgrid, np.nan)
    dl = np.full_like(vgrid, np.nan)
    du = np.full_like(vgrid, np.nan)
    for i, vv in enumerate(vgrid):
        intervals, _points = rk4_intervals_plot(float(vv))
        if intervals:
            if intervals[0][0] == 0:
                au[i] = intervals[0][1]
                if len(intervals) > 1:
                    dl[i], du[i] = intervals[1]
            else:
                dl[i], du[i] = intervals[0]
    ax.fill_between(vgrid, 0, au, where=np.isfinite(au), alpha=0.28, label="attached component")
    ax.fill_between(vgrid, dl, du, where=np.isfinite(dl), alpha=0.43, label="detached component")
    ax.plot(vgrid, au, linewidth=1.15)
    ax.plot(vgrid, dl, linewidth=1.15)
    ax.plot(vgrid, du, linewidth=1.15)
    return au, dl, du


def fig1_topology_conditioning():
    vm = float(sp.N(sp.sqrt(107 - 11 * sp.sqrt(33)) / 8, 17))
    vu = np.sqrt(3) / 2
    vp = float(sp.N(sp.sqrt(107 + 11 * sp.sqrt(33)) / 8, 17))
    vgrid = np.linspace(0, 4, 1801)
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.5), gridspec_kw={"width_ratios": [1.62, 1]})
    ax = axes[0]
    au, dl, du = exact_phase_fill(ax, vgrid)
    ax.axhline(2.785293563405282, ls="--", lw=0.9, label=r"population ceiling $\alpha_4$")
    for value, label in [(vm, r"$\varpi_-$"), (vu, r"$\sqrt{3}/2$"), (vp, r"$\varpi_+$")]:
        ax.axvline(value, ls=":", lw=0.9)
        ax.text(value + 0.025, 0.12, label, rotation=90, va="bottom", fontsize=8)
    ax.set(xlim=(0, 4), ylim=(0, 3), xlabel=r"frequency ratio $\varpi=|\omega|/\Gamma$", ylabel=r"step $x=\Gamma h$", title="(a) Exact RK4 channel-admissibility topology")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    # Inset is placed above the lower-left boundary and away from varpi_+.
    inset = ax.inset_axes([0.075, 0.50, 0.42, 0.36])
    vz = np.linspace(0.80, 0.88, 1101)
    exact_phase_fill(inset, vz)
    inset.axvline(vm, ls=":", lw=0.7); inset.axvline(vu, ls=":", lw=0.7)
    inset.set(xlim=(0.80, 0.88), ylim=(0, 2.9), xticks=[0.80, 0.83, 0.86, 0.88], yticks=[0, 1, 2])
    inset.tick_params(labelsize=7); inset.set_title("narrow reentrant band", fontsize=8)

    ax = axes[1]
    deltas = np.logspace(-12, -2, 120)
    separations = []
    for d in deltas:
        roots = _plotting_roots(vp + d)
        separations.append(roots[1] - roots[0] if len(roots) == 2 else np.nan)
    ax.loglog(deltas, separations, label=r"$x_+-x_-$")
    ax.loglog(deltas, np.sqrt(deltas) * separations[-1] / np.sqrt(deltas[-1]), ls="--", label=r"$O((\varpi-\varpi_+)^{1/2})$")
    ax.set(xlabel=r"distance above saddle node $\varpi-\varpi_+$", ylabel="detached-window width", title="(b) Double-root conditioning")
    ax.grid(which="both", alpha=0.18); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); save(fig, "Fig1_RK4_phase_diagram")
    with (DATA / "figure1_exact_branches.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["varpi", "attached_upper", "detached_lower", "detached_upper"]); w.writerows(zip(vgrid, au, dl, du))
    with (DATA / "figure1_conditioning.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["delta_varpi", "window_width"]); w.writerows(zip(deltas, separations))


def fig2_robustness_noncommuting():
    buffered = V["buffered_components"]
    noncommuting = V["noncommuting_transverse_field"]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5))
    ax = axes[0]
    labels = []
    for row in buffered:
        labels.append(fr"$\theta={row['theta']},\ \kappa={row['kappa']}$")
    ypos = np.arange(len(buffered))[::-1]
    for y, row in zip(ypos, buffered):
        for lo, hi in row["components"]:
            ax.hlines(y, lo, hi, lw=5); ax.plot([lo, hi], [y, y], "o", ms=3.5)
    ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=8)
    ax.set(xlim=(-0.02, 1.36), xlabel=r"$x=\Gamma h$ at $\varpi=2$", title="(a) Buffered phase-covariant cases")
    ax.grid(axis="x", alpha=0.18); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    ax = axes[1]
    ox = np.array([float(row["omega_x"]) for row in noncommuting])
    endpoints = []
    for row in noncommuting:
        comps = [[float(v) for v in pair] for pair in row["components"]]
        endpoints.append([comps[0][1], comps[1][0], comps[1][1]])
    endpoints = np.asarray(endpoints)
    ax.plot(ox, endpoints[:, 0], "o-", label="attached endpoint")
    ax.plot(ox, endpoints[:, 1], "s-", label="detached lower endpoint")
    ax.plot(ox, endpoints[:, 2], "^-", label="detached upper endpoint")
    ax.fill_between(ox, 0, endpoints[:, 0], alpha=0.18)
    ax.fill_between(ox, endpoints[:, 1], endpoints[:, 2], alpha=0.25)
    ax.set(xlabel=r"transverse coupling $\Omega_x/\Gamma$", ylabel=r"certified step boundary $x$", title="(b) Persistence under noncommuting Hamiltonian perturbation")
    ax.grid(alpha=0.18); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); save(fig, "Fig2_robustness")
    with (DATA / "figure2_noncommuting_endpoints.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["omega_x", "attached_upper", "detached_lower", "detached_upper"])
        for v, row in zip(ox, endpoints): w.writerow([v, *row])


def fig3_certified_guard_conditioning():
    conditioning = V["floating_point_conditioning"]
    rows = sorted(conditioning, key=lambda row: float(sp.Rational(row["x"])))
    x = np.array([float(sp.Rational(row["x"])) for row in rows])
    exact = np.array([float(row["exact_margin"]) for row in rows])
    dbl = np.array([float(row["binary64_margin"]) for row in rows])
    scaled = np.array([float(row["scaled_sign_polynomial_minus_C"]) for row in rows])
    lo = [float(sp.Rational(v)) for v in V["varpi2_detached_interval"]["x_minus_bracket"]]
    hi = [float(sp.Rational(v)) for v in V["varpi2_detached_interval"]["x_plus_bracket"]]
    xminus = sum(lo) / 2
    xplus = sum(hi) / 2

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.3))
    ax = axes[0]
    xs = np.linspace(0.02, 1.45, 900)
    q = 17.0
    C = (q * xs) ** 3 - 16 * (q * xs) ** 2 - 32 * (q - 6) * (q * xs) + 384 * (q - 4)
    signed = -C
    ax.plot(xs, signed, label=r"certified sign function $-\mathcal{C}_{17}(17x)$")
    ax.axhline(0, lw=0.8)
    ax.axvspan(xminus, xplus, alpha=0.24, label="CPTP component")
    seq = np.array([1.4 / 2**j for j in range(6)])
    seq_y = -((q * seq) ** 3 - 16 * (q * seq) ** 2 - 32 * (q - 6) * (q * seq) + 384 * (q - 4))
    ax.plot(seq, seq_y, "o", ms=4, label="halving proposals")
    for j, (xx, yy) in enumerate(zip(seq, seq_y)):
        ax.annotate(str(j), (xx, yy), xytext=(3, 4), textcoords="offset points", fontsize=7)
    ax.set_yscale("symlog", linthresh=100)
    ax.set(xlabel=r"$x=\Gamma h$ at $\varpi=2$", ylabel="signed polynomial (symlog)", title="(a) Geometry and step-reduction skip")
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    direct = np.abs(dbl)
    direct[direct == 0] = np.nan
    ax.loglog(x, np.abs(exact), "o-", label="exact $|M_4|$")
    ax.loglog(x, direct, "s--", label="binary64 direct subtraction")
    ax.axhline(1e-12, ls=":", lw=0.9, label=r"legacy absolute tolerance $10^{-12}$")
    false = np.array([0.004, 0.001, 1e-4, 1e-5])
    false_y = [abs(float(row["exact_margin_decimal"])) for row in V["blocking_false_acceptance_reproduction"]]
    ax.scatter(false, false_y, marker="x", s=45, label="legacy false accepts")
    ax.set(xlabel=r"$x$", ylabel="margin magnitude", title="(b) Cancellation defeats a fixed absolute tolerance")
    ax.grid(which="both", alpha=0.18)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, "Fig3_certified_guard")
    with (DATA / "figure3_conditioning.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["x", "exact_margin", "binary64_margin", "scaled_sign"])
        w.writerows(zip(x, exact, dbl, scaled))

def fig4_adaptive_rotating_frame():
    rows = [r for r in V["adaptive_benchmarks"] if r["scenario"] == "boundary"]
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 4.4))
    ax = axes[0]
    markers = {"error_only": "o", "certified_guard": "s", "rotating_frame": "^"}
    for row in rows:
        accepted = [s for s in row["steps"] if s["used_h"] > 0]
        cumulative = np.cumsum([s["used_h"] for s in accepted])
        ax.step(np.r_[0, cumulative], np.r_[[accepted[0]["used_h"]], [s["used_h"] for s in accepted]], where="pre", label=row["policy"].replace("_", " "))
        ax.plot(cumulative, [s["used_h"] for s in accepted], markers[row["policy"]], ms=4)
    ax.set(xlabel="integration time", ylabel="accepted step", title=r"(a) Generator-level PI benchmark, $\varpi=2$")
    ax.grid(alpha=0.18); ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    labels = [r["policy"].replace("_", "\n") for r in rows]
    errors = [r["global_normalized_choi_error"] for r in rows]
    defects = [max(0.0, -r["minimum_step_choi_eigenvalue"]) for r in rows]
    pos = np.arange(len(rows))
    ax.bar(pos - 0.16, errors, 0.32, label="global normalized-Choi error")
    # Only plot a CP-defect bar when it is scientifically resolved above roundoff.
    for i, defect in enumerate(defects):
        if defect > 1e-12:
            ax.bar(i + 0.16, defect, 0.32, label="largest step CP defect" if i == 0 else None)
        else:
            ax.annotate("certified CP", (i + 0.16, 2e-6), ha="center", va="bottom", rotation=90, fontsize=7)
    ax.set_yscale("log")
    ax.set_ylim(1e-6, 4e-3)
    ax.set_xticks(pos)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set(ylabel="normalized-Choi error / CP defect", title="(b) Accuracy and physicality")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", which="both", alpha=0.18)
    fig.tight_layout(); save(fig, "Fig4_adaptive_frame")
    with (DATA / "figure4_adaptive_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["policy","global_normalized_choi_error","max_cp_defect","accepted_steps","rhs_stage_evaluations","fallback_steps"])
        for r in rows: w.writerow([r["policy"],r["global_normalized_choi_error"],max(0,-r["minimum_step_choi_eigenvalue"]),r["stats"]["accepted_steps"],r["stats"]["rhs_stage_evaluations"],r["stats"]["fallback_steps"]])


def fig5_operational_refinement():
    mp.mp.dps = 80
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.25))
    ax=axes[0]
    vscan=np.linspace(.9,10,460); ratios=[]; absdef=[]
    for vv in vscan:
        d,e=normalized_choi_defect_and_error(1e-2,float(vv)); ratios.append(float(d/e) if e else np.nan); absdef.append(float(d))
    ax.plot(vscan,ratios,label=r"$\delta_{\rm CP}/E_J$")
    ax.set(xlabel=r"$\varpi$ at $x=10^{-2}$",ylabel="fraction of Bell-input trace error",title="(a) Physicality defect is an order-one fraction, not a channel norm")
    ax.set_ylim(0,1); ax.grid(alpha=.18)
    ax2=ax.twinx(); ax2.semilogy(vscan,absdef,ls="--",alpha=.65,label=r"absolute $\delta_{\rm CP}$"); ax2.set_ylabel("absolute negative mass")
    lines=ax.get_lines()+ax2.get_lines(); ax.legend(lines,[l.get_label() for l in lines],frameon=False,fontsize=8,loc="upper right")

    ax=axes[1]
    Ns=np.unique(np.logspace(1,5,80).astype(int))
    def fe(z): return 1+z
    def ssp(z): return 1+z+z**2/2+z**3/6
    for name,fun in [("Forward Euler",fe),("SSPRK(3,3)",ssp)]:
        vals=[]
        for N in Ns:
            h=mp.mpf(1)/N; A=fun(-h)**N; C=fun(-h/2)**N
            lam=(1+A-mp.sqrt((1-A)**2+4*abs(C)**2))/4
            vals.append(float(-lam))
        ax.loglog(1/Ns,vals,"o-",ms=2,lw=1,label=name)
    ax.set(xlabel=r"step $x_N=1/N$",ylabel=r"negative Bell probability magnitude $-p_{-,N}^{\rm RK}$",title="(b) Outward methods converge from outside the channel set")
    ax.grid(which="both",alpha=.18); ax.legend(frameon=False,fontsize=8)
    fig.tight_layout(); save(fig,"Fig5_operational_refinement")
    with (DATA / "figure5_operational.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["varpi","ratio","absolute_defect"]); w.writerows(zip(vscan,ratios,absdef))


def fig6_frequency_bands():
    fig,ax=plt.subplots(figsize=(7.3,3.25))
    vmax=2.3; vlo=np.sqrt(2-np.sqrt(33)/4); vr=np.sqrt(3)/2; vhi=np.sqrt(2+np.sqrt(33)/4)
    rows=[("Dormand-Prince 5",[(0,vlo,"inward"),(vlo,vhi,"outward"),(vhi,vmax,"inward")]),("classical RK4",[(0,vr,"inward"),(vr,vmax,"outward")])]
    for y,(name,segs) in enumerate(rows):
        for lo,hi,state in segs:
            ax.barh(y, hi-lo, left=lo, height=.52, alpha=.38 if state=="inward" else .72,
                    hatch="//" if state=="outward" else "")
            ax.text((lo+hi)/2,y,state,ha="center",va="center",fontsize=8)
    ax.set_yticks([0,1]); ax.set_yticklabels([r[0] for r in rows])
    for val,label in [(vlo,"DP5−"),(vr,"RK4"),(vhi,"DP5+")]: ax.axvline(val,ls=":",lw=.8); ax.text(val+.012,1.55,label,rotation=90,va="top",fontsize=7)
    ax.set(xlim=(0,vmax),ylim=(-.6,1.65),xlabel=r"frequency ratio $\varpi=|\omega|/\Gamma$",title="First-defect orientation bands on the one-way boundary")
    ax.grid(axis="x",alpha=.16); fig.tight_layout(); save(fig,"Fig6_frequency_bands")


def graphical_abstract():
    fig=plt.figure(figsize=(12.5,5.0)); gs=fig.add_gridspec(1,3,width_ratios=[1.25,.10,1.0])
    ax=fig.add_subplot(gs[0,0]); v=np.linspace(0,4,1200); exact_phase_fill(ax,v); ax.set(xlim=(0,4),ylim=(0,2.9),xlabel=r"representation frequency $|\omega|/\Gamma$",ylabel=r"step $\Gamma h$",title="Direct RK4 has disconnected physical windows")
    mid=fig.add_subplot(gs[0,1]); mid.axis("off"); mid.text(.5,.5,r"$\Longrightarrow$",ha="center",va="center",fontsize=26)
    ax3=fig.add_subplot(gs[0,2]); ax3.axis("off")
    ax3.text(.02,.86,"Certify, transform, or switch",fontsize=15,weight="bold")
    ax3.text(.02,.69,"Exact membership: two stability-function values\nfor a known phase-covariant qubit generator.",fontsize=12.5)
    ax3.text(.02,.48,"A PI error controller can accept non-CPTP steps.\nA certified guard returns PASS / FAIL / UNCERTAIN.",fontsize=12.5)
    ax3.text(.02,.27,"A rotating frame removes the high-frequency pathology\nin this benchmark; a transverse field does not remove\nthe disconnected topology.",fontsize=12.2)
    ax3.text(.02,.08,"Every boundary and sign change is regenerated\nfrom exact or certified source data.",fontsize=11.8)
    fig.tight_layout(); fig.savefig(FIG/"Graphical_Abstract.pdf",metadata=PDF_META); fig.savefig(FIG/"Graphical_Abstract.png",dpi=500,metadata=PNG_META); fig.savefig(FIG/"Graphical_Abstract.eps"); plt.close(fig)


if __name__ == "__main__":
    fig1_topology_conditioning()
    fig2_robustness_noncommuting()
    fig3_certified_guard_conditioning()
    fig4_adaptive_rotating_frame()
    fig5_operational_refinement()
    fig6_frequency_bands()
    graphical_abstract()
    print("Generated six v2.7 JCP figures, graphical abstract, and source data.")
