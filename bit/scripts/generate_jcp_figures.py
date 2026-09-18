#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import sympy as sp

from rk_choi_margin.topology import normalized_choi_defect_and_error
from rk_choi_margin.v26 import halving_sequence, rk4_cptp_guard_step, rk4_rotating_positive_roots

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures_jcp"
DATA = ROOT / "results" / "jcp_figures"
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

PDF_META = {
    "Title": "Runge-Kutta complete-positivity figures",
    "Author": "G. Blake Pierpoint, Olivier Bernard, Yichen Liu",
    "Creator": "generate_jcp_figures.py",
}


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight", metadata=PDF_META)
    fig.savefig(FIG / f"{stem}.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.eps", bbox_inches="tight")
    plt.close(fig)


def r4(z):
    return 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24


def phase_margin(theta, kappa, varpi, xs):
    a = r4(-xs)
    c = r4(-xs * (0.5 + kappa - 1j * varpi))
    return (a + theta * (1 - theta) * (1 - a) ** 2 - np.abs(c) ** 2).real


def rk4_intervals(varpi: float, alpha4: float = 2.785293563405282):
    """Plotting intervals from exact cubic roots; isolated points returned separately."""
    qminus = (123 - 11 * np.sqrt(33)) / 16
    qplus = (123 + 11 * np.sqrt(33)) / 16
    q = 1 + 4 * varpi**2
    eps = 1e-8
    roots = rk4_rotating_positive_roots(varpi)
    if q < qminus - eps:
        if not roots:
            return [], []
        return [(0.0, min(alpha4, roots[-1]))], []
    if abs(q - qminus) <= eps:
        # Connected interval; the extra double root is an internal tangency.
        upper = min(alpha4, (2 * np.sqrt(33) - 2) / qminus)
        return [(0.0, upper)], [(9 - np.sqrt(33)) / qminus]
    if q < 4 - eps:
        if len(roots) != 3:
            return [], []
        return [(0.0, roots[0]), (roots[1], min(alpha4, roots[2]))], []
    if abs(q - 4) <= eps:
        return [], [0.0, 2.0]
    if q < qplus - eps:
        return [], [0.0]
    if abs(q - qplus) <= eps:
        return [], [0.0, (9 + np.sqrt(33)) / qplus]
    if len(roots) != 2:
        return [], [0.0]
    return [(roots[0], min(alpha4, roots[1]))], [0.0]


def exact_phase_fill(ax, vgrid):
    attached_upper = np.full_like(vgrid, np.nan, dtype=float)
    detached_lower = np.full_like(vgrid, np.nan, dtype=float)
    detached_upper = np.full_like(vgrid, np.nan, dtype=float)
    rows = []
    for i, vv in enumerate(vgrid):
        intervals, points = rk4_intervals(float(vv))
        rows.append((float(vv), intervals, points))
        if intervals:
            if intervals[0][0] == 0:
                attached_upper[i] = intervals[0][1]
                if len(intervals) > 1:
                    detached_lower[i], detached_upper[i] = intervals[1]
            else:
                detached_lower[i], detached_upper[i] = intervals[0]
    ax.fill_between(vgrid, 0, attached_upper, where=np.isfinite(attached_upper), alpha=0.30, label="attached CPTP component")
    ax.fill_between(vgrid, detached_lower, detached_upper, where=np.isfinite(detached_lower), alpha=0.42, label="detached CPTP component")
    ax.plot(vgrid, attached_upper, linewidth=1.25)
    ax.plot(vgrid, detached_lower, linewidth=1.25)
    ax.plot(vgrid, detached_upper, linewidth=1.25)
    return rows, attached_upper, detached_lower, detached_upper


def figure1_phase_diagram():
    vm = float(sp.N(sp.sqrt(107 - 11 * sp.sqrt(33)) / 8, 16))
    vu = np.sqrt(3) / 2
    vp = float(sp.N(sp.sqrt(107 + 11 * sp.sqrt(33)) / 8, 16))
    vgrid = np.linspace(0, 4.0, 1601)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.35), gridspec_kw={"width_ratios": [1.65, 1.0]})
    ax = axes[0]
    rows, au, dl, du = exact_phase_fill(ax, vgrid)
    alpha4 = 2.785293563405282
    ax.axhline(alpha4, linestyle="--", linewidth=0.9, label=r"population ceiling $x=\alpha_4$")
    for value in (vm, vu, vp):
        ax.axvline(value, linestyle=":", linewidth=0.9)
    ax.scatter([vu, vp], [2.0, (9 + np.sqrt(33)) / ((123 + 11*np.sqrt(33))/16)], s=18, facecolors="white", edgecolors="black", zorder=5)
    ax.text(vm - 0.025, 2.92, r"$\varpi_-$", rotation=90, va="top", ha="right", fontsize=8)
    ax.text(vu + 0.020, 2.92, r"$\sqrt{3}/2$", rotation=90, va="top", ha="left", fontsize=8)
    ax.text(vp + 0.020, 2.92, r"$\varpi_+$", rotation=90, va="top", ha="left", fontsize=8)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 3)
    ax.set_xlabel(r"frequency ratio $\varpi=|\omega|/\Gamma$")
    ax.set_ylabel(r"dimensionless step $x=\Gamma h$")
    ax.set_title("(a) Exact RK4 CPTP topology")
    ax.legend(loc="upper right", frameon=False, fontsize=8)

    inset = ax.inset_axes([0.06, 0.10, 0.47, 0.46])
    vz = np.linspace(0.80, 0.88, 901)
    exact_phase_fill(inset, vz)
    inset.axvline(vm, linestyle=":", linewidth=0.7)
    inset.axvline(vu, linestyle=":", linewidth=0.7)
    inset.set_xlim(0.80, 0.88)
    inset.set_ylim(0, 2.9)
    inset.set_xticks([0.80, 0.83, 0.86, 0.88])
    inset.set_yticks([0, 1, 2])
    inset.tick_params(labelsize=7)
    inset.set_title("reentrant band", fontsize=8)

    ax = axes[1]
    vals = np.geomspace(vp + 1e-3, 100, 320)
    lower, upper = [], []
    for vv in vals:
        roots = rk4_rotating_positive_roots(float(vv))
        lower.append(roots[0] if len(roots) == 2 else np.nan)
        upper.append(roots[1] if len(roots) == 2 else np.nan)
    lower = np.asarray(lower); upper = np.asarray(upper)
    ax.semilogx(vals, lower * vals**2, label=r"$x_{-}\varpi^2$")
    ax.semilogx(vals, upper * vals, label=r"$x_{+}\varpi$")
    ax.axhline(3.0, linestyle="--", linewidth=0.9, label=r"limit $3$")
    ax.axhline(2 * np.sqrt(2), linestyle=":", linewidth=0.9, label=r"limit $2\sqrt{2}$")
    ax.set_xlabel(r"$\varpi$")
    ax.set_ylabel("scaled endpoint")
    ax.set_title("(b) Rotation-dominated limits")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.20)
    fig.tight_layout()
    save(fig, "Fig1_RK4_phase_diagram")

    with (DATA / "figure1_exact_intervals.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["varpi", "attached_upper", "detached_lower", "detached_upper"])
        for vv, aui, dli, dui in zip(vgrid, au, dl, du):
            writer.writerow([vv, aui, dli, dui])


def saddle_node_curve_varpi2():
    # Continue F=0 and dF/dx=0 at varpi=2.  This is the curve on which the
    # intermediate non-CPTP gap closes.  Exact rational examples are certified
    # separately; this curve is a high-precision visualization of the topology.
    mp.mp.dps = 60
    xs, ts, ks = [], [], []

    def rr(z):
        return 1 + z + z**2/2 + z**3/6 + z**4/24

    def f(xx, tt, kk):
        a = rr(-xx)
        c = rr(-xx*(mp.mpf('0.5') + kk - 2j))
        return a + tt*(1-tt)*(1-a)**2 - abs(c)**2

    def dfdx(xx, tt, kk):
        return mp.diff(lambda zz: f(zz, tt, kk), xx)

    xg = mp.mpf('0.561716884769793')
    tg = mp.mpf('0.0378962854729975')
    kcrit = mp.mpf('0.0124665355913517154')
    for kk in [kcrit * i / 160 for i in range(161)]:
        try:
            root = mp.findroot(lambda xx, tt: (f(xx, tt, kk), dfdx(xx, tt, kk)), (xg, tg), tol=mp.mpf('1e-38'), maxsteps=100)
            xg, tg = root[0], root[1]
            if tg < -mp.mpf('1e-10'):
                break
            xs.append(float(xg)); ts.append(max(0.0, float(tg))); ks.append(float(kk))
        except Exception:
            break
    return np.asarray(ks), np.asarray(ts), np.asarray(xs)


def figure2_robustness():
    cases = [
        (r"boundary $\theta=\kappa=0$", [(0.744782427896, 1.270334626974)]),
        (r"$\kappa=10^{-3}$", [(0.0, 0.253348869116), (0.740695999652, 1.269542711140)]),
        (r"$\theta=10^{-3}$", [(0.0, 0.122268903470), (0.742558312554, 1.270544370136)]),
        (r"$\theta=\kappa=10^{-3}$", [(0.0, 0.262033333161), (0.738364627214, 1.269753570827)]),
        (r"$\kappa=10^{-2}$", [(0.0, 0.519308218820), (0.681621610068, 1.262280154980)]),
        (r"$\theta=10^{-2}$", [(0.0, 0.276930275811), (0.720896805059, 1.272391394935)]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3), gridspec_kw={"width_ratios": [1.1, 1.0]})
    ax = axes[0]
    ypos = np.arange(len(cases))[::-1]
    for y, (label, intervals) in zip(ypos, cases):
        for lo, hi in intervals:
            ax.hlines(y, lo, hi, linewidth=5)
            ax.plot([lo, hi], [y, y], "o", markersize=3.5)
    ax.set_yticks(ypos)
    ax.set_yticklabels([c[0] for c in cases], fontsize=8)
    ax.set_xlim(-0.02, 1.34)
    ax.set_xlabel(r"$x=\Gamma h$ at $\varpi=2$")
    ax.set_title("(a) Certified components persist off the boundary")
    ax.grid(axis="x", alpha=0.18)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    ks, ts, saddle_x = saddle_node_curve_varpi2()
    ax = axes[1]
    ax.fill_between(ks, 0, ts, alpha=0.28, label="two positive-step components")
    ax.plot(ks, ts, linewidth=1.4, label=r"gap-closing curve $F=\partial_xF=0$")
    cert_points = [(0.001,0.0),(0.0,0.001),(0.001,0.001),(0.01,0.0),(0.0,0.01)]
    # plot coordinates as (kappa, theta)
    for theta, kappa in cert_points:
        ax.plot(kappa, theta, "o", markersize=4)
    ax.set_xlim(0, max(ks)*1.04)
    ax.set_ylim(0, max(ts)*1.05)
    ax.set_xlabel(r"dephasing ratio $\kappa=\gamma_\phi/\Gamma$")
    ax.set_ylabel(r"thermal fraction $\theta$")
    ax.set_title("(b) Numerically continued saddle-node boundary")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "Fig2_robustness")

    with (DATA / "figure2_saddle_node_curve.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle); writer.writerow(["kappa", "theta_critical", "double_root_x"])
        writer.writerows(zip(ks, ts, saddle_x))


def figure3_operational_error():
    mp.mp.dps = 80
    steps = np.geomspace(1e-4, 2e-1, 150)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax = axes[0]
    for vv in (1.0, 2.0, 5.0):
        defects = []
        for xx in steps:
            d, _ = normalized_choi_defect_and_error(float(xx), vv)
            defects.append(float(d) if d > 0 else np.nan)
        ax.loglog(steps, defects, label=fr"$\varpi={vv:g}$")
    ax.set_xlabel(r"$x=\Gamma h$")
    ax.set_ylabel(r"negative Born probability magnitude $\delta_{\rm CP}$")
    ax.set_title("(a) Ancilla-assisted positivity defect")
    ax.legend(frameon=False)
    ax.grid(which="both", alpha=0.20)

    ax = axes[1]
    vscan = np.linspace(0.9, 10.0, 460)
    ratios=[]
    for vv in vscan:
        d,e = normalized_choi_defect_and_error(1e-2, float(vv))
        ratios.append(float(d/e) if (e and d>0) else np.nan)
    ratios=np.asarray(ratios)
    ax.plot(vscan, ratios)
    d2,e2=normalized_choi_defect_and_error(1e-2,2.0)
    r2=float(d2/e2)
    ax.plot(2.0,r2,"o")
    ax.annotate(fr"$\varpi=2$: {r2:.3f}", xy=(2.0,r2), xytext=(2.8,0.80), arrowprops=dict(arrowstyle="->", linewidth=0.8), fontsize=8)
    ax.set_ylim(0,1.0)
    ax.set_xlabel(r"frequency ratio $\varpi$ at $x=10^{-2}$")
    ax.set_ylabel(r"$\delta_{\rm CP}/E_J$")
    ax.set_title("(b) Physicality defect as a fraction of total channel error")
    ax.grid(alpha=0.20)
    fig.tight_layout()
    save(fig,"Fig3_operational_error")

    with (DATA / "figure3_operational_error.csv").open("w", newline="", encoding="utf-8") as handle:
        writer=csv.writer(handle); writer.writerow(["x","varpi","defect","choi_trace_distance","ratio"])
        for vv in (1.0,2.0,5.0):
            for xx in steps:
                d,e=normalized_choi_defect_and_error(float(xx),vv)
                writer.writerow([xx,vv,mp.nstr(d,22),mp.nstr(e,22),mp.nstr(d/e,18) if e else "nan"])
    with (DATA / "figure3_ratio_vs_frequency.csv").open("w", newline="", encoding="utf-8") as handle:
        writer=csv.writer(handle); writer.writerow(["varpi","ratio_at_x_1e-2"])
        writer.writerows(zip(vscan,ratios))


def lower_eig(A,C):
    return (1 + A - mp.sqrt((1-A)**2 + 4*abs(C)**2))/4


def figure4_frequency_and_refinement():
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.1))
    ax=axes[0]
    vmax=2.3
    v_dp_lo=np.sqrt(2-np.sqrt(33)/4)
    v_rk4=np.sqrt(3)/2
    v_dp_hi=np.sqrt(2+np.sqrt(33)/4)
    rows=[("classical RK4",[(0,v_rk4,"inward"),(v_rk4,vmax,"outward")]),
          ("Dormand-Prince 5",[(0,v_dp_lo,"inward"),(v_dp_lo,v_dp_hi,"outward"),(v_dp_hi,vmax,"inward")])]
    for y,(name,segs) in enumerate(rows[::-1]):
        for lo,hi,state in segs:
            ax.barh(y,hi-lo,left=lo,height=.48,alpha=.38 if state=="inward" else .70)
            ax.text((lo+hi)/2,y,state,ha="center",va="center",fontsize=8)
    ax.set_yticks([0,1]); ax.set_yticklabels([rows[1][0],rows[0][0]])
    for val,label in [(v_dp_lo,r"DP5 lower"),(v_rk4,r"RK4"),(v_dp_hi,r"DP5 upper")]:
        ax.axvline(val,linestyle=":",linewidth=.8)
        ax.text(val+.015,1.42,label,rotation=90,va="top",fontsize=7)
    ax.set_xlim(0,vmax); ax.set_ylim(-.6,1.55)
    ax.set_xlabel(r"frequency ratio $\varpi=|\omega|/\Gamma$")
    ax.set_title("(a) Small-step Choi-boundary orientation")
    ax.grid(axis="x",alpha=.15)

    mp.mp.dps=80
    Ns=np.unique(np.logspace(1,5,80).astype(int))
    ax=axes[1]
    methods=[("Forward Euler",lambda z:1+z),("SSPRK(3,3)",lambda z:1+z+z**2/2+z**3/6)]
    for name,fun in methods:
        vals=[]
        for N in Ns:
            h=mp.mpf(1)/N; A=fun(-h)**N; C=fun(-h/2)**N
            vals.append(float(-lower_eig(A,C)))
        ax.loglog(1/Ns,vals,marker="o",markersize=2,linewidth=1.0,label=name)
    ax.set_xlabel(r"step $x_N=1/N$")
    ax.set_ylabel(r"negative Born probability magnitude $-p_{-,N}^{\rm RK}$")
    ax.set_title("(b) Convergence from outside the channel set")
    ax.legend(frameon=False); ax.grid(which="both",alpha=.20)
    fig.tight_layout(); save(fig,"Fig4_frequency_refinement")


def figure5_controller():
    varpi=2.0
    lo,hi=rk4_rotating_positive_roots(varpi)
    x0=1.4
    seq=halving_sequence(x0,6)
    decision=rk4_cptp_guard_step(x0,varpi)
    fig,ax=plt.subplots(figsize=(7.2,3.35))
    ax.axvspan(lo,hi,alpha=.30,label="CPTP interval")
    y=0.55
    for k,xx in enumerate(seq):
        ax.plot(xx,y,"o",markersize=5)
        ax.text(xx,y+.06,str(k),ha="center",fontsize=7)
        if k<len(seq)-1:
            ax.annotate("",xy=(seq[k+1],y),xytext=(xx,y),arrowprops=dict(arrowstyle="->",linewidth=.8))
    if decision.accepted_x is not None:
        ax.plot(decision.accepted_x,0.20,"s",markersize=6,label="CPTP-aware choice")
        ax.annotate("guarded choice",xy=(decision.accepted_x,.20),xytext=(1.08,.08),arrowprops=dict(arrowstyle="->",linewidth=.8),fontsize=8)
    ax.axvline(lo,linestyle=":",linewidth=.8); ax.axvline(hi,linestyle=":",linewidth=.8)
    ax.text(lo,.86,r"$x_{-}$",ha="center",fontsize=8); ax.text(hi,.86,r"$x_{+}$",ha="center",fontsize=8)
    ax.text(.15,.64,"naive halving",fontsize=8)
    ax.set_xlim(0,1.48); ax.set_ylim(0,1.0); ax.set_yticks([])
    ax.set_xlabel(r"dimensionless step $x=\Gamma h$ at $\varpi=2$")
    ax.set_title("A monotone step-reduction rule can skip the physical window")
    ax.legend(frameon=False,loc="upper left",fontsize=8)
    ax.grid(axis="x",alpha=.15)
    fig.tight_layout(); save(fig,"Fig5_controller")
    with (DATA / "figure5_controller.csv").open("w",newline="",encoding="utf-8") as handle:
        w=csv.writer(handle); w.writerow(["quantity","value"]); w.writerow(["x_minus",lo]); w.writerow(["x_plus",hi]);
        for k,xx in enumerate(seq): w.writerow([f"halving_{k}",xx])
        w.writerow(["guarded_choice",decision.accepted_x])


def graphical_abstract():
    fig=plt.figure(figsize=(12.5,5.0))
    gs=fig.add_gridspec(1,3,width_ratios=[1.25,.10,1.0])
    ax=fig.add_subplot(gs[0,0])
    vgrid=np.linspace(0,4,1200)
    exact_phase_fill(ax,vgrid)
    ax.set_xlim(0,4); ax.set_ylim(0,2.9)
    ax.set_xlabel(r"rotation / relaxation $|\omega|/\Gamma$")
    ax.set_ylabel(r"step $\Gamma h$")
    ax.set_title("RK4 physical-channel windows can be disconnected")
    ax2=fig.add_subplot(gs[0,1]); ax2.axis("off"); ax2.text(.5,.5,r"$\Longrightarrow$",ha="center",va="center",fontsize=26)
    ax3=fig.add_subplot(gs[0,2]); ax3.axis("off")
    ax3.text(.02,.84,"Exact CPTP guard",fontsize=15,weight="bold")
    ax3.text(.02,.68,r"test $a\leq1$, $A_\theta\geq0$, $D_\theta\geq0$, $F_{R,\theta}\geq0$",fontsize=13)
    ax3.text(.02,.48,"A smaller step can leave the physical set;\nmonotone error control is not a physicality guarantee.",fontsize=13)
    ax3.text(.02,.25,r"high rotation: $3\Gamma/\omega^2\lesssim h\lesssim2\sqrt{2}/|\omega|$",fontsize=13.5,weight="bold")
    ax3.text(.02,.08,"Outside the window, a Bell-pair witness predicts\na negative Born probability.",fontsize=12.5)
    fig.tight_layout()
    fig.savefig(FIG / "Graphical_Abstract.pdf", metadata=PDF_META)
    fig.savefig(FIG / "Graphical_Abstract.png", dpi=450)
    fig.savefig(FIG / "Graphical_Abstract.eps")
    plt.close(fig)


if __name__ == "__main__":
    figure1_phase_diagram()
    figure2_robustness()
    figure3_operational_error()
    figure4_frequency_and_refinement()
    figure5_controller()
    graphical_abstract()
    print("Generated five JCP figures, graphical abstract, and source-data tables.")
