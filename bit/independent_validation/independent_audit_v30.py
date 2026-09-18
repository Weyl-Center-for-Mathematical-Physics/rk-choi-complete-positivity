#!/usr/bin/env python3
from __future__ import annotations

"""Standalone analytic and direct-Choi audit for v3.0.

This script deliberately does not import rk_choi_margin.  It reconstructs the
phase-covariant criterion, RK4 topology, Richardson candidate, noncommuting
Choi samples, and frame comparison independently.
"""

import json
from pathlib import Path
import random

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "independent_validation" / "independent_audit_v30_report.json"


def r4(z):
    return 1 + z + z**2/2 + z**3/6 + z**4/24


def vec(m):
    return np.asarray(m, dtype=complex).reshape(4, order="F")


def unvec(v):
    return np.asarray(v, dtype=complex).reshape((2, 2), order="F")


def gksl_super(gamma=1.0, theta=1e-3, kappa=1e-3, omega_z=2.0, omega_x=0.1):
    sm = np.array([[0,1],[0,0]],dtype=complex)
    spm = sm.conj().T
    sz = np.array([[1,0],[0,-1]],dtype=complex)
    sx = np.array([[0,1],[1,0]],dtype=complex)
    H = -omega_z*sz/2 + omega_x*sx/2
    gd = gamma*(1-theta); gu=gamma*theta; gp=gamma*kappa
    basis=[]
    for j in range(2):
        for i in range(2):
            E=np.zeros((2,2),complex); E[i,j]=1; basis.append(E)
    cols=[]
    for X in basis:
        Y=-1j*(H@X-X@H)
        for rate,L in [(gd,sm),(gu,spm)]:
            A=L.conj().T@L
            Y += rate*(L@X@L.conj().T-(A@X+X@A)/2)
        Y += gp/2*(sz@X@sz-X)
        cols.append(vec(Y))
    return np.column_stack(cols)


def choi_from_super(S):
    J=np.zeros((4,4),complex)
    for i in range(2):
        for j in range(2):
            E=np.zeros((2,2),complex); E[i,j]=1
            Phi=unvec(S@vec(E))
            J += np.kron(E,Phi)
    return (J+J.conj().T)/2


def rk4_super(L,h):
    Z=h*L; I=np.eye(L.shape[0],dtype=complex)
    return I+Z+Z@Z/2+Z@Z@Z/6+Z@Z@Z@Z/24


def main():
    checks=[]
    def check(name, condition, detail=""):
        checks.append({"name":name,"pass":bool(condition),"detail":detail})
        if not condition: raise AssertionError(name+": "+detail)

    x,q,y,s=sp.symbols("x q y s",real=True)
    v=sp.symbols("v",real=True)
    R=lambda z:1+z+z**2/sp.Integer(2)+z**3/sp.Integer(6)+z**4/sp.Integer(24)

    a=R(-x); c=R(-x*(sp.Rational(1,2)-sp.I*v))
    M=sp.factor(sp.expand(a-c*sp.conjugate(c)))
    C=y**3-16*y**2-32*(q-6)*y+384*(q-4)
    target=sp.factor(-x**5*(1+4*v**2)*C.subs({q:1+4*v**2,y:(1+4*v**2)*x})/147456)
    check("exact RK4 rotating factorization",sp.simplify(M-target)==0)
    disc=sp.factor(sp.discriminant(C,y))
    check("cubic discriminant",sp.simplify(disc-16384*(q-4)*(8*q**2-123*q+348))==0,str(disc))
    resultant=sp.factor(sp.resultant(C,C.subs(y,2*y),y))
    check("halving resultant",sp.simplify(resultant-6291456*(q-4)*(72*q**3-2071*q**2+12552*q-21744))==0,str(resultant))

    poly=sp.Poly(72*q**3-2071*q**2+12552*q-21744,q,domain=sp.QQ)
    qplus=(sp.Integer(123)+11*sp.sqrt(33))/16
    roots=[r for r in sp.nroots(poly.as_expr(),n=70,maxsteps=500) if abs(sp.im(r))<sp.Float('1e-60') and sp.re(r)>qplus]
    check("unique detached halving threshold root",len(roots)==1,str(roots))
    varpi_half=sp.sqrt((sp.re(roots[0])-1)/4)
    check("halving threshold value",abs(float(varpi_half)-2.2482546045149805)<1e-14,str(varpi_half))

    # Random direct-Choi check.
    rng=random.Random(20260810); mp.mp.dps=100; mismatches=0
    for _ in range(15000):
        xx=mp.mpf(rng.randint(1,300000))/100000
        th=mp.mpf(rng.randint(0,1000))/1000
        kap=mp.mpf(rng.randint(0,200))/1000
        vv=mp.mpf(rng.randint(0,5000))/1000
        aa=r4(-xx); cc=r4(-xx*(mp.mpf('.5')+kap-1j*vv))
        A=1-th*(1-aa); B=th*(1-aa); Cc=(1-th)*(1-aa); D=th+(1-th)*aa
        F=A*D-abs(cc)**2
        scalar=(1-aa>=0 and A>=0 and D>=0 and F>=0)
        d=mp.sqrt((A-D)**2+4*abs(cc)**2)
        eig=[B,Cc,(A+D-d)/2,(A+D+d)/2]
        spectral=min(eig)>=0
        if scalar!=spectral: mismatches+=1
    check("15000 scalar-vs-Choi comparisons",mismatches==0,f"mismatches={mismatches}")

    # Candidate-map non-nesting.
    def pc_status(total, substeps, vv):
        hh=mp.mpf(total)/substeps
        aa=r4(-hh)**substeps
        cc=r4(-hh*(mp.mpf('.5')-1j*mp.mpf(vv)))**substeps
        return aa <= 1 and aa >= 0 and aa-abs(cc)**2 >= 0
    check("H=1 full pass fine fail", pc_status(1,1,2) and not pc_status(1,2,2))
    check("H=2 full fail fine pass", not pc_status(2,1,2) and pc_status(2,2,2))

    # Richardson stability function and exact certificate.
    Rext=sp.expand((16*R(s/2)**2-R(s))/15)
    expected=(1+s+s**2/2+s**3/6+s**4/24+s**5/120+s**6/864+s**7/8640+s**8/138240)
    check("Richardson stability function",sp.expand(Rext-expected)==0)
    eta6=sp.expand(Rext-sp.exp(s)).series(s,0,7).removeO().coeff(s,6)
    check("Richardson first defect",eta6==-sp.Rational(1,4320),str(eta6))
    G6=q*(q**2-18*q+48)/32
    check("Richardson nonrotating boundary coefficient",sp.factor(eta6*G6.subs(q,1))==-sp.Rational(31,138240))

    ac=R(-x); cc=R(-x/2)
    af=R(-x/2)**2; cf=R(-x/4)**2
    ae=sp.factor((16*af-ac)/15); ce=sp.factor((16*cf-cc)/15)
    mext=sp.factor(ae-ce**2)
    p10=(x**10-64*x**9+2304*x**8-59392*x**7+1183744*x**6-19169280*x**5+258932736*x**4-2960916480*x**3+19888865280*x**2-97391738880*x+280850595840)
    target_ext=sp.factor(-x**6*p10/sp.Integer(1252412463513600))
    check("Richardson margin factorization", sp.factor(mext-target_ext)==0)
    check("Richardson leading coefficient", sp.expand(mext).coeff(x,6)==-sp.Rational(31,138240))
    p10_positive_roots=sorted(float(sp.re(r)) for r in sp.nroots(p10,n=50,maxsteps=500) if abs(float(sp.im(r)))<1e-20 and float(sp.re(r))>0)
    check("Richardson first positive margin root",abs(p10_positive_roots[0]-6.034523958649)<1e-10,str(p10_positive_roots[:2]))

    # Binary and decimal boundary inputs differ.
    text="1.270334626973889"
    f=float(text); bf=sp.Rational(*f.as_integer_ratio()); df=sp.Rational(text)
    def exact_margin(xx):
        return sp.factor(R(-xx)-R(-xx*(sp.Rational(1,2)-2*sp.I))*R(-xx*(sp.Rational(1,2)+2*sp.I)))
    check("binary-vs-decimal certification reversal", exact_margin(bf)<0 and exact_margin(df)>0)

    # Noncommuting map-level Choi signs, independently reconstructed.
    L=gksl_super(); noncomm=[]
    for h,expected_sign in [(0.1,1),(0.5,-1),(0.9,1),(1.4,-1)]:
        eig=np.linalg.eigvalsh(choi_from_super(rk4_super(L,h)))
        sign=1 if eig[0]>0 else -1 if eig[0]<0 else 0
        noncomm.append({"h":h,"lambda_min":float(eig[0]),"sign":sign})
        check(f"noncommuting sign h={h}",sign==expected_sign,str(eig[0]))

    # Frame dependence.
    h=mp.mpf('.7'); aa=r4(-h); clab=r4(-h*(mp.mpf('.5')-2j)); crot=r4(-h/2)
    lab=aa-abs(clab)**2; rot=aa-abs(crot)**2
    check("frame dependence at h=.7",lab<0 and rot>0,f"lab={lab}, rotating={rot}")

    report={
        "implementation":"standalone: no rk_choi_margin imports",
        "checks":checks,
        "passed":sum(c["pass"] for c in checks),
        "total":len(checks),
        "randomized_samples":15000,
        "randomized_mismatches":mismatches,
        "noncommuting_samples":noncomm,
        "varpi_half":str(sp.N(varpi_half,30)),
    }
    OUT.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(f"Independent audit v3.0: {report['passed']}/{report['total']} checks passed; 15000 randomized mismatches={mismatches}.")
    print(f"Wrote {OUT}")


if __name__=='__main__':
    main()
