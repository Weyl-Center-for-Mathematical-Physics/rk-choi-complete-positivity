#!/usr/bin/env python3
from __future__ import annotations

"""Independent v2.8 audit.

This script intentionally does not import ``rk_choi_margin``.  It reconstructs
all formulas, Choi matrices, reduction thresholds, and noncommuting samples in a
separate implementation.
"""

import json
import math
import random
from itertools import combinations, product
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "independent_validation" / "independent_audit_v28_report.json"


def r4(z):
    return 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24


def phase_choi(a, c, theta):
    A = 1 - theta * (1 - a)
    B = theta * (1 - a)
    C = (1 - theta) * (1 - a)
    D = theta + (1 - theta) * a
    return np.array([[A,0,0,c],[0,B,0,0],[0,0,C,0],[np.conjugate(c),0,0,D]],dtype=np.complex128)


def vecf(M):
    return np.asarray(M,dtype=np.complex128).reshape(-1,order="F")


def choi_from_super(S, d=2):
    J=np.zeros((d*d,d*d),dtype=np.complex128)
    for j,k in product(range(d),repeat=2):
        E=np.zeros((d,d),dtype=np.complex128); E[j,k]=1
        O=(S@vecf(E)).reshape((d,d),order="F")
        J[j*d:(j+1)*d,k*d:(k+1)*d]=O
    return (J+J.conj().T)/2


def dissipator_super(L):
    I=np.eye(2,dtype=np.complex128); Q=L.conj().T@L
    return np.kron(np.conjugate(L),L)-np.kron(I,Q)/2-np.kron(Q.T,I)/2


def gksl_super(theta=.001,kappa=.001,omega_z=2,omega_x=.1):
    I=np.eye(2,dtype=np.complex128); sx=np.array([[0,1],[1,0]],complex); sz=np.diag([1,-1]).astype(complex)
    sm=np.array([[0,1],[0,0]],complex); spm=sm.T
    H=-omega_z*sz/2+omega_x*sx/2
    L=-1j*(np.kron(I,H)-np.kron(H.T,I))
    L+=(1-theta)*dissipator_super(sm)+theta*dissipator_super(spm)+kappa/2*dissipator_super(sz)
    return L


def rk4_super(L,h):
    Z=h*L; I=np.eye(L.shape[0],dtype=complex)
    return I+Z+Z@Z/2+Z@Z@Z/6+Z@Z@Z@Z/24


def main():
    checks=[]
    def check(name, condition, detail=""):
        checks.append({"name":name,"pass":bool(condition),"detail":detail})
        if not condition: raise AssertionError(name+": "+detail)

    x,q,y=sp.symbols("x q y",real=True)
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
        disc=mp.sqrt((A-D)**2+4*abs(cc)**2)
        eig=[B,Cc,(A+D-disc)/2,(A+D+disc)/2]
        spectral=min(eig)>=0
        if scalar!=spectral: mismatches+=1
    check("15000 scalar-vs-Choi comparisons",mismatches==0,f"mismatches={mismatches}")

    # The four false accepts have exact negative margins.
    false_rows=[]
    for text in ('0.004','0.001','0.0001','0.00001'):
        xx=sp.Rational(text); qq=sp.Integer(17); yy=qq*xx
        exact=sp.factor(-xx**5*qq*(yy**3-16*yy**2-32*(qq-6)*yy+384*(qq-4))/147456)
        false_rows.append({"x":text,"exact_margin":str(exact),"negative":bool(exact<0)})
    check("small-step false-accept signs",all(r["negative"] for r in false_rows),str(false_rows))

    # Candidate-map calculus: full and two-half candidates are non-nested.
    def pc_status(total, substeps, vv):
        hh=mp.mpf(total)/substeps
        aa=r4(-hh)**substeps
        cc=r4(-hh*(mp.mpf(".5")-1j*mp.mpf(vv)))**substeps
        return aa <= 1 and aa >= 0 and aa-abs(cc)**2 >= 0
    check("H=1 full pass fine fail", pc_status(1,1,2) and not pc_status(1,2,2))
    check("H=2 full fail fine pass", not pc_status(2,1,2) and pc_status(2,2,2))

    # Exact Richardson factorization.
    ac=R(-x); cc=R(-x/2)
    af=R(-x/2)**2; cf=R(-x/4)**2
    ae=sp.factor((16*af-ac)/15); ce=sp.factor((16*cf-cc)/15)
    mext=sp.factor(ae-ce**2)
    p10=(x**10-64*x**9+2304*x**8-59392*x**7+1183744*x**6-19169280*x**5+258932736*x**4-2960916480*x**3+19888865280*x**2-97391738880*x+280850595840)
    target_ext=sp.factor(-x**6*p10/sp.Integer(1252412463513600))
    check("Richardson margin factorization", sp.factor(mext-target_ext)==0)
    check("Richardson leading coefficient", sp.expand(mext).coeff(x,6)==-sp.Rational(31,138240))
    alpha4=sp.N(sp.nroots(x**3-4*x**2+12*x-24)[0],50)
    p10_roots=sp.Poly(p10,x,domain=sp.QQ).intervals()
    check("Richardson P10 positive through alpha4", all(not (sp.Rational(lo)<alpha4 and sp.Rational(hi)>0) for (lo,hi),mult in p10_roots))

    # Executed binary input and displayed decimal lie on opposite sides.
    text="1.270334626973889"
    f=float(text); bf=sp.Rational(*f.as_integer_ratio()); df=sp.Rational(text)
    def exact_margin(xx):
        return sp.factor(R(-xx)-R(-xx*(sp.Rational(1,2)-2*sp.I))*R(-xx*(sp.Rational(1,2)+2*sp.I)))
    check("binary-vs-decimal certification reversal", exact_margin(bf)<0 and exact_margin(df)>0)

    # Noncommuting direct construction, independent of package implementation.
    L=gksl_super(); noncomm=[]
    for h,expected in [(0.1,1),(0.5,-1),(0.9,1),(1.4,-1)]:
        eig=np.linalg.eigvalsh(choi_from_super(rk4_super(L,h)))
        sign=1 if eig[0]>0 else -1 if eig[0]<0 else 0
        noncomm.append({"h":h,"lambda_min":float(eig[0]),"sign":sign})
        check(f"noncommuting sign h={h}",sign==expected,str(eig[0]))

    # Rotating-frame RK4: CP status is that of nonrotating dissipative RK4.
    # At h=.7, lab-frame varpi=2 is non-CP while de-rotated RK4 is CP.
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
        "false_accept_rows":false_rows,
        "noncommuting_samples":noncomm,
        "varpi_half":str(sp.N(varpi_half,30)),
    }
    OUT.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(f"Independent audit v2.8: {report['passed']}/{report['total']} checks passed; 15000 randomized mismatches={mismatches}.")
    print(f"Wrote {OUT}")

if __name__=='__main__': main()
