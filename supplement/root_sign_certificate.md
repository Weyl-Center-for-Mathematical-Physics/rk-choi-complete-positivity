# Supplement S1: Exact boundary polynomials, root isolation, and sign charts

This supplement records the exact root-isolation and sign-chart calculations used to determine the pure-amplitude-damping CPTP-admissible set for each of the eight Runge-Kutta formulas examined.

For a stability function \(R\), define

\[
a(x)=R(-x),
\qquad
M_R(x)=R(-x)-R(-x/2)^2,
\qquad x\ge 0.
\]

The direct one-step map is CPTP exactly when the Runge-Kutta stage equations are defined and

\[
a(x)\le 1,
\qquad
M_R(x)\ge 0.
\]

The second inequality implies \(a(x)\ge 0\), because \(M_R(x)=a(x)-R(-x/2)^2\).

For each method, the reproduction script derives \(R\) from exact rational Butcher coefficients, forms primitive integer numerator polynomials for \(a=0\), \(a=1\), and \(M_R=0\), verifies every listed positive-root bracket using a Sturm count, and evaluates all interval signs at exact rational sample points. Decimal bracket endpoints are interpreted as exact rational numbers. The complete scalar factors, denominators, multiplicities, samples, and truth values are stored in `data/certified_boundaries.json`.

## S1.1 Summary

| Method | Stage domain for \(x\ge 0\) | Certified CPTP-admissible set |
|---|---|---|
| Forward Euler | all \(x\ge 0\) | \(\{0\}\) |
| Heun RK2 | all \(x\ge 0\) | \([0,2]\) |
| SSPRK(3,3) | all \(x\ge 0\) | \(\{0\}\) |
| Classical RK4 | all \(x\ge 0\) | \([0,\alpha_4]\) |
| Dormand-Prince 5 principal formula | all \(x\ge 0\) | \([0,\alpha_5]\) |
| Dormand-Prince embedded 4 formula | all \(x\ge 0\) | \(\{0\}\cup[\beta_-,\beta_+]\) |
| Backward Euler | all \(x\ge 0\) | \([0,\infty)\) |
| Implicit midpoint | all \(x\ge 0\) | \(\{0\}\) |

## S1.2 Forward Euler

\[
R(s)=1+s.
\]

The full-step and half-step stage-resolvent determinants are both \(1\).

Primitive boundary numerators:

\[
P_0(x)=x-1,
\qquad
P_1(x)=x,
\qquad
P_M(x)=x^2.
\]

| Boundary | Certified root | Approximation | Mechanism |
|---|---:|---:|---|
| \(a=0\) | \(1\) | 1 | population floor |

| Open interval | \(\operatorname{sgn}a\) | \(\operatorname{sgn}(a-1)\) | \(\operatorname{sgn}M_R\) | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| \((0,1)\) | + | - | - | yes | no |
| \((1,\infty)\) | - | - | - | yes | no |

Therefore,

\[
\mathcal C_{\mathrm{FE}}=\{0\}.
\]

The first exponential defect is

\[
(m,\eta_m)=\left(2,-\frac12\right),
\]

and the leading margin coefficient is \(-1/4\).

## S1.3 Heun RK2

\[
R(s)=1+s+\frac{s^2}{2}.
\]

The full-step and half-step stage-resolvent determinants are both \(1\).

Primitive boundary numerators:

\[
P_0(x)=x^2-2x+2,
\qquad
P_1(x)=x(x-2),
\qquad
P_M(x)=x^3(x-8).
\]

| Boundary | Certified root | Approximation | Mechanism |
|---|---:|---:|---|
| \(a=1\) | \(2\) | 2 | population ceiling |
| \(M_R=0\) | \(8\) | 8 | coherence-population block |

| Open interval | \(\operatorname{sgn}a\) | \(\operatorname{sgn}(a-1)\) | \(\operatorname{sgn}M_R\) | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| \((0,2)\) | + | - | + | yes | yes |
| \((2,8)\) | + | + | + | yes | no |
| \((8,\infty)\) | + | + | - | yes | no |

Therefore,

\[
\mathcal C_{\mathrm{RK2}}=[0,2].
\]

The first exponential defect is

\[
(m,\eta_m)=\left(3,-\frac16\right),
\]

and the leading margin coefficient is \(1/8\).

## S1.4 SSPRK(3,3)

\[
R(s)=1+s+\frac{s^2}{2}+\frac{s^3}{6}.
\]

The full-step and half-step stage-resolvent determinants are both \(1\).

Primitive boundary numerators:

\[
P_0(x)=x^3-3x^2+6x-6,
\qquad
P_1(x)=x(x^2-3x+6),
\qquad
P_M(x)=x^4(x^2-12x+84).
\]

| Boundary | Certified rational bracket | Approximation | Mechanism |
|---|---|---:|---|
| \(a=0\) | \((1.596071637983,1.596071637984)\) | 1.596071637983322 | population floor |

| Open interval | \(\operatorname{sgn}a\) | \(\operatorname{sgn}(a-1)\) | \(\operatorname{sgn}M_R\) | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| \((0,1.596071637983322)\) | + | - | - | yes | no |
| \((1.596071637983322,\infty)\) | - | - | - | yes | no |

Therefore,

\[
\mathcal C_{\mathrm{SSP3}}=\{0\}.
\]

The first exponential defect is

\[
(m,\eta_m)=\left(4,-\frac1{24}\right),
\]

and the leading margin coefficient is \(-7/192\).

## S1.5 Classical RK4

\[
R(s)=1+s+\frac{s^2}{2}+\frac{s^3}{6}+\frac{s^4}{24}.
\]

The full-step and half-step stage-resolvent determinants are both \(1\).

Primitive boundary numerators:

\[
P_0(x)=x^4-4x^3+12x^2-24x+24,
\]

\[
P_1(x)=x(x^3-4x^2+12x-24),
\]

\[
P_M(x)=x^5(x^3-16x^2+160x-1152).
\]

| Label | Boundary | Certified rational bracket | Approximation | Mechanism |
|---|---|---|---:|---|
| \(\alpha_4\) | \(a=1\) | \((2.785293563405,2.785293563406)\) | 2.785293563405282 | population ceiling |
| \(\mu_4\) | \(M_R=0\) | \((10.982425466293,10.982425466294)\) | 10.982425466293273 | coherence-population block |

| Open interval | \(\operatorname{sgn}a\) | \(\operatorname{sgn}(a-1)\) | \(\operatorname{sgn}M_R\) | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| \((0,\alpha_4)\) | + | - | + | yes | yes |
| \((\alpha_4,\mu_4)\) | + | + | + | yes | no |
| \((\mu_4,\infty)\) | + | + | - | yes | no |

Therefore,

\[
\mathcal C_{\mathrm{RK4}}=[0,\alpha_4].
\]

The first exponential defect is

\[
(m,\eta_m)=\left(5,-\frac1{120}\right),
\]

and the leading margin coefficient is \(1/128\).

## S1.6 Dormand-Prince 5 principal formula

\[
R(s)=1+s+\frac{s^2}{2}+\frac{s^3}{6}+\frac{s^4}{24}+\frac{s^5}{120}+\frac{s^6}{600}.
\]

The full-step and half-step stage-resolvent determinants are both \(1\).

Primitive boundary numerators:

\[
P_0(x)=x^6-5x^5+25x^4-100x^3+300x^2-600x+600,
\]

\[
P_1(x)=x(x^5-5x^4+25x^3-100x^2+300x-600),
\]

\[
P_M(x)=x^6\left(x^6-20x^5+300x^4-3600x^3+35600x^2-294400x-396800\right).
\]

| Label | Boundary | Certified rational bracket | Approximation | Mechanism |
|---|---|---|---:|---|
| \(\alpha_5\) | \(a=1\) | \((3.306567892634,3.306567892635)\) | 3.306567892634947 | population ceiling |
| \(\mu_5\) | \(M_R=0\) | \((13.174093462800,13.174093462801)\) | 13.174093462800933 | coherence-population block |

| Open interval | \(\operatorname{sgn}a\) | \(\operatorname{sgn}(a-1)\) | \(\operatorname{sgn}M_R\) | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| \((0,\alpha_5)\) | + | - | + | yes | yes |
| \((\alpha_5,\mu_5)\) | + | + | + | yes | no |
| \((\mu_5,\infty)\) | + | + | - | yes | no |

Therefore,

\[
\mathcal C_{\mathrm{DP5}}=[0,\alpha_5].
\]

The first exponential defect is

\[
(m,\eta_m)=\left(6,\frac1{3600}\right),
\]

and the leading margin coefficient is \(31/115200\).

## S1.7 Dormand-Prince embedded fourth-order formula

\[
R(s)=1+s+\frac{s^2}{2}+\frac{s^3}{6}+\frac{s^4}{24}
+\frac{1097}{120000}s^5
+\frac{161}{120000}s^6
+\frac1{24000}s^7.
\]

The full-step and half-step stage-resolvent determinants are both \(1\).

Primitive boundary numerators:

\[
\begin{aligned}
P_0(x)={}&5x^7-161x^6+1097x^5-5000x^4+20000x^3\\
&-60000x^2+120000x-120000,
\end{aligned}
\]

\[
P_1(x)=x\left(5x^6-161x^5+1097x^4-5000x^3+20000x^2-60000x+120000\right),
\]

\[
\begin{aligned}
P_M(x)=x^5\big(&25x^9-3220x^8+147564x^7-3225872x^6+48214544x^5\\
&-576320000x^4+5721600000x^3-37719040000x^2\\
&+16752640000x+178790400000\big).
\end{aligned}
\]

| Label | Boundary | Certified rational bracket | Approximation | Mechanism |
|---|---|---|---:|---|
| \(\beta_-\) | \(M_R=0\) | \((3.068856548050,3.068856548051)\) | 3.068856548050381 | coherence-population block |
| \(\beta_+\) | \(a=1\) | \((4.384986320801,4.384986320802)\) | 4.384986320801944 | population ceiling |
| \(\delta_M\) | \(M_R=0\) | \((15.298767807585,15.298767807586)\) | 15.298767807585427 | coherence-population block |
| \(\delta_1\) | \(a=1\) | \((24.727754318760,24.727754318761)\) | 24.727754318760310 | population ceiling |
| \(\delta_0\) | \(a=0\) | \((24.727895030129,24.727895030130)\) | 24.727895030129697 | population floor |

| Open interval | \(\operatorname{sgn}a\) | \(\operatorname{sgn}(a-1)\) | \(\operatorname{sgn}M_R\) | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| \((0,\beta_-)\) | + | - | - | yes | no |
| \((\beta_-,\beta_+)\) | + | - | + | yes | yes |
| \((\beta_+,\delta_M)\) | + | + | + | yes | no |
| \((\delta_M,\delta_1)\) | + | + | - | yes | no |
| \((\delta_1,\delta_0)\) | + | - | - | yes | no |
| \((\delta_0,\infty)\) | - | - | - | yes | no |

Therefore,

\[
\mathcal C_{\mathrm{DP4}}=\{0\}\cup[\beta_-,\beta_+].
\]

The first exponential defect is

\[
(m,\eta_m)=\left(5,\frac{97}{120000}\right),
\]

and the leading margin coefficient is \(-97/128000\).

The later roots \(\delta_M\), \(\delta_1\), and \(\delta_0\) do not create a second positive admissible component: the population or margin inequality fails throughout every interval after \(\beta_+\).

## S1.8 Backward Euler

\[
R(s)=\frac{1}{1-s}.
\]

The stage-resolvent determinants are

\[
1+x
\qquad\text{and}\qquad
1+\frac{x}{2},
\]

which are positive for every \(x\ge 0\).

Primitive boundary numerators:

\[
P_0(x)=1,
\qquad
P_1(x)=x,
\qquad
P_M(x)=x^2.
\]

There are no positive boundary roots or stage-resolvent singularities.

| Open interval | \(\operatorname{sgn}a\) | \(\operatorname{sgn}(a-1)\) | \(\operatorname{sgn}M_R\) | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| \((0,\infty)\) | + | - | + | yes | yes |

Therefore,

\[
\mathcal C_{\mathrm{BE}}=[0,\infty).
\]

The first exponential defect is

\[
(m,\eta_m)=\left(2,\frac12\right),
\]

and the leading margin coefficient is \(1/4\).

## S1.9 Implicit midpoint

\[
R(s)=\frac{1+s/2}{1-s/2}.
\]

The stage-resolvent determinants are

\[
1+\frac{x}{2}
\qquad\text{and}\qquad
1+\frac{x}{4},
\]

which are positive for every \(x\ge 0\).

Primitive boundary numerators:

\[
P_0(x)=x-2,
\qquad
P_1(x)=x,
\qquad
P_M(x)=x^3.
\]

| Boundary | Certified root | Approximation | Mechanism |
|---|---:|---:|---|
| \(a=0\) | \(2\) | 2 | population floor |

| Open interval | \(\operatorname{sgn}a\) | \(\operatorname{sgn}(a-1)\) | \(\operatorname{sgn}M_R\) | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| \((0,2)\) | + | - | - | yes | no |
| \((2,\infty)\) | - | - | - | yes | no |

Therefore,

\[
\mathcal C_{\mathrm{IM}}=\{0\}.
\]

The first exponential defect is

\[
(m,\eta_m)=\left(3,\frac1{12}\right),
\]

and the leading margin coefficient is \(-1/16\).

## S1.10 Reproducibility record

Running

```bash
python src/reproduce_pla_results.py
```

reconstructs the stability functions from the exact tableaux, verifies the reference formulas and first defects, derives all boundary polynomials, certifies the 14 listed positive roots, reconstructs every interval sign chart, derives all eight admissible sets, and regenerates the data and figures.

Running

```bash
python src/reproduce_pla_results.py --check
```

performs the same calculations in a temporary directory and compares every tracked data and figure file byte-for-byte with the repository copy. The checks use explicit exceptions and remain active under optimized Python execution.
