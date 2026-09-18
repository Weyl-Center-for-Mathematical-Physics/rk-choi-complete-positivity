# Supplemental Material: Exact boundary polynomials, root isolation, and sign charts

This certificate reconstructs the pure-amplitude-damping CPTP-admissible set for each audited Runge-Kutta formula from exact rational data. For every method it records the primitive numerator polynomials associated with `a(x)=R(-x)`, `a(x)-1`, and `M_R(x,0)=a(x)-R(-x/2)^2`; the operational stage-resolvent factors at the full and half step; all nonnegative real roots; rational isolating intervals; and the exact sign on every intervening open interval. Finite decimal endpoints in the brackets are exact rational numbers. Root counts were certified by Sturm sequences (`Poly.count_roots`) and signs were evaluated at exact rational sample points.

The one-step map is CPTP exactly when it lies in the operational RK domain and satisfies `a(x) <= 1` and `M_R(x,0) >= 0`. The second inequality already forces `a(x) >= 0`.

## S1.1 Summary

| Method | Stage domain on x >= 0 | Certified admissible set |
|---|---|---|
| Forward Euler | all x >= 0 | `$\{0\}$` |
| Heun RK2 | all x >= 0 | `$[0,2]$` |
| SSPRK(3,3) | all x >= 0 | `$\{0\}$` |
| Classical RK4 | all x >= 0 | `$[0,\alpha_4]$` |
| Dormand-Prince 5 principal formula | all x >= 0 | `$[0,\alpha_5]$` |
| Dormand-Prince embedded 4 formula | all x >= 0 | `$\{0\}\cup[\beta_-,\beta_+]$` |
| Backward Euler | all x >= 0 | `$[0,\infty)$` |
| Implicit midpoint | all x >= 0 | `$\{0\}$` |

## S1.2 Forward Euler

Stability function: `$s + 1`.
Exact population factor: `$1 - x`.
Exact Choi margin: `$-x**2/4`.

Stage-resolvent factors on the full and half steps are `1` and `1`, respectively.

Primitive boundary numerators (overall nonzero scalar factors are recorded in the machine-readable JSON):

- `P_0(x) = x - 1` for `a(x)=0`; square-free part `x - 1`.
- `P_1(x) = x` for `a(x)=1`; square-free part `x`.
- `P_M(x) = x**2` for `M_R(x,0)=0`; square-free part `x`.

| Label | Certified rational bracket | Approximation | Source / active mechanism |
|---|---|---:|---|
| - | `1` | 1 | population floor (a=0) |

| Open interval | sgn a | sgn(a-1) | sgn M | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| `(0, 1)` | + | - | - | yes | no |
| `(1, infinity)` | - | - | - | yes | no |

**Certified result:** `$\{0\}`. The first exponential defect is `(m, eta_m)=(2, -1/2)`, with predicted local margin coefficient `-1/4`.

## S1.3 Heun RK2

Stability function: `$s**2/2 + s + 1`.
Exact population factor: `$(x**2 - 2*x + 2)/2`.
Exact Choi margin: `$-x**3*(x - 8)/64`.

Stage-resolvent factors on the full and half steps are `1` and `1`, respectively.

Primitive boundary numerators (overall nonzero scalar factors are recorded in the machine-readable JSON):

- `P_0(x) = x**2 - 2*x + 2` for `a(x)=0`; square-free part `x**2 - 2*x + 2`.
- `P_1(x) = x**2 - 2*x` for `a(x)=1`; square-free part `x**2 - 2*x`.
- `P_M(x) = x**4 - 8*x**3` for `M_R(x,0)=0`; square-free part `x**2 - 8*x`.

| Label | Certified rational bracket | Approximation | Source / active mechanism |
|---|---|---:|---|
| - | `2` | 2 | population ceiling (a=1) |
| - | `8` | 8 | coherence-population block (M=0) |

| Open interval | sgn a | sgn(a-1) | sgn M | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| `(0, 2)` | + | - | + | yes | yes |
| `(2, 8)` | + | + | + | yes | no |
| `(8, infinity)` | + | + | - | yes | no |

**Certified result:** `$[0,2]`. The first exponential defect is `(m, eta_m)=(3, -1/6)`, with predicted local margin coefficient `1/8`.

## S1.4 SSPRK(3,3)

Stability function: `$s**3/6 + s**2/2 + s + 1`.
Exact population factor: `$-(x**3 - 3*x**2 + 6*x - 6)/6`.
Exact Choi margin: `$-x**4*(x**2 - 12*x + 84)/2304`.

Stage-resolvent factors on the full and half steps are `1` and `1`, respectively.

Primitive boundary numerators (overall nonzero scalar factors are recorded in the machine-readable JSON):

- `P_0(x) = x**3 - 3*x**2 + 6*x - 6` for `a(x)=0`; square-free part `x**3 - 3*x**2 + 6*x - 6`.
- `P_1(x) = x**3 - 3*x**2 + 6*x` for `a(x)=1`; square-free part `x**3 - 3*x**2 + 6*x`.
- `P_M(x) = x**6 - 12*x**5 + 84*x**4` for `M_R(x,0)=0`; square-free part `x**3 - 12*x**2 + 84*x`.

| Label | Certified rational bracket | Approximation | Source / active mechanism |
|---|---|---:|---|
| - | `(1.596071637983321, 1.596071637983322)` | 1.596071637983322 | population floor (a=0) |

| Open interval | sgn a | sgn(a-1) | sgn M | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| `(0, 1.596071637983322)` | + | - | - | yes | no |
| `(1.596071637983322, infinity)` | - | - | - | yes | no |

**Certified result:** `$\{0\}`. The first exponential defect is `(m, eta_m)=(4, -1/24)`, with predicted local margin coefficient `-7/192`.

## S1.5 Classical RK4

Stability function: `$s**4/24 + s**3/6 + s**2/2 + s + 1`.
Exact population factor: `$(x**4 - 4*x**3 + 12*x**2 - 24*x + 24)/24`.
Exact Choi margin: `$-x**5*(x**3 - 16*x**2 + 160*x - 1152)/147456`.

Stage-resolvent factors on the full and half steps are `1` and `1`, respectively.

Primitive boundary numerators (overall nonzero scalar factors are recorded in the machine-readable JSON):

- `P_0(x) = x**4 - 4*x**3 + 12*x**2 - 24*x + 24` for `a(x)=0`; square-free part `x**4 - 4*x**3 + 12*x**2 - 24*x + 24`.
- `P_1(x) = x**4 - 4*x**3 + 12*x**2 - 24*x` for `a(x)=1`; square-free part `x**4 - 4*x**3 + 12*x**2 - 24*x`.
- `P_M(x) = x**8 - 16*x**7 + 160*x**6 - 1152*x**5` for `M_R(x,0)=0`; square-free part `x**4 - 16*x**3 + 160*x**2 - 1152*x`.

| Label | Certified rational bracket | Approximation | Source / active mechanism |
|---|---|---:|---|
| \alpha_4 | `(2.785293563405281, 2.785293563405282)` | 2.785293563405282 | population ceiling (a=1) |
| \mu_4 | `(10.982425466293273, 10.982425466293274)` | 10.982425466293273 | coherence-population block (M=0) |

| Open interval | sgn a | sgn(a-1) | sgn M | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| `(0, 2.785293563405282)` | + | - | + | yes | yes |
| `(2.785293563405282, 10.982425466293273)` | + | + | + | yes | no |
| `(10.982425466293273, infinity)` | + | + | - | yes | no |

**Certified result:** `$[0,\alpha_4]`. The first exponential defect is `(m, eta_m)=(5, -1/120)`, with predicted local margin coefficient `1/128`.

## S1.6 Dormand-Prince 5 principal formula

Stability function: `$s**6/600 + s**5/120 + s**4/24 + s**3/6 + s**2/2 + s + 1`.
Exact population factor: `$(x**6 - 5*x**5 + 25*x**4 - 100*x**3 + 300*x**2 - 600*x + 600)/600`.
Exact Choi margin: `$-x**6*(x**6 - 20*x**5 + 300*x**4 - 3600*x**3 + 35600*x**2 - 294400*x - 396800)/1474560000`.

Stage-resolvent factors on the full and half steps are `1` and `1`, respectively.

Primitive boundary numerators (overall nonzero scalar factors are recorded in the machine-readable JSON):

- `P_0(x) = x**6 - 5*x**5 + 25*x**4 - 100*x**3 + 300*x**2 - 600*x + 600` for `a(x)=0`; square-free part `x**6 - 5*x**5 + 25*x**4 - 100*x**3 + 300*x**2 - 600*x + 600`.
- `P_1(x) = x**6 - 5*x**5 + 25*x**4 - 100*x**3 + 300*x**2 - 600*x` for `a(x)=1`; square-free part `x**6 - 5*x**5 + 25*x**4 - 100*x**3 + 300*x**2 - 600*x`.
- `P_M(x) = x**12 - 20*x**11 + 300*x**10 - 3600*x**9 + 35600*x**8 - 294400*x**7 - 396800*x**6` for `M_R(x,0)=0`; square-free part `x**7 - 20*x**6 + 300*x**5 - 3600*x**4 + 35600*x**3 - 294400*x**2 - 396800*x`.

| Label | Certified rational bracket | Approximation | Source / active mechanism |
|---|---|---:|---|
| \alpha_5 | `(3.306567892634946, 3.306567892634947)` | 3.306567892634947 | population ceiling (a=1) |
| \mu_5 | `(13.174093462800933, 13.174093462800934)` | 13.174093462800933 | coherence-population block (M=0) |

| Open interval | sgn a | sgn(a-1) | sgn M | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| `(0, 3.306567892634947)` | + | - | + | yes | yes |
| `(3.306567892634947, 13.174093462800933)` | + | + | + | yes | no |
| `(13.174093462800933, infinity)` | + | + | - | yes | no |

**Certified result:** `$[0,\alpha_5]`. The first exponential defect is `(m, eta_m)=(6, 1/3600)`, with predicted local margin coefficient `31/115200`.

## S1.7 Dormand-Prince embedded 4 formula

Stability function: `$s**7/24000 + 161*s**6/120000 + 1097*s**5/120000 + s**4/24 + s**3/6 + s**2/2 + s + 1`.
Exact population factor: `$-(5*x**7 - 161*x**6 + 1097*x**5 - 5000*x**4 + 20000*x**3 - 60000*x**2 + 120000*x - 120000)/120000`.
Exact Choi margin: `$-x**5*(25*x**9 - 3220*x**8 + 147564*x**7 - 3225872*x**6 + 48214544*x**5 - 576320000*x**4 + 5721600000*x**3 - 37719040000*x**2 + 16752640000*x + 178790400000)/235929600000000`.

Stage-resolvent factors on the full and half steps are `1` and `1`, respectively.

Primitive boundary numerators (overall nonzero scalar factors are recorded in the machine-readable JSON):

- `P_0(x) = 5*x**7 - 161*x**6 + 1097*x**5 - 5000*x**4 + 20000*x**3 - 60000*x**2 + 120000*x - 120000` for `a(x)=0`; square-free part `5*x**7 - 161*x**6 + 1097*x**5 - 5000*x**4 + 20000*x**3 - 60000*x**2 + 120000*x - 120000`.
- `P_1(x) = 5*x**7 - 161*x**6 + 1097*x**5 - 5000*x**4 + 20000*x**3 - 60000*x**2 + 120000*x` for `a(x)=1`; square-free part `5*x**7 - 161*x**6 + 1097*x**5 - 5000*x**4 + 20000*x**3 - 60000*x**2 + 120000*x`.
- `P_M(x) = 25*x**14 - 3220*x**13 + 147564*x**12 - 3225872*x**11 + 48214544*x**10 - 576320000*x**9 + 5721600000*x**8 - 37719040000*x**7 + 16752640000*x**6 + 178790400000*x**5` for `M_R(x,0)=0`; square-free part `25*x**10 - 3220*x**9 + 147564*x**8 - 3225872*x**7 + 48214544*x**6 - 576320000*x**5 + 5721600000*x**4 - 37719040000*x**3 + 16752640000*x**2 + 178790400000*x`.

| Label | Certified rational bracket | Approximation | Source / active mechanism |
|---|---|---:|---|
| \beta_- | `(3.06885654805038, 3.068856548050381)` | 3.068856548050381 | coherence-population block (M=0) |
| \beta_+ | `(4.384986320801944, 4.384986320801945)` | 4.384986320801944 | population ceiling (a=1) |
| \delta_M | `(15.298767807585426, 15.298767807585427)` | 15.298767807585427 | coherence-population block (M=0) |
| \delta_1 | `(24.72775431876031, 24.727754318760311)` | 24.72775431876031 | population ceiling (a=1) |
| \delta_0 | `(24.727895030129697, 24.727895030129698)` | 24.727895030129697 | population floor (a=0) |

| Open interval | sgn a | sgn(a-1) | sgn M | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| `(0, 3.068856548050381)` | + | - | - | yes | no |
| `(3.068856548050381, 4.384986320801944)` | + | - | + | yes | yes |
| `(4.384986320801944, 15.298767807585427)` | + | + | + | yes | no |
| `(15.298767807585427, 24.72775431876031)` | + | + | - | yes | no |
| `(24.72775431876031, 24.727895030129697)` | + | - | - | yes | no |
| `(24.727895030129697, infinity)` | - | - | - | yes | no |

**Certified result:** `$\{0\}\cup[\beta_-,\beta_+]`. The first exponential defect is `(m, eta_m)=(5, 97/120000)`, with predicted local margin coefficient `-97/128000`.

## S1.8 Backward Euler

Stability function: `$-1/(s - 1)`.
Exact population factor: `$1/(x + 1)`.
Exact Choi margin: `$x**2/((x + 1)*(x + 2)**2)`.

Stage-resolvent factors on the full and half steps are `x + 1` and `(x + 2)/2`, respectively.

Primitive boundary numerators (overall nonzero scalar factors are recorded in the machine-readable JSON):

- `P_0(x) = 1` for `a(x)=0`; square-free part `1`.
- `P_1(x) = x` for `a(x)=1`; square-free part `x`.
- `P_M(x) = x**2` for `M_R(x,0)=0`; square-free part `x`.

No positive real boundary roots or stage-resolvent singularities occur.

| Open interval | sgn a | sgn(a-1) | sgn M | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| `(0, infinity)` | + | - | + | yes | yes |

**Certified result:** `$[0,\infty)`. The first exponential defect is `(m, eta_m)=(2, 1/2)`, with predicted local margin coefficient `1/4`.

## S1.9 Implicit midpoint

Stability function: `$-(s + 2)/(s - 2)`.
Exact population factor: `$-(x - 2)/(x + 2)`.
Exact Choi margin: `$-2*x**3/((x + 2)*(x + 4)**2)`.

Stage-resolvent factors on the full and half steps are `(x + 2)/2` and `(x + 4)/4`, respectively.

Primitive boundary numerators (overall nonzero scalar factors are recorded in the machine-readable JSON):

- `P_0(x) = x - 2` for `a(x)=0`; square-free part `x - 2`.
- `P_1(x) = x` for `a(x)=1`; square-free part `x`.
- `P_M(x) = x**3` for `M_R(x,0)=0`; square-free part `x`.

| Label | Certified rational bracket | Approximation | Source / active mechanism |
|---|---|---:|---|
| - | `2` | 2 | population floor (a=0) |

| Open interval | sgn a | sgn(a-1) | sgn M | RK domain | CPTP |
|---|---:|---:|---:|:---:|:---:|
| `(0, 2)` | + | - | - | yes | no |
| `(2, infinity)` | - | - | - | yes | no |

**Certified result:** `$\{0\}`. The first exponential defect is `(m, eta_m)=(3, 1/12)`, with predicted local margin coefficient `-1/16`.

## S1.10 Dormand-Prince embedded fourth-order no-omission check

The embedded fourth-order formula has, in increasing positive order, the first margin zero `beta_-`, the first population-ceiling zero `beta_+`, a second margin zero `delta_M`, a second population-ceiling zero `delta_1`, and finally the population-floor zero `delta_0`. The sign chart shows: `M<0` on `(0,beta_-)`; both inequalities hold only on `(beta_-,beta_+)`; `a>1` on `(beta_+,delta_1)`; `M<0` after `delta_M`; and `a<0` after `delta_0`. Hence no later root creates another admissible component.

## S1.11 Reproducibility statement

The companion JSON contains exact expressions, rational brackets, multiplicities, exact rational sign samples, and interval truth values. The generator script reconstructs this document directly from the exact Butcher tableaux in `src/rk_choi_margin/methods.py`.
