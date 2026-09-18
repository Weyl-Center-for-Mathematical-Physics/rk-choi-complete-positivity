# Phase-covariant symbolic verification

This report independently reconstructs the exact Choi determinant, local buffer coefficients, coherent-precession boundary law including the cancellation thresholds, and backward-Euler domain used in the article.

- PASS: trace identity on E00
- PASS: trace identity on E11
- PASS: thermal Choi determinant identity
- PASS: exact phase-covariant margin
- PASS: exact margin is independent of coherent phase
- PASS: fixed-frequency bidirectional local coefficient
- PASS: order-two frequency cancellation
- PASS: Forward-Euler local coefficient
- PASS: frequency boundary law: Forward Euler
- PASS: frequency boundary law: Heun RK2
- PASS: frequency boundary law: SSPRK(3,3)
- PASS: frequency boundary law: Classical RK4
- PASS: frequency boundary law: Dormand-Prince 5 principal formula
- PASS: frequency boundary law: Dormand-Prince embedded 4 formula
- PASS: frequency boundary law: Backward Euler
- PASS: frequency boundary law: Implicit midpoint
- PASS: SSPRK(3,3) leading cancellation threshold
- PASS: classical RK4 leading cancellation threshold
- PASS: SSPRK(3,3) exact threshold factorization
- PASS: classical RK4 exact threshold factorization
- PASS: SSPRK(3,3) threshold is locally inward
- PASS: classical RK4 threshold is locally outward
- PASS: backward-Euler full zero-temperature domain
- PASS: backward-Euler margin identity

## Exact semigroup margin

```text
-(theta**2*exp(2*x) - 2*theta**2*exp(x) + theta**2 - theta*exp(2*x) + 2*theta*exp(x) - theta - exp(x) + exp(-2*u)*exp(x))*exp(-2*x)
```

## General fixed-frequency quadratic coefficient

```text
(8*r2*varpi**2 + 2*r2 - 4*theta**2 + 4*theta - 4*varpi**2 - 1)/4
```
