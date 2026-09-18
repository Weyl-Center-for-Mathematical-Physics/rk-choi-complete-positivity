# Section 3 symbolic verification

The checks below independently reconstruct the exact Choi reduction and the local asymptotic laws used in Section 3.

- PASS: Choi characteristic polynomial factors into the scalar and block modes
- PASS: rationalized lower-eigenvalue identity is exact
- PASS: block discriminant vanishes at the stated cusp point
- PASS: RK2 population-exit counterexample
- PASS: two-rate local gradient is (0,2)
- PASS: pure-dephasing slice is one minus the squared stability factor
- PASS: backward Euler margin is nonnegative on the full quadrant
- PASS: signed first-defect coefficient for m=2
- PASS: signed first-defect coefficient for m=3
- PASS: signed first-defect coefficient for m=4
- PASS: signed first-defect coefficient for m=5
- PASS: signed first-defect coefficient for m=6
- PASS: signed first-defect coefficient for m=7
- PASS: signed first-defect coefficient for m=8

## Exact Choi matrix

```text
Matrix([[1, 0, 0, c], [0, 0, 0, 0], [0, 0, 1 - a, 0], [c, 0, 0, a]])
```

## Backward-Euler two-rate margin

```text
(4*u**2 + 4*u*z + 8*u + z**2)/((z + 1)*(2*u + z + 2)**2)
```
