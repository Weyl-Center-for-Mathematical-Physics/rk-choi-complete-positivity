# Extrapolation corollary evidence

All checks pass: True

| Method | p | n | eta_{p+1} | eta_{p+2} | defect order | eta^E | predicted | orientation at varpi=0 |
|---|---|---|---|---|---|---|---|---|
| Forward Euler | 1 | 2 | -1/2 | -1/6 | 3 | -1/6 | -1/6 | CP |
| Forward Euler | 1 | 3 | -1/2 | -1/6 | 3 | -1/9 | -1/9 | CP |
| Forward Euler | 1 | 4 | -1/2 | -1/6 | 3 | -1/12 | -1/12 | CP |
| Heun RK2 | 2 | 2 | -1/6 | -1/24 | 4 | -1/48 | -1/48 | non-CP |
| Heun RK2 | 2 | 3 | -1/6 | -1/24 | 4 | -1/96 | -1/96 | non-CP |
| Heun RK2 | 2 | 4 | -1/6 | -1/24 | 4 | -1/160 | -1/160 | non-CP |
| SSPRK(3,3) | 3 | 2 | -1/24 | -1/120 | 5 | -1/420 | -1/420 | CP |
| SSPRK(3,3) | 3 | 3 | -1/24 | -1/120 | 5 | -1/1170 | -1/1170 | CP |
| SSPRK(3,3) | 3 | 4 | -1/24 | -1/120 | 5 | -1/2520 | -1/2520 | CP |
| Classical RK4 | 4 | 2 | -1/120 | -1/720 | 6 | -1/4320 | -1/4320 | non-CP |
| Classical RK4 | 4 | 3 | -1/120 | -1/720 | 6 | -1/17280 | -1/17280 | non-CP |
| Classical RK4 | 4 | 4 | -1/120 | -1/720 | 6 | -1/48960 | -1/48960 | non-CP |
| Dormand-Prince 5 principal formula | 5 | 2 | 1/3600 | -1/5040 | 7 | 1/130200 | 1/130200 | non-CP |
| Dormand-Prince 5 principal formula | 5 | 3 | 1/3600 | -1/5040 | 7 | 1/762300 | 1/762300 | non-CP |
| Dormand-Prince 5 principal formula | 5 | 4 | 1/3600 | -1/5040 | 7 | 1/2864400 | 1/2864400 | non-CP |
| Dormand-Prince embedded 4 formula | 4 | 2 | 97/120000 | -17/360000 | 6 | 77/2700000 | 77/2700000 | CP |
| Dormand-Prince embedded 4 formula | 4 | 3 | 97/120000 | -17/360000 | 6 | 77/10800000 | 77/10800000 | CP |
| Dormand-Prince embedded 4 formula | 4 | 4 | 97/120000 | -17/360000 | 6 | 77/30600000 | 77/30600000 | CP |
| Backward Euler | 1 | 2 | 1/2 | 5/6 | 3 | -1/6 | -1/6 | CP |
| Backward Euler | 1 | 3 | 1/2 | 5/6 | 3 | -1/9 | -1/9 | CP |
| Backward Euler | 1 | 4 | 1/2 | 5/6 | 3 | -1/12 | -1/12 | CP |
| Implicit midpoint | 2 | 2 | 1/12 | 1/12 | 5 | -1/320 | 0 | degenerate |
| Implicit midpoint | 2 | 3 | 1/12 | 1/12 | 5 | -1/720 | 0 | degenerate |
| Implicit midpoint | 2 | 4 | 1/12 | 1/12 | 5 | -1/1280 | 0 | degenerate |
