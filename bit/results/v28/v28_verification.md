# v2.8 candidate-map verification

## Exact binary input semantics
- binary float status at `1.270334626973889`: **FAIL**
- exact decimal status at `1.270334626973889`: **PASS**

## Candidate-map mismatch at varpi=2
- H=1: coarse PASS, two-half FAIL, Richardson FAIL
- H=2: coarse FAIL, two-half PASS, Richardson FAIL

## Richardson extrapolation
- exact leading coefficient: `-31/138240`
- the degree-ten factor has no positive root before the classical RK4 population ceiling
- therefore the negative-weight extrapolate is non-CPTP while its coarse and fine constituents are CPTP on their common positive interval

## Adaptive benchmark
See `adaptive_benchmark_v28.csv` for the complete tolerance sweep.
