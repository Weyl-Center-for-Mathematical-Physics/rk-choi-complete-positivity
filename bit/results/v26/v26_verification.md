# Version 2.6 exact verification

- RK4 rotating cubic: `-32*q*y + 384*q + y**3 - 16*y**2 + 192*y - 1536`
- Discriminant: `16384*(q - 4)*(8*q**2 - 123*q + 348)`
- q_- = 3.738113180505105, varpi_- = 0.8273622514511260
- q = 4, varpi = sqrt(3)/2
- q_+ = 11.63688681949489, varpi_+ = 1.630712023894386
- Correct DP5 factor: G6 = `q*(q**2 - 18*q + 48)/32`
- Exact off-boundary cases: 5, each with three positive margin roots.
- Bell-state error rows: 9; the physicality defect never exceeds the total normalized-Choi trace error.
- Controller counterexample: exact detached interval at varpi=2 with a halving sequence that skips it.
