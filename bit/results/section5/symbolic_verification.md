# Section 5 verification

## Exact repeated-step factorization

- N=1: verified exactly
- N=2: verified exactly
- N=3: verified exactly
- N=5: verified exactly
- N=8: verified exactly

## Fixed-horizon asymptotics at X=1

### Forward Euler
- m = 2
- K = -0.25
- observed M_N/x^(m-1) = -0.0919716145099619269
- predicted = -0.0919698602928605804
- observed lambda_-/x^(m-1) = -0.0672368128431265676
- predicted = -0.0672353553424987802

### SSPRK(3,3)
- m = 4
- K = -0.0364583333333333333
- observed M_N/x^(m-1) = -0.0134126221344932752
- predicted = -0.0134122712927088346
- observed lambda_-/x^(m-1) = -0.00980541247334399702
- predicted = -0.00980515598744773878

### Dormand-Prince embedded 4 formula
- m = 5
- K = -0.0007578125
- observed M_N/x^(m-1) = -0.000278792944093613479
- predicted = -0.000278783639012733634
- observed lambda_-/x^(m-1) = -0.000203813973441151484
- predicted = -0.000203807170881949427

## Logical-separation checks

- Forward Euler at x=1: both scalar multipliers satisfy absolute stability, but M=-1/4.
- SSPRK(3,3) admissible set: [(0.0, 0.0)]
- Classical RK4 admissible set: [(0.0, 2.7852935634053)]
- Embedded DP4 at x=4: a=847/1875, c=91/750, M=245819/562500 > 0.
