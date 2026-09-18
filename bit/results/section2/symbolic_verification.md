# Section 2 symbolic verification

- PASS: L(E00)=0
- PASS: L(E11-E00)=-gamma(E11-E00)
- PASS: L(E01)=-(gamma/2+gamma_phi)E01
- PASS: L(E10)=-(gamma/2+gamma_phi)E10
- PASS: trace_row*T=trace_row
- PASS: T(E00)=E00
- PASS: T(E11)=(1-a)E00+aE11

The exact numerical superoperator in vectorization order (E00,E01,E10,E11) is

```text
Matrix([[1, 0, 0, 1 - a], [0, c, 0, 0], [0, 0, c, 0], [0, 0, 0, a]])
```
