# Version 2.7 certified verification

## Blocking guard defect

- x=0.004: exact M=-5.86498858453777777777777777778e-13; legacy tolerance accepts=True; certified=FAIL.
- x=0.001: exact M=-5.74830413652350531684027777778e-16; legacy tolerance accepts=True; certified=FAIL.
- x=0.0001: exact M=-5.75451839287675788031684027778e-21; legacy tolerance accepts=True; certified=FAIL.
- x=0.00001: exact M=-5.75513934408496660163031684028e-26; legacy tolerance accepts=True; certified=FAIL.

## Halving threshold

- resultant: `6291456*(q - 4)*(72*q**3 - 2071*q**2 + 12552*q - 21744)`
- varpi_1/2: 2.24825460451498048047320960260

## Noncommuting persistence

- Omega_x/Gamma=0: components [('0', '0.2620333331608279'), ('0.7383646272134431', '1.269753570826389')]
- Omega_x/Gamma=0.05: components [('0', '0.2624189364007533'), ('0.7377443064252767', '1.269432714005483')]
- Omega_x/Gamma=0.1: components [('0', '0.2635769098745156'), ('0.7358838130069857', '1.268472478716630')]
- Omega_x/Gamma=0.2: components [('0', '0.2682224579252424'), ('0.7284582232649768', '1.264665614699492')]

## Adaptive benchmark summary

- boundary / error_only: accepted=4, min lambda=-0.00158394, global Bell/Choi error=0.00155257.
- boundary / certified_guard: accepted=4, min lambda=-6.19018e-28, global Bell/Choi error=1.2503e-05.
- boundary / rotating_frame: accepted=2, min lambda=0, global Bell/Choi error=0.00034661.
- buffered / error_only: accepted=4, min lambda=-0.00094355, global Bell/Choi error=0.00155012.
- buffered / certified_guard: accepted=4, min lambda=7.8429e-05, global Bell/Choi error=0.000508925.
- buffered / rotating_frame: accepted=2, min lambda=0.00045116, global Bell/Choi error=0.000346291.
- near_saddle / error_only: accepted=3, min lambda=-0.00118494, global Bell/Choi error=0.000883554.
- near_saddle / certified_guard: accepted=3, min lambda=-1.2398e-27, global Bell/Choi error=2.70249e-05.
- near_saddle / rotating_frame: accepted=2, min lambda=-2.28336e-28, global Bell/Choi error=0.00034661.
