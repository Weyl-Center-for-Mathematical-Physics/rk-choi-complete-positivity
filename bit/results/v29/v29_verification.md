# v2.9 candidate-map verification

## Candidate-map topology
- full RK4 interval at varpi=2: `[{'lower': '15274297615272439/20508402243624044', 'upper': '25714392766941956/20242219822187449', 'lower_decimal': 0.7447824278959195, 'upper_decimal': 1.270334626973889, 'left_sources': ['F'], 'right_sources': ['F']}]`
- two-half-step interval: `[{'lower': '15274297615272439/10254201121812022', 'upper': '125112718618773356/49244000738926945', 'lower_decimal': 1.489564855791839, 'upper_decimal': 2.540669253947778, 'left_sources': ['F'], 'right_sources': ['F']}]`
- the intervals are disjoint: **True**

## Richardson extrapolation
- stability function: `s**8/138240 + s**7/8640 + s**6/864 + s**5/120 + s**4/24 + s**3/6 + s**2/2 + s + 1`
- first exponential defect: m=6, eta=-1/4320
- exact margin leading coefficient: `-31/138240`
- global remote component: `[{'lower': '82720409432394497/13707859973582984', 'upper': '858388633779760739/132895441092770412', 'lower_decimal': 6.03452395864917, 'upper_decimal': 6.459127767825721, 'left_sources': ['F'], 'right_sources': ['1-a']}]`

## Executed-candidate controller
- projection component searches: `1`
- every accumulated map in the projection benchmark has verified provenance
