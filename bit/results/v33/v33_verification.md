# v3.3 external-review closeout verification

## Candidate-map topology
- full RK4 interval at varpi=2: `[{'lower': '1019017290566537485/1368208019414954846', 'upper': '2006012581244003614/1579121389474047577', 'lower_decimal': 0.7447824278959195, 'upper_decimal': 1.270334626973889, 'left_sources': ['F'], 'right_sources': ['F']}]`
- two-half-step interval: `[{'lower': '3660348255086624578/2457327212611240667', 'upper': '4012025162488007228/1579121389474047577', 'lower_decimal': 1.489564855791839, 'upper_decimal': 2.540669253947778, 'left_sources': ['F'], 'right_sources': ['F']}]`
- nontrivial positive components are disjoint

## Richardson extrapolation
- first defect: m=6, eta=-1/4320
- global remote component: `[{'lower': '544428513983366551531/90218966353269235844', 'upper': '2637367531516119310/408316358851639967', 'lower_decimal': 6.03452395864917, 'upper_decimal': 6.459127767825721, 'left_sources': ['F'], 'right_sources': ['1-a']}]`

## Adaptive evidence
- boundary guard/fallback action totals: `{'accepted_steps': 15, 'direct_accepts': 0, 'projected_accepts': 0, 'rotating_accepts': 15, 'strang_accepts': 0, 'exact_fallback_accepts': 0, 'component_searches': 0, 'fallback_steps': 15}`
- buffered useful case: direct=3, projected=2, error=0.0002853084392580805
- near-saddle useful case: projected=1, searches=1, error=0.00021027592790007218
- all accumulated records have verified provenance

## Comparator scope
- noncommuting dense-qubit RK4/Strang data are a common-representation reference only; no universal scalability claim is made.
