#!/usr/bin/env python3
from __future__ import annotations

"""Independent package-level executed-candidate contract checks."""

import json
from pathlib import Path
import numpy as np
import sympy as sp

from rk_choi_margin.candidates import (
    CandidateSpec,
    build_candidate_execution,
    certify_candidate,
    direct_candidate,
    equal_substep_candidate,
    phase_covariant_superoperator_numeric,
    verify_execution_provenance,
)
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import (
    _rk4_step_doubling_evaluation,
    certified_candidate_guard_step,
    run_adaptive_channel_benchmark,
)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'independent_validation'/'executed_candidate_contract_v29_report.json'
checks=[]

def add(name, condition, detail=''):
    checks.append({'name':name,'pass':bool(condition),'detail':detail})
    if not condition:
        raise AssertionError(f'{name}: {detail}')

# Projection and exact-binary recertification over a range of safety fractions.
for kind,proposed,substeps in [('direct',1.4,1),('equal-substeps',2.8,2)]:
    for exponent in (2,10,20,30):
        decision=certified_candidate_guard_step(
            method='rk4',candidate_kind=kind,substeps=substeps,proposed_x=proposed,
            theta=0,kappa=0,varpi=2,safety_fraction=sp.Rational(1,10**exponent)
        )
        add(f'{kind} projection pass 1e-{exponent}',decision.accepted_execution is not None and decision.cp_status is CPStatus.PASS)
        execution=decision.accepted_execution
        add(f'{kind} provenance 1e-{exponent}',verify_execution_provenance(execution))
        candidate=(direct_candidate('rk4',total_x=execution.executed_float,theta=0,kappa=0,varpi=2)
                   if kind=='direct' else equal_substep_candidate('rk4',total_x=execution.executed_float,substeps=2,theta=0,kappa=0,varpi=2))
        add(f'{kind} executed recertificate 1e-{exponent}',certify_candidate(candidate).status is CPStatus.PASS)

# Coarse/fine/extrapolated identity and distinct hashes.
for frame in ('lab','rotating'):
    for h in (0.5,1.0,2.0):
        ev=_rk4_step_doubling_evaluation(h=h,theta=0,kappa=0,varpi=2,frame=frame)
        add(f'{frame} coarse provenance h={h}',verify_execution_provenance(ev.coarse_execution))
        add(f'{frame} fine provenance h={h}',verify_execution_provenance(ev.fine_execution))
        add(f'{frame} fine matrix identity h={h}',np.linalg.norm(phase_covariant_superoperator_numeric(ev.bundle.fine)-ev.fine_matrix)<1e-13)
        if ev.extrapolated_execution is not None:
            add(f'{frame} extrapolated provenance h={h}',verify_execution_provenance(ev.extrapolated_execution))
            add(f'{frame} distinct hashes h={h}',len({ev.coarse_execution.provenance_hash,ev.fine_execution.provenance_hash,ev.extrapolated_execution.provenance_hash})==3)

# Candidate-guard benchmark must accumulate only provenance-verified maps.
for final_time,initial_h,tol in ((2.0,1.0,1e-2),(2.8,2.8,10.0),(5.6,5.6,2.0)):
    result=run_adaptive_channel_benchmark(final_time=final_time,initial_h=initial_h,tolerance=tol,theta=0,kappa=0,varpi=2 if final_time<5 else 1,policy='candidate_guard',fallback='rotating_frame')
    add(f'benchmark provenance {final_time}',all(row.provenance_verified for row in result.records))
    add(f'benchmark physicality {final_time}',all(row.min_choi_eigenvalue>-1e-12 for row in result.records))
    add(f'benchmark identity {final_time}',all(row.candidate_kind==row.certified_kind or row.action.startswith('fallback') for row in result.records))

report={'checks':checks,'passed':sum(row['pass'] for row in checks),'total':len(checks)}
OUT.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(f"Executed-candidate contract v2.9: {report['passed']}/{report['total']} checks passed.")
print(f'Wrote {OUT}')
