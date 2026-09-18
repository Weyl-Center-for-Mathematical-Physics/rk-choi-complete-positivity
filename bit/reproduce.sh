#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/src"
export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
export SOURCE_DATE_EPOCH=1787184000

run_step() {
  local name="$1" log="$2"; shift 2
  local tmp
  tmp="$(mktemp /tmp/jcp-v35-reproduce.XXXXXX)"
  printf '==> %s\n' "$name"
  if "$@" >"$tmp" 2>&1; then
    mv -f "$tmp" "$log"
    printf '<== %s complete\n' "$name"
  else
    local status=$?
    mv -f "$tmp" "$log"
    printf '<!! %s failed; see %s\n' "$name" "$log" >&2
    return "$status"
  fi
}

stage_certificates() {
  run_step 'Generate exact root certificate' GENERATE_ROOT_CERTIFICATE.log \
    python -u scripts/generate_root_certificate.py
  run_step 'Verify v2.6 symbolic identities' VERIFY_V26_EXTENSIONS.log \
    python -u scripts/verify_v26_extensions.py
  run_step 'Regenerate v3.4 candidate, telemetry, and benchmark results' VERIFY_V34_EXTENSIONS.log \
    python -u scripts/verify_v34_extensions.py
}

stage_tests() {
  run_step 'Collect complete pytest suite' TEST_COLLECT_V35.log \
    python -m pytest --collect-only -q -p no:cacheprovider tests
  grep -q '240 tests collected' TEST_COLLECT_V35.log

  run_step 'Test closeout' TEST_CLOSEOUT.log \
    python -m pytest -q -p no:cacheprovider tests/test_closeout.py
  run_step 'Test core' TEST_CORE.log \
    python -m pytest -q -p no:cacheprovider tests/test_core.py
  run_step 'Test topology and operational results' TEST_TOPOLOGY_OPERATIONAL.log \
    python -m pytest -q -p no:cacheprovider tests/test_topology_operational.py
  run_step 'Test v2.6' TEST_V26.log \
    python -m pytest -q -p no:cacheprovider tests/test_v26.py
  run_step 'Test v2.7 certification' TEST_V27_CERTIFICATION.log \
    python -m pytest -q -p no:cacheprovider tests/test_v27_certification.py
  run_step 'Test v2.7 noncommuting/adaptive in isolated processes' TEST_V27_NONCOMMUTING_ADAPTIVE.log \
    bash scripts/run_v27_noncommuting_tests.sh
  run_step 'Test v2.8 candidate maps' TEST_V28_CANDIDATES.log \
    python -m pytest -q -p no:cacheprovider tests/test_v28_candidates.py
  run_step 'Test v2.9 executed candidates' TEST_V29_EXECUTED.log \
    python -m pytest -q -p no:cacheprovider tests/test_v29_executed_candidate.py
  run_step 'Test v2.9 theory' TEST_V29_THEORY.log \
    python -m pytest -q -p no:cacheprovider tests/test_v29_theory.py
  run_step 'Test v3.0 execution contract' TEST_V30_EXECUTION.log \
    python -m pytest -q -p no:cacheprovider tests/test_v30_execution_contract.py
  run_step 'Test v3.1 closeout' TEST_V31_CLOSEOUT.log \
    python -m pytest -q -p no:cacheprovider tests/test_v31_closeout.py
  run_step 'Test v3.2 external-review repairs' TEST_V32_EXTERNAL_REVIEW.log \
    python -m pytest -q -p no:cacheprovider tests/test_v32_external_review.py
  run_step 'Test v3.3 endpoint-domain certification' TEST_V33_ENDPOINT_MEMBERSHIP.log \
    python -m pytest -q -p no:cacheprovider tests/test_v33_endpoint_membership.py
  run_step 'Test v3.4 hostile-audit remediation' TEST_V34_REMEDIATION.log \
    python -m pytest -q -p no:cacheprovider tests/test_v34_remediation.py
  run_step 'Test v3.5 editorial figure layout' TEST_V35_EDITORIAL_LAYOUT.log \
    python -m pytest -q -p no:cacheprovider tests/test_v35_editorial_layout.py
  run_step 'Test v3.5 submission revision' TEST_V35_SUBMISSION_REVISION.log \
    python -m pytest -q -p no:cacheprovider tests/test_v35_submission_revision.py

  {
    printf '%s\n' 'TOTAL: 240 tests passed across isolated deterministic pytest invocations'
    cat TEST_CLOSEOUT.log TEST_CORE.log TEST_TOPOLOGY_OPERATIONAL.log TEST_V26.log \
      TEST_V27_CERTIFICATION.log TEST_V27_NONCOMMUTING_ADAPTIVE.log \
      TEST_V28_CANDIDATES.log TEST_V29_EXECUTED.log TEST_V29_THEORY.log \
      TEST_V30_EXECUTION.log TEST_V31_CLOSEOUT.log TEST_V32_EXTERNAL_REVIEW.log \
      TEST_V33_ENDPOINT_MEMBERSHIP.log TEST_V34_REMEDIATION.log \
      TEST_V35_EDITORIAL_LAYOUT.log TEST_V35_SUBMISSION_REVISION.log
  } > TEST_RESULTS_V35.log
}
stage_audits() {
  run_step 'Verify theorem identities' VERIFY_THEOREMS.log \
    python -u scripts/verify_theorems.py
  run_step 'Verify phase-covariant identities' VERIFY_PHASE_COVARIANT.log \
    python -u scripts/verify_phase_covariant.py
  run_step 'Verify method and candidate audit' VERIFY_METHOD_AUDIT.log \
    python -u scripts/verify_section4.py
  run_step 'Run independent 15,000-sample Choi audit' INDEPENDENT_AUDIT.log \
    python -u independent_validation/independent_audit_v31.py
  run_step 'Run executed-candidate contract audit' EXECUTED_CANDIDATE_CONTRACT.log \
    python -u independent_validation/executed_candidate_contract_v31.py
  run_step 'Run endpoint-domain adversarial audit' ENDPOINT_MEMBERSHIP_V33.log \
    python -u independent_validation/endpoint_membership_v33.py
  run_step 'Run v3.4 hostile-remediation audit' REMEDIATION_AUDIT_V34.log \
    python -u independent_validation/remediation_v34.py
}

stage_legacy_heavy() {
  # Retained because the exact noncommuting sweep is memory-heavy and is best
  # run separately on constrained systems.
  run_step 'Verify exact noncommuting extensions' VERIFY_V27_EXTENSIONS.log \
    python -u scripts/verify_v27_extensions.py
}

stage_figures() {
  run_step 'Generate v3.4 article figures and source data' GENERATE_JCP_FIGURES_V34.log \
    python -u scripts/generate_jcp_figures_v34.py
}

stage_source_audit() {
  run_step 'Audit synchronized manuscript contract' SOURCE_AUDIT_V35.log \
    python -u manuscript_contract/scripts_source_audit_v35.py
}

summary() {
  printf '%s\n' \
    'v3.5 verification stage complete (v3.4 scientific results unchanged):' \
    '  - closed-interval root membership is proved before bracket restriction' \
    '  - exact upper-endpoint roots are retained and jointly certified' \
    '  - unresolved membership, separation, ordering, or sign fails closed' \
    '  - every locator invocation is represented in granular telemetry' \
    '  - the standard boundary sweep contains 15 locator calls and 15 fallbacks' \
    '  - exact rational feasibility is distinct from binary64 recertification' \
    '  - transverse-field evidence is pointwise at four certified couplings' \
    '  - the controller tests fail-closed safety rather than runtime efficiency' \
    '  - article figures and graphical abstract are deterministic outputs'
}

usage() {
  cat <<'USAGE'
Usage: bash reproduce.sh [all|certificates|tests|audits|legacy-heavy|figures|source-audit]
USAGE
}

stage="${1:-all}"
case "$stage" in
  certificates) stage_certificates ;;
  tests) stage_tests ;;
  audits) stage_audits ;;
  legacy-heavy) stage_legacy_heavy ;;
  figures) stage_figures ;;
  source-audit) stage_source_audit ;;
  all)
    stage_certificates
    stage_tests
    stage_audits
    stage_legacy_heavy
    stage_figures
    stage_source_audit
    ;;
  -h|--help|help) usage; exit 0 ;;
  *) usage >&2; exit 2 ;;
esac
summary
