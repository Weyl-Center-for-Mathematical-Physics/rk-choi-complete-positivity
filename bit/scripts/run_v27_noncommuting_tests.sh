#!/usr/bin/env bash
set -euo pipefail
export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
nodes=(
  'tests/test_v27_noncommuting_adaptive.py::test_noncommuting_transverse_field_preserves_disconnected_components[0]'
  'tests/test_v27_noncommuting_adaptive.py::test_noncommuting_transverse_field_preserves_disconnected_components[0.05]'
  'tests/test_v27_noncommuting_adaptive.py::test_noncommuting_transverse_field_preserves_disconnected_components[0.1]'
  'tests/test_v27_noncommuting_adaptive.py::test_noncommuting_transverse_field_preserves_disconnected_components[0.2]'
  'tests/test_v27_noncommuting_adaptive.py::test_noncommuting_strict_sign_samples_are_exact'
  'tests/test_v27_noncommuting_adaptive.py::test_adaptive_error_only_can_accept_non_cptp_steps'
  'tests/test_v27_noncommuting_adaptive.py::test_certified_guard_or_rotating_frame_avoid_false_physical_steps'
  'tests/test_v27_noncommuting_adaptive.py::test_buffered_guard_recomputes_error_at_actual_candidate'
)
passed=0
for node in "${nodes[@]}"; do
  python -m pytest -q -p no:cacheprovider "$node"
  passed=$((passed+1))
done
printf '%s\n' "$passed v2.7 noncommuting/adaptive tests passed in isolated processes"
