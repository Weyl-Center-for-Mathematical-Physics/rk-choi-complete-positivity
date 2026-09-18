"""Contract tests for reproduce_bit.py (Task 13): CLI handling, fail-fast behaviour, output-path confinement,
manifest completeness."""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import reproduce_bit as rb  # noqa: E402

PY = sys.executable


def test_unknown_stage_is_rejected():
    with pytest.raises(SystemExit) as exc:
        rb.main(["bogus-stage"])
    assert exc.value.code == 2


def test_list_prints_exact_commands(capsys):
    assert rb.main(["all", "--list"]) == 0
    out = capsys.readouterr().out
    for stage in rb.STAGES:
        assert f"[{stage}]" in out
    assert "scripts/generate_bit_figures.py" in out and "scripts/verify_bit_claims.py" in out
    for line in out.splitlines():
        if line.startswith("  "):
            assert line.strip().startswith(PY), line


def test_log_dir_must_be_inside_archive(tmp_path):
    assert rb.main(["figures", "--log-dir", str(tmp_path)]) == 2


def test_manuscript_contract_needs_a_manuscript_directory(tmp_path):
    assert rb.main(["manuscript-contract", "--manuscript", str(tmp_path)]) == 2


def test_fail_fast_stops_at_first_failure_and_logs_it(monkeypatch):
    log_dir = rb.RESULTS / "logs" / "_pytest_tmp"
    sentinel = ROOT / "results" / "bit_revision" / "_should_not_exist.txt"
    if sentinel.exists():
        sentinel.unlink()

    def fake(stage, manuscript, submission, audit_json):
        return [
            ("ok", [PY, "-c", "print('first ok')"]),
            ("boom", [PY, "-c", "import sys; print('failing'); sys.exit(7)"]),
            ("never", [PY, "-c", f"open(r'{sentinel}', 'w').write('x')"]),
        ]

    monkeypatch.setattr(rb, "stage_commands", fake)
    try:
        rc = rb.main(["certificates", "--log-dir", str(log_dir)])
        assert rc == 7
        logs = sorted(log_dir.glob("*_certificates.log"))
        assert logs, "no log written"
        text = logs[-1].read_text(encoding="utf-8")
        assert "first ok" in text and "exit 7" in text and "==> never" not in text
        assert not sentinel.exists()
    finally:
        shutil.rmtree(log_dir, ignore_errors=True)
        if sentinel.exists():
            sentinel.unlink()


def test_manifest_lists_every_artifact_with_correct_hashes(tmp_path):
    records = [{"stage": "figures", "status": "PASS", "commands": [{"title": "t", "argv": [], "exit": 0, "seconds": 0.0, "tests_passed": 3}]}]
    target = tmp_path / "manifest.json"
    path = rb.write_manifest(records, target)
    m = json.loads(path.read_text(encoding="utf-8"))
    assert m["stage_counts"] == {"PASS": 1, "FAIL": 0, "SKIPPED": 0}
    assert m["tests_passed"] == 3
    expected = set()
    for pattern in rb.ARTIFACT_GLOBS:
        for p in ROOT.glob(pattern):
            if p.is_file():
                expected.add(p.relative_to(ROOT).as_posix())
    assert set(m["artifacts"]) == expected and expected
    for rel, entry in m["artifacts"].items():
        assert entry["sha256"] == hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
    for rel in rb.FROZEN_INPUTS:
        assert rel in m["frozen_inputs_sha256"], rel
        assert m["frozen_inputs_sha256"][rel] == hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def test_default_log_dir_and_manifest_are_confined_to_the_archive():
    assert rb._inside(rb.DEFAULT_LOG_DIR, rb.ROOT)
    assert rb._inside(rb.MANIFEST, rb.ROOT)
    assert not rb._inside(rb.ROOT.parent / "manuscript", rb.ROOT)


def test_manifest_binds_current_computation_sources(tmp_path):
    path = rb.write_manifest([], tmp_path / "manifest.json")
    m = json.loads(path.read_text())
    hashes = m.get("source_inputs_sha256", {})
    for rel in ("reproduce_bit.py", "pyproject.toml", "scripts/verify_bit_claims.py", "src/rk_choi_margin/__init__.py",
                "manuscript_contract/scripts_source_audit_v35.py"):
        assert hashes.get(rel) == hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def test_all_fails_when_collection_exceeds_executed_tests(tmp_path, monkeypatch):
    monkeypatch.setattr(rb, "STAGES", ["tests"])
    monkeypatch.setattr(rb, "stage_commands", lambda *args: [
        ("collection", [PY, "-c", "print('2 tests collected')", "--collect-only", "pytest"]),
        ("run", [PY, "-c", "print('1 passed')", "pytest"]),
    ])
    monkeypatch.setattr(rb, "write_manifest", lambda *args, **kwargs: tmp_path / "manifest.json")
    assert rb.main(["all"]) == 1


def test_standalone_all_accounts_for_explicit_skipped_workspace_tests(tmp_path, monkeypatch):
    monkeypatch.setattr(rb, "ROOT", tmp_path)
    monkeypatch.setattr(rb, "DEFAULT_LOG_DIR", tmp_path / "logs")
    monkeypatch.setattr(rb, "STAGES", ["tests"])
    monkeypatch.setattr(rb, "stage_commands", lambda *args: [
        ("collection", [PY, "-c", "print('2 tests collected')", "--collect-only", "pytest"]),
        ("run", [PY, "-c", "print('1 passed, 1 skipped')", "pytest"]),
    ])
    monkeypatch.setattr(rb, "write_manifest", lambda *args, **kwargs: tmp_path / "manifest.json")
    assert rb.main(["all"]) == 0


@pytest.mark.parametrize("kind", ["frozen", "manuscript", "source"])
def test_all_rejects_midrun_input_changes_and_retains_start_hash(tmp_path, monkeypatch, kind):
    code = tmp_path / "code"
    code.mkdir()
    manuscript = tmp_path / "manuscript"
    manuscript.mkdir()
    (manuscript / "main.tex").write_text("original")
    (code / "frozen.txt").write_text("original")
    (code / "reproduce_bit.py").write_text("original")
    targets = {"frozen": code / "frozen.txt", "manuscript": manuscript / "main.tex", "source": code / "reproduce_bit.py"}
    target = targets[kind]
    initial = hashlib.sha256(target.read_bytes()).hexdigest()
    monkeypatch.setattr(rb, "ROOT", code)
    monkeypatch.setattr(rb, "DEFAULT_LOG_DIR", code / "logs")
    monkeypatch.setattr(rb, "FROZEN_INPUTS", ["frozen.txt"])
    monkeypatch.setattr(rb, "STAGES", ["certificates"])
    monkeypatch.setattr(rb, "stage_commands", lambda *args: [
        ("mutation", [PY, "-c", f"from pathlib import Path; Path({str(target)!r}).write_text('changed')"]),
    ])
    real_write = rb.write_manifest
    monkeypatch.setattr(rb, "write_manifest", lambda records, **kwargs: real_write(records, code / "manifest.json", **kwargs))
    assert rb.main(["all"]) == 1
    m = json.loads((code / "manifest.json").read_text())
    field = {"frozen": "frozen_inputs_sha256", "manuscript": "manuscript_inputs_sha256", "source": "source_inputs_sha256"}[kind]
    assert m[field][target.name] == initial
