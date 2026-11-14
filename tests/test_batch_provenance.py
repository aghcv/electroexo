"""Exercise the batch process boundary without a scheduler or container runtime."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BATCH = ROOT / "examples" / "slurm" / "run_demo.sbatch"


@pytest.fixture
def batch_environment(tmp_path: Path) -> dict[str, str]:
    package = tmp_path / "installed package" / "electro_exocytosis"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("", encoding="utf-8")
    (package / "__main__.py").write_text(
        "import os, pathlib, sys\n"
        "pathlib.Path('cli_started').write_text(str(pathlib.Path.cwd()))\n"
        "sys.exit(int(os.environ.get('FAKE_SIMULATOR_EXIT', '0')))\n",
        encoding="utf-8",
    )
    return {
        **os.environ,
        "ELECTROEXO_PYTHON": sys.executable,
        "ELECTROEXO_OUTPUT_ROOT": str(tmp_path / "batch results"),
        "PYTHONPATH": str(package.parent),
        "SLURM_SUBMIT_DIR": str(ROOT),
        "SLURM_JOB_ID": "test-123",
        "SLURM_CPUS_PER_TASK": "1",
    }


def test_batch_isolates_repeated_runs_and_records_actual_execution(batch_environment):
    for _ in range(2):
        result = subprocess.run(
            ["bash", str(BATCH)], cwd=ROOT, env=batch_environment,
            capture_output=True, text=True, check=False,
        )
        assert result.returncode == 0, result.stderr

    runs = list(Path(batch_environment["ELECTROEXO_OUTPUT_ROOT"]).glob("run-*"))
    assert len(runs) == 2
    for run in runs:
        provenance = json.loads((run / "run_provenance.json").read_text())
        assert provenance["exit_code"] == 0
        assert provenance["scheduler"]["SLURM_JOB_ID"] == "test-123"
        assert provenance["working_directory"] == str(run)
        assert (run / "cli_started").read_text() == str(run)
        assert provenance["package_path"].startswith(batch_environment["PYTHONPATH"])
        assert provenance["scenario_sha256"] == hashlib.sha256(
            (run / "scenario.yml").read_bytes()
        ).hexdigest()
        assert provenance["command"][-2:] == ["--out", str(run)]
        assert (run / "launcher_sources" / "run_demo.sbatch").is_file()
        assert provenance["git_commit"]
        assert provenance["thread_limits"]["OPENBLAS_NUM_THREADS"] == "1"


def test_batch_preserves_simulator_failure(batch_environment):
    batch_environment["FAKE_SIMULATOR_EXIT"] = "23"
    result = subprocess.run(
        ["bash", str(BATCH)], cwd=ROOT, env=batch_environment,
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 23, result.stderr
    run, = Path(batch_environment["ELECTROEXO_OUTPUT_ROOT"]).glob("run-*")
    provenance = json.loads((run / "run_provenance.json").read_text())
    assert provenance["exit_code"] == 23
    assert provenance["finished_at"]
