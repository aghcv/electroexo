from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "examples" / "slurm" / "wahab_python.sh"


@pytest.fixture
def launcher_environment(tmp_path: Path) -> dict[str, str]:
    environment = tmp_path / "project environment"
    (environment / "bin").mkdir(parents=True)
    (environment / "pyvenv.cfg").write_text("version = 3.12\n", encoding="utf-8")
    (environment / "bin" / "activate").touch()
    executables = tmp_path / "bin"
    executables.mkdir()
    module = executables / "module"
    module.write_text(
        '#!/bin/sh\nprintf "%s\\n" "$@" > "$TEST_MODULE_LOG"\n'
        'exit "${TEST_MODULE_EXIT_CODE:-0}"\n',
        encoding="utf-8",
    )
    module.chmod(0o755)
    crun = executables / "crun"
    crun.write_text(
        f"#!{sys.executable}\n"
        "import json, os, sys\n"
        'with open(os.environ["TEST_CRUN_LOG"], "w") as handle:\n'
        "    json.dump(sys.argv[1:], handle)\n"
        'sys.exit(int(os.environ.get("TEST_CRUN_EXIT_CODE", "0")))\n',
        encoding="utf-8",
    )
    crun.chmod(0o755)
    env = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith("BASH_FUNC_") and key not in ("BASH_ENV", "ENV", "SHELLOPTS")
    }
    env.update(
        PATH=f"{executables}{os.pathsep}{os.environ['PATH']}",
        ELECTROEXO_ENV=str(environment),
        TEST_MODULE_LOG=str(tmp_path / "module.log"),
        TEST_CRUN_LOG=str(tmp_path / "crun.json"),
    )
    env.pop("ELECTROEXO_WAHAB_MODULE", None)
    return env


def test_wahab_launcher_forwards_arguments(launcher_environment: dict[str, str]) -> None:
    env = launcher_environment
    arguments = ["-m", "electro_exocytosis", "run", "scenario with spaces.yml", "--no-plots"]

    result = subprocess.run([str(LAUNCHER), *arguments], env=env, capture_output=True, text=True)

    assert result.returncode == 0, result.stderr
    assert Path(env["TEST_MODULE_LOG"]).read_text().splitlines() == [
        "load", "container_env", "python3/2026.0-py312"
    ]
    assert json.loads(Path(env["TEST_CRUN_LOG"]).read_text()) == [
        "-p", env["ELECTROEXO_ENV"], "python", *arguments
    ]


def test_wahab_launcher_rejects_missing_environment(
    launcher_environment: dict[str, str], tmp_path: Path
) -> None:
    env = launcher_environment
    env["ELECTROEXO_ENV"] = str(tmp_path / "missing")

    result = subprocess.run([str(LAUNCHER), "--version"], env=env, capture_output=True, text=True)

    assert result.returncode == 2
    assert "environment is missing or incomplete" in result.stderr
    assert not Path(env["TEST_MODULE_LOG"]).exists()
    assert not Path(env["TEST_CRUN_LOG"]).exists()


def test_wahab_launcher_resolves_default_environment_outside_checkout(
    launcher_environment: dict[str, str], tmp_path: Path
) -> None:
    env = launcher_environment
    checkout = tmp_path / "checkout"
    launcher = checkout / "examples" / "slurm" / LAUNCHER.name
    launcher.parent.mkdir(parents=True)
    shutil.copy2(LAUNCHER, launcher)
    environment = checkout / ".venv"
    shutil.copytree(env.pop("ELECTROEXO_ENV"), environment)

    result = subprocess.run(
        [str(launcher), "--version"], cwd=tmp_path, env=env, capture_output=True, text=True
    )

    assert result.returncode == 0, result.stderr
    assert json.loads(Path(env["TEST_CRUN_LOG"]).read_text()) == [
        "-p", str(environment), "python", "--version"
    ]


def test_wahab_batch_forwards_scenario_and_output_override(tmp_path: Path) -> None:
    checkout = tmp_path / "checkout"
    slurm = checkout / "examples" / "slurm"
    slurm.mkdir(parents=True)
    batch = slurm / "run_wahab_demo.sbatch"
    shutil.copy2(REPO / "examples" / "slurm" / batch.name, batch)
    (slurm / "run_demo.sbatch").write_text(
        '#!/bin/sh\nprintf "%s\\n" "$ELECTROEXO_PYTHON" "$ELECTROEXO_OUTPUT_ROOT" "$@"\n',
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["SLURM_SUBMIT_DIR"] = str(checkout)
    env["ELECTROEXO_OUTPUT_ROOT"] = str(tmp_path / "custom results")

    result = subprocess.run(
        [str(batch), "scenario with spaces.yml"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [
        str(slurm / "wahab_python.sh"),
        env["ELECTROEXO_OUTPUT_ROOT"],
        "scenario with spaces.yml",
    ]


def test_wahab_launcher_propagates_container_failure(
    launcher_environment: dict[str, str],
) -> None:
    env = launcher_environment
    env["TEST_CRUN_EXIT_CODE"] = "17"
    env["ELECTROEXO_WAHAB_MODULE"] = "python3/2025.1-py312"

    result = subprocess.run([str(LAUNCHER), "--version"], env=env, capture_output=True, text=True)

    assert result.returncode == 17
    assert Path(env["TEST_MODULE_LOG"]).read_text().splitlines()[-1] == "python3/2025.1-py312"


def test_wahab_launcher_stops_after_module_failure(
    launcher_environment: dict[str, str],
) -> None:
    env = launcher_environment
    env["TEST_MODULE_EXIT_CODE"] = "9"

    result = subprocess.run([str(LAUNCHER), "--version"], env=env, capture_output=True, text=True)

    assert result.returncode == 9
    assert not Path(env["TEST_CRUN_LOG"]).exists()
