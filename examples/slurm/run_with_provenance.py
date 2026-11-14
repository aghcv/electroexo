"""Run the installed simulator outside the checkout and retain deployment provenance."""

from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def git_output(repo: Path, *arguments: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo), *arguments],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return result.stdout.rstrip("\n") if result.returncode == 0 else None


def run(scenario: Path, outdir: Path, repo: Path) -> int:
    scenario, outdir, repo = scenario.resolve(), outdir.resolve(), repo.resolve()
    command = [
        sys.executable, "-m", "electro_exocytosis", "run", str(scenario),
        "--no-plots", "--out", str(outdir),
    ]
    image = os.environ.get("APPTAINER_CONTAINER") or os.environ.get("SINGULARITY_CONTAINER")
    status = git_output(repo, "status", "--porcelain", "--untracked-files=all")
    patch = git_output(repo, "diff", "--binary", "HEAD")
    if patch is not None:
        (outdir / "source_changes.patch").write_text(patch + "\n", encoding="utf-8")
    launcher_dir = outdir / "launcher_sources"
    launcher_dir.mkdir()
    for name in (
        "run_demo.sbatch", "run_wahab_demo.sbatch", "wahab_python.sh",
        "run_with_provenance.py", "check_demo.py",
    ):
        source = repo / "examples" / "slurm" / name
        if source.is_file():
            shutil.copy2(source, launcher_dir / name)
    # Resolve the package in the same working directory used by the child CLI.
    os.chdir(outdir)
    spec = importlib.util.find_spec("electro_exocytosis")
    provenance = {
        "started_at": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "working_directory": str(outdir),
        "scenario_sha256": hashlib.sha256(scenario.read_bytes()).hexdigest(),
        "repository": str(repo),
        "git_commit": git_output(repo, "rev-parse", "HEAD"),
        "git_status": status,
        "git_dirty": bool(status) if status is not None else None,
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "package_path": spec.origin if spec else None,
        "packages": dict(sorted(
            (dist.metadata["Name"], dist.version)
            for dist in importlib.metadata.distributions()
            if dist.metadata["Name"]
        )),
        "loaded_modules": os.environ.get("LOADEDMODULES"),
        "container_image": image,
        "container_image_resolved": str(Path(image).resolve()) if image else None,
        "scheduler": {
            name: os.environ.get(name) for name in (
                "SLURM_JOB_ID", "SLURM_JOB_PARTITION", "SLURM_JOB_ACCOUNT",
                "SLURM_CPUS_PER_TASK", "SLURM_SUBMIT_DIR",
            )
        },
        "thread_limits": {
            name: os.environ.get(name)
            for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        },
        "exit_code": None,
    }
    provenance_file = outdir / "run_provenance.json"

    def save() -> None:
        provenance_file.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    save()
    try:
        result = subprocess.run(command, cwd=outdir, check=False)
        exit_code = result.returncode if result.returncode >= 0 else 128 - result.returncode
    except OSError as exc:
        provenance["error"] = str(exc)
        exit_code = 127
    provenance["exit_code"] = exit_code
    provenance["finished_at"] = datetime.now(timezone.utc).isoformat()
    save()
    return exit_code


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: run_with_provenance.py SCENARIO OUTPUT_DIRECTORY CHECKOUT")
    sys.exit(run(*(Path(argument) for argument in sys.argv[1:])))
