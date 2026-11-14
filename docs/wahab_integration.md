# Wahab integration and smoke test

Use Wahab's Python 3.12 container to build a project environment, then run
the installed simulator through Slurm. The first acceptance target is the
30-second `release_demo`, with seven output times from 0 to 30 seconds.
These checks establish that installation, scheduling, and artifact retrieval
work; the model defaults remain exploratory and uncalibrated.

## Create the environment

Start in the checkout root. Record the revision and any local changes before
building an environment for a shared run.

```bash
git rev-parse HEAD
git status --short
module load container_env python3/2026.0-py312
crun python -m venv "$PWD/.venv"
crun -p "$PWD/.venv" python -m pip install \
    -c requirements/cluster-py312.txt '.[dev]' build
crun -p "$PWD/.venv" python -m pip check
```

Create the environment once and install dependencies before submitting jobs.
The constraints pin the direct runtime dependencies; retain the resolved
dependency versions with each run. On nodes without network access, prepare
wheels on a compatible networked node and install with `--no-index
--find-links /path/to/wheels`.

The environment is isolated from the container's preinstalled Python
packages. Its Python symlinks can point into `/opt/conda`, which exists inside
the container. Run it through
[`wahab_python.sh`](../examples/slurm/wahab_python.sh), rather than executing
`.venv/bin/python` on the host:

```bash
examples/slurm/wahab_python.sh --version
```

The launcher loads `container_env` and `python3/2026.0-py312` and forwards
Python arguments through `crun`. Set `ELECTROEXO_ENV` to an absolute path to
use another project environment. Set `ELECTROEXO_WAHAB_MODULE` to select a
different reviewed module, and recreate the environment for that runtime.
Wahab's [Python documentation](https://wiki.hpc.odu.edu/en/Software/Python)
describes the container system.

## Run development checks

Use a compute allocation for the full test suite. Its time and memory budget
is separate from the small demo job:

```bash
srun --partition=main --nodes=1 --ntasks=1 --cpus-per-task=1 \
    --mem=4G --time=00:30:00 bash -c \
    'examples/slurm/wahab_python.sh -m ruff check . &&
     examples/slurm/wahab_python.sh -m pytest'
```

Add your approved `--account` or other site options where required. Adjust
the allocation to measured suite requirements. The launcher tests use fake
modules and a fake container command, so they also run in ordinary Linux CI.

## Submit the demo

Submit from the checkout root:

```bash
job_id=$(sbatch --parsable examples/slurm/run_wahab_demo.sbatch)
squeue --jobs="$job_id"
sacct --jobs="$job_id" --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

The job requests one CPU, 1 GB RAM, five minutes, and the `main` partition.
Slurm command-line options override those defaults. It limits numerical
library threads to the CPU allocation and disables plots. A different
scenario can be supplied as the first script argument; the demo checker
below is specifically for the unchanged `release_demo`.

By default, results go into a unique directory beneath
`/scratch/$USER/electroexo`. Set `ELECTROEXO_OUTPUT_ROOT` before submission to
use another writable results location. The Slurm log `slurm-<job-id>.out`
prints the exact run directory. Wahab's scratch storage is temporary and is
not backed up; copy retained results to your project's durable storage.
See [ODU storage guidance](https://wiki.hpc.odu.edu/en/Cluster/Storage).

## Verify the result

Require the top-level job to reach `COMPLETED` with `ExitCode` `0:0`, then
check the run directory reported in its log:

```bash
examples/slurm/wahab_python.sh examples/slurm/check_demo.py /path/to/run
```

The checker requires the standard JSON, CSV, parameter, and abbreviation
artifacts; the scenario must be `release_demo`, with seven time points at
0, 5, 10, 15, 20, 25, and 30 seconds. Numerical outputs must be finite.
Model warnings are reported for review. Warnings about uncalibrated models
or disabled extracellular loss mechanisms are expected for this demo and do
not by themselves indicate a deployment failure. A `quality_pass` value of
false is also a model result, rather than a scheduler failure.

Compare against a retained demo from the same input and parameter snapshot
using explicit numerical tolerances:

```bash
examples/slurm/wahab_python.sh examples/slurm/check_demo.py /path/to/run \
    --reference /path/to/reference --rtol 1e-7 --atol 1e-10
```

Choose and record tolerances suitable for the software and hardware under
comparison. Timestamps are not numerical acceptance criteria.

Each batch run retains the copied scenario, resolved parameters, and
`run_provenance.json` alongside the outputs. Provenance records the command,
scenario hash, Git revision and working-tree state, Python and dependency
versions, imported package location, loaded modules and container identity,
scheduler job ID, and process exit status. Retain the Slurm log as well. A
dirty checkout identifies an uncommitted build; review its recorded changes
before using it as a shared reference. Execution from the run directory also
checks that the installed package can load its packaged defaults without
importing the checkout accidentally.

Before promoting the setup, submit the demo twice and check that both jobs
produce distinct run directories. Test a copied scenario with
`simulation.output_dt_s: 0`: it must fail with a nonzero exit and must not
produce a successful artifact set. A missing project environment must fail
with a clear launcher error. Keep failed-run logs and provenance for
diagnosis.

After these checks pass, run a representative longer scenario and measure
runtime, peak memory, and output size before selecting sweep resources.
Implement service job states, retries, and cancellation around this tested
execution path before connecting a language-model interface, following the
[developer handoff](developer_handoff.md#natural-language-interface).
