# Developer handoff

## Runtime and ownership

| Area | Entry point | Responsibility |
| --- | --- | --- |
| Configuration | `config.SimulationScenario`, `io.readers.load_scenario` | Parse scenario YAML, validate typed fields and time settings, supply defaults |
| Parameters | `io.readers.load_default_parameters`, `merge_parameters`, `parameter_registry` | Resolve sparse overrides, expose units and runtime parameter status |
| Simulation | `simulation.Simulation(...).run()` | Resolve pulse descriptors, integrate the state equations, return a `SimulationResult` |
| Observation | `experimental_bridge`, `models.extracellular_kinetics`, `models.ev_size_observation` | Convert between scales, propagate extracellular stocks, map release to measured size bands |
| Sensitivity | `experiment_design` | Model-agnostic Morris screening, Sobol estimation, local curvature, and ranking utilities |
| Outputs | `io.writers`, `visualization` | Export stable artifacts and generate optional plots |
| Entry points | `electro-exocytosis`, `python -m electro_exocytosis run` | Run an input file and write outputs to a specified directory |

The module paths above are relative to `electro_exocytosis`. Keep the
scientific equations, parameter defaults, and output conventions in this
package. A service or interface should call the same entry points as a
command-line user.

```python
from electro_exocytosis.io.readers import load_scenario
from electro_exocytosis.simulation import Simulation

scenario = load_scenario("examples/scenario_demo.yml")
result = Simulation(scenario).run()
print(result.summary)
```

See the [model assumptions](model_assumptions.md), [size observation
model](size_resolved_ev_observation.md), and [parameter registry
guide](parameter_defaults_and_snapshots.md) when changing a model. The
calibration and example settings remain exploratory. Record units, source,
bounds, and the expected effect on outputs with each scientific change.

## Cluster operation

1. Clone `stable` from the repository. Record the checked-out commit with
   `git rev-parse HEAD`. Use that commit when recreating the environment.
2. Create a Python 3.12 virtual environment and install with
   `python -m pip install -c requirements/cluster-py312.txt .`.
   For compute nodes without network access, build or download wheels on a
   compatible networked node, then install from the wheel directory using
   `--no-index --find-links`.
3. Run the headless demo and inspect `summary.json`, `run_metadata.json`, and
   the time-series files. The demo should finish with seven time points.
4. From the checkout root, set `ELECTROEXO_PYTHON` to the virtual environment's
   Python executable. Set `ELECTROEXO_OUTPUT_ROOT` to a writable results or
   scratch directory. Submit `sbatch examples/slurm/run_demo.sbatch`. Add the
   site's account or partition options as required.
5. Check job state with `squeue` and inspect the Slurm log and run directory.
   If the job fails, check the Python executable, input validation errors,
   output permissions, and scheduler memory and time limits. The example
   requests one CPU and limits numerical library threads to the allocation.

For published or shared runs, keep the immutable input, resolved parameter
snapshot, commit ID, Python and dependency versions, command, scheduler job
ID, and output together. Run numerical comparisons with explicit tolerances
and record the software and platform used. Build packages on a node
compatible with the target environment.

## Development workflow

Work on a feature branch and open a pull request against `stable`. Run
`python -m ruff check .` and `python -m pytest` before review. The workflow
builds the wheel, runs tests, and exercises the installed CLI on Linux.
Python tests use synthetic inputs and temporary directories. Add a focused
test for each behavioral or scientific change and keep the small demo usable
as an integration test. Keep generated outputs, local environment files,
study data, and credentials in their designated external storage.

Treat scenario fields, file outputs, and the Python entry points as
interfaces used by other developers. Document any deliberate compatibility
change in the pull request and add a migration example. When editing a
scientific model, compare a representative baseline and edge cases against
the previous revision before promoting it on the cluster.

## Natural-language interface

Build the interface as a separate service around the deterministic simulator.
The language model proposes a typed edit. Application code validates the
edit, presents the resolved change, and owns job submission and output
retrieval.

### Scenario edit contract

Start with a reviewed demo template and a small allowlist of editable fields.
A first version can support changing the output interval and duration. For
example, a structured proposal can carry a numeric value in the canonical unit.

```json
{
  "template_id": "release_demo",
  "changes": [
    {"field": "simulation.output_dt_s", "value": 2.0}
  ]
}
```

- Define a strict schema for proposals. Validate allowed field paths and value
  types, reject unknown keys and nonfinite values, and convert explicit user
  units to schema units. Ask for clarification when the unit or requested
  change is ambiguous.
- Apply the change to a copy of a server-owned template and revalidate the
  entire object with `SimulationScenario.model_validate`. The existing
  scenario loader has permissive handling of unknown fields, so enforce the
  edit allowlist in the interface.
- Enforce deployment-specific duration, output-point, memory, and queue
  budgets before allocating simulation work. Set these limits with the
  cluster administrators. Give the user the resolved before-and-after values,
  units, and estimated output size for review.
- Serialize validated settings with `yaml.safe_dump` into a new per-run
  directory. Keep parameter overrides behind a separate reviewed interface
  backed by the parameter registry.

### Execution and results

- Use a fixed worker entry point with a known Python environment and a fixed
  argument vector. The worker owns paths, job names, scheduler options, and
  environment settings. Scope its permissions to the template directory and
  its run directory.
- Start with a local worker for a small demo, then add a Slurm adapter.
  Represent runs with durable states such as `validated`, `queued`,
  `running`, `succeeded`, and `failed`. Use an idempotency key for submission
  retries and retain the scheduler job ID for cancellation and status checks.
- Keep the validated YAML, merged parameters, software revision, and artifact
  paths with each run. Expose a small allowlist of outputs from the run
  directory. Return warnings and observed failure messages with concise
  recovery steps. Restrict logs to the structured run information needed
  for debugging and provenance.
- Use headless execution for jobs. Render plots from completed results as a
  separate operation when the interface needs them.

### Implementation milestones

1. Build a deterministic edit-and-validate function. Test field selection,
   units, invalid values, unchanged templates, and resource limits.
2. Execute one validated demo locally. Test artifact paths, failed validation,
   repeated submissions, timeouts, and per-run directory isolation.
3. Add a Slurm adapter and test its state transitions with a fake scheduler.
   Verify an end-to-end demo on a development allocation.
4. Connect the language model to the structured edit endpoint. Use synthetic
   requests to verify ambiguous-unit handling, diff review, duplicate-job
   protection, failure reporting, and retrieval of permitted artifacts.
