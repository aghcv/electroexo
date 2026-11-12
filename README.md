# ElectroExo

ElectroExo is a Python package for simulating extracellular vesicle release
driven by nanosecond pulsed electric fields. A YAML scenario defines the pulse
train, exposure, cell state, time grid, and extracellular medium. The simulator
exports intracellular state, vesicle release and concentration, and quality
and manufacturing proxies. The mechanistic model and parameters are
exploratory research defaults.

## Run a demo

Clone the `stable` branch for a cluster or shared development checkout.

```bash
git clone --branch stable --single-branch https://github.com/aghcv/electroexo.git
cd electroexo
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -c requirements/cluster-py312.txt .
electro-exocytosis examples/scenario_demo.yml --no-plots --out results/demo
```

The demo covers 30 seconds of simulation and produces seven output time
points. Remove `--no-plots` for plots. The package also supports Python 3.10
or newer through its unconstrained project dependencies. The cluster
constraints record the direct dependencies used for the Python 3.12 baseline.

| Example | Purpose |
| --- | --- |
| [Quick demo](examples/scenario_demo.yml) | Short, headless verification run |
| [Baseline scenario](examples/scenario_baseline.yaml) | Two-hour cell-based release simulation |
| [Direct EV scenario](examples/scenario_direct_EV_engineering.yaml) | Exercise the direct-engineering mode |
| [Extracellular decline](examples/scenario_extracellular_decline.yml) | Pair with `--parameters-file examples/parameters_extracellular_decline.yml` |
| [Observation bridge](examples/experimental_bridge_demo.yml) | Illustrative concentration-to-cell conversion |
| [Parameter sweep](examples/solution_space_analysis/README.md) | Generate scenarios and compare model outputs |
| [Slurm job](examples/slurm/run_demo.sbatch) | Run the quick demo on a compute node |

The `examples/` directory also contains scripts for the dosimetry, membrane,
ion-transport, repair, and release modules, along with plotting and spatial
examples. Optional 3D visualization requires the `viz` extra.

## Run on Slurm

From the checkout root, submit the short example using the installed Python
environment. Set `ELECTROEXO_OUTPUT_ROOT` to a scratch directory when the
cluster provides one.

```bash
export ELECTROEXO_PYTHON="$PWD/.venv/bin/python"
sbatch examples/slurm/run_demo.sbatch
```

The job prints a unique output directory and copies the exact scenario there.
Pass a scenario path as the first job argument to use another input. See the
[developer handoff](docs/developer_handoff.md#cluster-operation) for deployment
and troubleshooting.

## Outputs

| Artifact | Contents |
| --- | --- |
| `summary.json` | Headline simulation metrics |
| `state_timeseries.csv` | Intracellular state and pathway diagnostics |
| `ev_outputs.csv` | Release rates, cumulative outputs, extracellular stocks, and assay concentrations |
| `parameters_used.yaml` | Resolved parameter snapshot |
| `run_metadata.json` | Scenario name, time-point count, and warnings |
| `abbreviations.*` | Glossary for output identifiers |
| `plots/` | Optional plots |

## Develop

```bash
python -m pip install -c requirements/cluster-py312.txt -e '.[dev]'
python -m ruff check .
python -m pytest
```

Use the [input schema](docs/input_schema.md) for units and file structure.
The [developer handoff](docs/developer_handoff.md) maps the implementation and
describes how to build a natural-language interface around validated scenario
edits. The [model assumptions](docs/model_assumptions.md) and
[extracellular kinetics](docs/extracellular_ev_kinetics.md) describe the
scientific model.
