# electro_exocytosis

`electro_exocytosis` is a research-oriented Python package for simulating nanosecond pulsed electric field (nsPEF) driven extracellular vesicle (EV) release. It turns a structured scenario definition into a multi-layer dynamical simulation spanning pulse delivery, electrodynamics, calcium signaling, stress remodeling, EV subtype release, cargo/potency proxies, injury gating, and manufacturing/QC outputs.

## Scientific motivation

nsPEF exposures can perturb plasma and organelle membranes, trigger calcium mobilization, alter ROS and bioenergetics, and reshape vesicle release pathways. These mechanisms matter both for understanding electro-exocytosis as a biological phenomenon and for designing EV engineering workflows that bias yield, subtype balance, and product quality.

This package provides a first professional implementation of that conceptual framework. Version 0.1.0 is intentionally conservative: it offers a coherent, testable, end-to-end simulator while clearly marking all mechanistic modules as placeholders pending literature-backed parameterization.

## Modular solver architecture

The model is organized around eight biological/computational layers plus a cross-cutting cell-state modifier module:

1. **Pulse delivery, exposure geometry, dosimetry**
2. **Plasma membrane and organelle electrodynamics**
3. **Ion transport, Ca2+ mobilization, ROS, bioenergetics**
4. **Ca2+-dependent remodeling and repair**
5. **EV biogenesis and subtype release**
6. **Cargo sorting, EV composition, potency**
7. **Injury, debris, quality gate**
8. **Manufacturing, isolation, QC interface**
9. **Cross-cutting cell-state and disease modifiers**

The simulator resolves nanosecond pulses through descriptors, then integrates a stable ODE system over seconds to hours.

## Installation

```bash
pip install -e ".[dev]"
```

## Quick start

Run the baseline example:

```bash
python -m electro_exocytosis run examples/scenario_baseline.yaml --out results/baseline
```

Other examples:

```bash
python -m electro_exocytosis run examples/scenario_high_dose.yaml --out results/high_dose
python -m electro_exocytosis run examples/scenario_direct_EV_engineering.yaml --out results/direct_ev
```

Compare the pulse/dosimetry model choices directly:

```bash
python examples/compare_dosimetry_models.py --out results/dosimetry_model_comparison
```

This writes CSV summaries and temperature-profile plots for `legacy`,
`joule_adiabatic`, and `joule_lumped_thermal` across several nsPEF scenarios.

Compare reduced membrane/organelle electrodynamics responses directly:

```bash
python examples/compare_membrane_electrodynamics.py --out results/electrodynamics_model_comparison
```

This writes a CSV summary plus pulse-width and cell-radius sensitivity plots
for the reduced Schwan charging, permeability, and pore-density proxies.

Compare ion-transport, organelle, ROS, and ATP responses directly:

```bash
python examples/compare_ion_transport_bioenergetics.py --out results/ion_transport_bioenergetics_comparison
```

This writes summary and time-series CSV files plus manuscript-style plots for
pore-mediated Ca2+ entry, ER release, Na/K/Cl osmotic perturbation,
mitochondrial depolarization, ROS, ATP, and damage coupling.

Compare Ca2+-dependent remodeling and membrane repair responses directly:

```bash
python examples/compare_remodeling_repair.py --out results/remodeling_repair_comparison
```

This writes summary and time-series CSV files plus plots for local Ca2+
microdomain strength, PS externalization, calpain, annexin, lysosomal repair,
actin remodeling, resealing, and repair-associated shedding.

Compare EV-biogenesis, docking, fusion, budding, and apoptotic-release regimes directly:

```bash
python examples/compare_ev_biogenesis_release.py --out results/ev_biogenesis_comparison
```

This writes summary and time-series CSV files plus plots for MVB pool size,
ILV load, docked MVBs, secretory versus lysosomal routing bias, and subtype
resolved EV release.

Generate a three-scenario, five-layer storytelling overview:

```bash
python examples/storyboard_multilayer_overview.py --out results/multilayer_storyboard
```

This writes a paper-style multi-panel figure with three scenarios as rows and
Layers 1-5 as columns, showing how dose, charging, Ca2+/ROS/ATP, repair, and
EV subtype outputs compose into coherent scenario-level behavior.

Map EV-release regimes across amplitude and pulse number:

```bash
python examples/map_ev_release_regimes.py --out results/ev_release_regime_map
```

This writes a parameter-sweep heatmap for small-EV, medium/large-EV, and
apoptotic-body output together with viability across a pulse-dose grid.

## Output files

Each run writes a directory containing:

- `summary.json` – compact run summary and headline metrics
- `state_timeseries.csv` – time series for Ca_i, Ca_ER, Ca_mito, mitochondrial potential, ROS, ATP, Na_i, K_i, Cl_i, osmotic stress, ion-flux diagnostics, remodeling/repair diagnostics, and damage
- `ev_outputs.csv` – EV rates, cumulative counts, viability, quality gate
- `parameters_used.yaml` – merged parameter set used for the run
- `run_metadata.json` – metadata, warnings, timestamps
- `*.png` – plots unless `--no-plots` is used

Framework-generated figures use manuscript-oriented defaults: standard 16:9
landscape layout and 1200 dpi PNG export, unless a specific analysis overrides
the figure standard.
For overlaid line plots, the default style is color-blind friendly and
print-friendly: up to three series use monochrome line styles and markers; plots
with more than three series switch to a color-blind-safe palette while retaining
distinct line formats.

## Placeholder status in v0.1.0

All scientific modules are placeholders in this release. The code structure is real, the numerical workflow is functional, and the interfaces are designed for extension, but the equations and constants should be treated as exploratory defaults rather than validated mechanistic truth.

Placeholder-heavy modules include:

- pulse and dosimetry scaling
- membrane/organelle electrodynamics
- ion transport, calcium mobilization, mitochondrial stress, ROS, and ATP coupling
- remodeling and repair logic
- EV subtype release kinetics remain reduced but now include explicit MVB
  maturation, ILV loading, docking, fusion, budding, and apoptotic-commitment
  states
- cargo/potency proxies
- injury and purity gates
- manufacturing/QC transforms
- cell-state modifiers

## Evidence workbook integration

The repository includes an Excel evidence workbook and PDF source material. The `EvidenceLoader` reads the workbook sheets into pandas DataFrames and can summarize module coverage and placeholder status. Future versions will use curated literature evidence to replace placeholder constants, geometry factors, gating equations, and coupling strengths with module-specific parameter sets.

The first full-text PDF pass is summarized in
`docs/fulltext_calibration_opportunities.md`. The corresponding structured
targets are available through `EvidenceLoader.get_calibration_targets()` and in
`electro_exocytosis/evidence/calibration_targets.csv`.

## Disclaimer

This package is an exploratory research simulation. It is **not** experimentally validated, **not** intended for clinical or regulatory use, and should not be used to make high-stakes biological, medical, or manufacturing decisions without extensive literature review and experimental calibration.
