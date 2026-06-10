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

## Output files

Each run writes a directory containing:

- `summary.json` – compact run summary and headline metrics
- `state_timeseries.csv` – time series for Ca_i, Ca_ER, ROS, ATP, damage
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
- calcium/ROS/ATP coupling
- remodeling and repair logic
- EV subtype release kinetics
- cargo/potency proxies
- injury and purity gates
- manufacturing/QC transforms
- cell-state modifiers

## Evidence workbook integration

The repository includes an Excel evidence workbook and PDF source material. The `EvidenceLoader` reads the workbook sheets into pandas DataFrames and can summarize module coverage and placeholder status. Future versions will use curated literature evidence to replace placeholder constants, geometry factors, gating equations, and coupling strengths with module-specific parameter sets.

## Disclaimer

This package is an exploratory research simulation. It is **not** experimentally validated, **not** intended for clinical or regulatory use, and should not be used to make high-stakes biological, medical, or manufacturing decisions without extensive literature review and experimental calibration.
