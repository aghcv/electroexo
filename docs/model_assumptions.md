# Model assumptions

## Layer 1: Pulse delivery, exposure geometry, dosimetry
Current implementation uses simple field conversion, energy-density scaling, and fixed geometry correction factors. These are placeholders and need literature-backed exposure-specific dosimetry relations.

## Layer 2: Plasma membrane and organelle electrodynamics
A simplified Schwan-style transmembrane voltage estimate is used with fixed organelle fractions and sigmoid permeability. Placeholder support is needed for membrane charging, pore formation, and organelle-specific coupling.

## Layer 3: Ion transport, Ca2+ mobilization, ROS, bioenergetics
The ODE system uses decaying Ca2+ release, first-order uptake, ROS production from Ca2+ overload, and ATP depletion from oxidative stress. All constants and couplings are placeholders.

## Layer 4: Ca2+-dependent remodeling and repair
PS exposure, calpain, annexins, actin disruption, and repair are represented by Hill functions of cytosolic Ca2+. Literature support is needed for kinetics, thresholds, and recovery structure.

## Layer 5: EV biogenesis and subtype release
sEV, mlEV, and apoptotic body release are driven by Ca2+, PS, and damage proxies. This is a placeholder release architecture pending subtype-specific evidence.

## Layer 6: Cargo sorting, composition, potency
Cargo enrichment and potency score are linear proxies from Ca2+, ROS, and ATP. Literature is needed to map stimulation history to protein, RNA, lipid, and functional potency outputs.

## Layer 7: Injury, debris, quality gate
Damage accumulation is phenomenological, with Hill-based apoptosis and necrosis fractions plus debris-based purity. Placeholder thresholds should be replaced by evidence-supported injury metrics.

## Layer 8: Manufacturing, isolation, QC
Isolation efficiency, purity factor, batch consistency, and scalability are static coefficients. Workflow-specific isolation and QC evidence is needed.

## Cross-cutting cell state
Cell-state modifiers scale calcium handling, baseline EV release, and stress sensitivity. Future work should map these modifiers to curated cell-type and disease-state evidence.

## Overall status
Version 0.1.0 is structurally complete but scientifically provisional. Every mechanistic layer still requires targeted literature review and parameter replacement.
