# Model assumptions

## Layer 1: Pulse delivery, exposure geometry, dosimetry
Current implementation uses simple field conversion, energy-density scaling, and fixed geometry correction factors. These are placeholders and need literature-backed exposure-specific dosimetry relations.

## Layer 2: Plasma membrane and organelle electrodynamics
A simplified Schwan-style transmembrane voltage estimate is used with fixed organelle fractions and sigmoid permeability. Placeholder support is needed for membrane charging, pore formation, and organelle-specific coupling.

## Layer 3: Ion transport, Ca2+ mobilization, ROS, bioenergetics
The ODE system now separates pore-mediated extracellular Ca2+ influx, ER release and SERCA refilling, PMCA/NCX clearance, mitochondrial Ca2+ uptake/release, mitochondrial depolarization, Na+/K+/Cl- perturbation, osmotic stress, ROS generation, and ATP production/consumption. The structure follows the Table A5/A6 submodules, but the constants and closure relations remain provisional and should be calibrated against module-specific literature and experiments.

## Layer 4: Ca2+-dependent remodeling and repair
The reduced implementation now maps cytosolic Ca2+, pore activation, osmotic stress, and mitochondrial stress into a local Ca2+ microdomain proxy, scramblase/PS exposure, flippase suppression, calpain activity, annexin recruitment, lysosomal repair activity, actomyosin tension, actin disruption, resealing state, and repair-associated shedding. This follows the Table A7/A8 submodule structure, but the parameters remain provisional and should be calibrated against local Ca2+ imaging, annexin/ESCRT recruitment, PS externalization, calpain inhibition, cytoskeletal remodeling, and membrane-repair assays.

## Layer 5: EV biogenesis and subtype release
Layer 5 now uses a reduced pool-balance release architecture with explicit MVB pool, ILV load, docked-MVB, budding-precursor, and apoptotic-commitment states. These states are driven by Rab-conversion/maturation, ESCRT-dependent and ceramide-assisted ILV loading, secretory-versus-lysosomal routing, Ca2+-sensitive fusion, and budding/blebbing proxies. The structure is literature-backed and aligned to Table A9/A10, but the parameters remain provisional and should be calibrated against subtype-resolved secretion, trafficking, and inhibition data.

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

## Full-text calibration targets
The first full-text PDF pass identified experimental targets that map directly
to current placeholder parameters. See `docs/fulltext_calibration_opportunities.md`
for the rationale and `electro_exocytosis/evidence/calibration_targets.csv` for
the machine-readable target table. These targets should be used as priors or
constraints for exposure, electrodynamics, Ca2+, ROS/ATP, remodeling, and repair
states before fitting EV yield, cargo, and quality outputs to project-specific
in vitro data.
