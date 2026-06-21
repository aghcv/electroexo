# Full-text calibration opportunities

This note summarizes the first pass through the full-text PDF corpus in
`/Users/aghorban/Downloads/electro-exocytosis/files`. The goal is not to claim
that the model is now calibrated. Instead, it identifies experiments that map
cleanly onto parameters already present in the computational framework.

The companion machine-readable table is
`electro_exocytosis/evidence/calibration_targets.csv`.

## Highest-priority calibration targets

| Layer | Target source | Useful reported data | Current model targets |
|---|---|---|---|
| A1 | Orlacchio et al. 2023 | Pulse duration, PRR, buffer conductivity, absorbed energy density, temperature, and viability in 3D spheroids | `waveform_energy_factor`, `thermal_retention_factor`, `dose_to_damage_scale` |
| A3 | Pakhomov et al. 2007 | Patch-clamp conductance increase and minute-scale recovery after 60 ns pulses | `pore_conductance_scale`, `tau_pore_reseal_s`, `pulse_number_permeability_scale` |
| A3 | Nesin et al. 2012 | nsPEF leak currents plus voltage-gated Ca2+/Na+ current inhibition after 300-600 ns pulses | `VGCC_activation_scale`, `ion_channel_inhibition_scale`, `leak_current_proxy` |
| A5 | Semenov et al. 2013 | External Ca2+ and ER-release thresholds, response slopes, and CICR threshold after 60 ns pulses | `Ca_pore_gain`, `ER_release_threshold`, `ER_release_gain`, `CICR_threshold_uM` |
| A5 | Bagalkot et al. 2018 and Yun et al. 2024 | Pulse-duration and field-strength effects on Ca2+ influx in chromaffin cells | `Ca_pore_gain`, `VGCC_activation_scale`, `pulse_width_Ca_scale` |
| A5 | Nuccitelli et al. 2013 | Pulse-number-dependent ROS response and Ca2+/antioxidant sensitivity | `ROS_pulse_factor_s`, `ROS_Ca_factor`, `tau_ROS_s` |
| A5 | Radzeviciute-Valciuke et al. 2024 | ATP depletion under nanosecond calcium electroporation protocols | `ATP_depletion_factor`, `ATP_damage_threshold`, `Ca_toxicity_scale` |
| A7 | Bhattacharya et al. 2022 | Fast and slow membrane resealing constants, Ca2+-dependent improvement of 90% resealing, dye-uptake reduction | `tau_repair_s`, `pore_reseal_fast_fraction`, `resealing_annexin_weight`, `resealing_lysosome_weight` |
| A7 | Muratori et al. 2017 and Silkunas et al. 2026 | TMEM16F/PS dependence, PS-puncta delay, site density, and Ca2+ dependence | `K_scramblase_uM`, `PS_max`, `PS_exposure_tau_s`, `PS_site_density` |
| A7 | Hellwich et al. 2026 | Ca2+-dependent actin remodeling and swelling after 200 x 300 ns pulses | `actin_calpain_weight`, `actin_osmotic_weight`, `microdomain_osmotic_gain` |
| A7/A9 | Williams et al. 2023 and 2025 | Ca2+-dependent exosome secretion and calpain/annexin microvesicle shedding during membrane repair | `K_sEV_Ca_uM`, `sEV_Ca_scale`, `K_calpain_uM`, `shedding_rate_scale` |

## How to use these data

The most immediate fitting targets are time constants, thresholds, and relative
response ratios rather than absolute EV yield. Examples include resealing
time constants from Bhattacharya et al., Ca2+ thresholds from Semenov et al.,
conductance recovery from Pakhomov et al., and PS/blebbing ratios from Muratori
et al. These can constrain the upstream signaling and repair states before the
model is asked to predict EV yield or cargo.

The review draft frames electro-exocytosis as a biomanufacturing problem: nsPEF
should be optimized for EV yield, subtype balance, cargo control, potency,
purity, and producer-cell viability. Therefore, full-text calibration should be
used as a staged hierarchy. First tune exposure, permeability, Ca2+, repair, and
injury states against literature assays. Then use the ODU group's in vitro EV
data to fit release, subtype, cargo, and quality-control outputs.

## What not to overfit yet

Many papers support inclusion of a submodule but do not provide enough
extractable numerical data for parameter fitting. Reviews and mechanistic
studies of MISEV, sphingolipids, annexin curvature, lysosomal repair, and
mitochondrial physiology are useful for model structure and manuscript
justification, but they should not be treated as direct calibration datasets.

Similarly, calcium electroporation and cancer-ablation protocols are valuable
for ATP, ROS, and injury boundaries, but they often operate near cytotoxic
conditions. Those data should tune stress and failure modes, not the preferred
non-lethal electro-exocytosis operating window.
