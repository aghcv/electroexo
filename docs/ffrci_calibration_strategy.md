# Calibration Strategy for the FFRCI Electro Exocytosis Data

## Decision

The shared data can support a useful first calibration, but not a full fit of the current multilayer model. The Exoid file provides a three-time-point response surface for particle concentration and size distribution. The RNA data provide cellular-state evidence and an EGTA perturbation series that may help distinguish calcium-dependent from calcium-independent responses. The correct first target is a reduced observation and EV-release model with strong literature priors on upstream biology. Fitting all current parameters would be nonidentifiable and would produce a deceptively precise result.

Formal fitting should wait until the exposure settings, replicate structure, controls, collection protocol, and Exoid concentration export are confirmed. A collaborator-facing question set is available in `docs/ffrci_collaborator_questions.docx` and its editable Markdown source.

## Data Currently Available

The Exoid CSV contains 32 sample distributions: 26 treatment samples and 6 controls. The labeled treatment design is 1p20kV, 3p40kV, and 5p40kV at 0.5, 1, and 3 hours. Most cells have three labeled replicates; 3p40kV at 1 hour has two. The summed size-bin concentrations are shown below, pending confirmation that summing bins recovers total concentration.

| Labeled exposure | 0.5 hours | 1 hour | 3 hours |
|---|---:|---:|---:|
| 1p20kV | 6.80e9 | 4.59e9 | 8.06e9 |
| 3p40kV | 8.59e9 | 1.37e10 | 9.84e9 |
| 5p40kV | 6.45e9 | 6.33e9 | 2.70e9 |

The `Sham2` and `sham media` controls have summed concentrations of 6.85e9 and 3.71e9 particles per milliliter, respectively, but neither label contains a harvest time. Every replicate within a given condition and time has effectively the same summed concentration while retaining a different size distribution. Concentration therefore appears to be condition-level or normalized, rather than an independent replicate measurement.

The RNA matrix contains 78,986 gene rows and paired CPM and reported-count columns for 33 samples. It includes CD4 and HUVEC samples plus CD4 EGTA conditions. The provided CD4 sham versus 3p40kV table has 236 differential-expression rows at FDR below 0.05. Thirteen of 50 enrichment rows have FDR q-values below 0.05. The prominent enriched programs include TNFA signaling through NFKB, cholesterol homeostasis, hypoxia, apoptosis, IL2 STAT5 signaling, and p53. These signals support a stress and membrane-remodeling response, but do not directly measure a kinetic rate constant.

## Immediate Structural Mismatch

The current simulator reports cumulative EV release. A cumulative state cannot decline. The Exoid concentration can decline across independently harvested cultures because the measured supernatant pool reflects production, uptake, degradation, aggregation, loss during processing, changing volume, and the number of viable producer cells.

The repository now includes `ExperimentalObservationBridge`, which corrects the static population/volume unit mismatch and records cell count, volume, viability basis, recovery, dilution, and background assumptions. The three FFRI analyses have been repeated on a particle-equivalents-per-initial-cell basis and are summarized in `docs/ffrci_single_cell_reanalysis.md`. This bridge is necessary but is not the dynamic supernatant observation layer described below:

```text
dN_m/dt = viable_cells(t) * release_rate_m(t) - loss_rate_m * N_m(t)
measured_concentration_m = recovery_m * N_m(t) / collection_volume + background_m
```

Here `m` denotes a particle class or a reduced size-distribution component. Initially, use one total-particle pool and at most one broad size-shift component. Do not equate diameter bins with exosome, microvesicle, or apoptotic-body biogenesis. MISEV2023 recommends using EV terminology when biogenesis has not been demonstrated and emphasizes source quantity, processing, detection limits, and orthogonal characterization ([Welsh et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10850029/)).

## Parameter Importance Vectors

For parameter `i`, define an importance vector across experimental conditions and observables using local log elasticities:

```text
S_ij = d log y_j / d log theta_i
```

The vector should include predicted particle-pool concentration at each time, broad size-distribution summaries, viability, and any measured calcium, ROS, ATP, or repair endpoints. Scale each element by measurement uncertainty so a noisy observable does not dominate. The weighted sensitivity matrix then supports three decisions:

1. Rank individual parameters by the root mean square of their weighted elasticity vector.
2. Use singular vectors or eigenvectors of the sensitivity Gram matrix to identify parameter combinations that the design constrains.
3. Remove parameters with negligible influence on all measured outputs and combine parameters whose sensitivity vectors are nearly collinear.

The repository now contains a provisional local screen across 35 parameters. Under unverified assumptions of 60 ns square pulses, 1 Hz, cuvette geometry, 1.6 S/m, and interpreting 20 kV and 40 kV as kV/cm, the most influential EV-output parameters include `baseline_mlEV_rate`, `k_budding_s`, `baseline_sEV_rate`, `baseline_AB_rate`, and `k_ILV_release_s`. Important upstream state parameters include `tau_pore_reseal_s`, `K_scramblase_uM`, `microdomain_pore_gain`, and `PMCA_Vmax_uM_s`. These rankings describe the current equations around their defaults. They do not show that the data can identify those parameters.

After metadata correction, repeat the screen in two passes. Use Morris screening over broad prior ranges for all plausible parameters, followed by variance-based sensitivity on the reduced set. Sensitivity must be calculated for the actual measurement model, not only latent simulation states. Systems-biology models commonly constrain a few combinations much more strongly than individual parameters, so uncertainty in predictions is more informative than a single best-fit vector ([Gutenkunst et al. 2007](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.0030189)).

## Parameter Tiers

### Measured Inputs and Fixed Definitions

The pulse waveform, pulse width, pulse count, repetition rate, measured voltage or field, electrode spacing, sample geometry, medium conductivity, extracellular calcium, temperature, cell count, collection volume, and harvest time should be experimental inputs. Physical constants and unit conversions can be fixed. These values should not be inferred from EV output.

### Literature Informed Priors

Biological values from other cell types should generally be informative priors or bounds, not exact fixed constants. High-priority examples already represented in the repository evidence table include:

- `tau_pore_reseal_s` and conductance recovery from nsPEF patch-clamp studies.
- calcium-influx and ER-release thresholds and response slopes from 60 ns experiments ([Semenov et al. 2013](https://pubmed.ncbi.nlm.nih.gov/23220180/)).
- fast and slow membrane-repair time constants and calcium dependence. Bhattacharya and colleagues reported approximately 13 to 15 second fast and roughly 70 second slow components with calcium, with a slower component of roughly 110 seconds or more without calcium ([Bhattacharya et al. 2022](https://pubmed.ncbi.nlm.nih.gov/34838875/)).
- TMEM16F and phosphatidylserine dependence for the sign and plausible magnitude of scramblase-to-budding effects ([Muratori et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5702676/)).

These studies justify model structure and prior ranges. Their numerical values should not be transferred unchanged to primary CD4 cells because cell size, activation state, membrane composition, and calcium handling differ.

### Fit From the FFRCI Data

The first fit should contain no more than roughly four to six effective parameters, depending on the confirmed number of biological replicates. Recommended candidates are:

- one baseline particle-release scale, normalized by viable producer cells;
- one exposure-to-release gain shared across conditions;
- one release adaptation or exhaustion time constant;
- one extracellular loss or uptake time constant;
- one injury-dependent suppression term if viability data are available;
- one broad size-shift parameter if the Exoid size distributions are independently replicated.

The current model has several detailed EV rates. They should be grouped into effective parameters for this first fit. Individual `baseline_sEV_rate`, `baseline_mlEV_rate`, `k_ILV_release_s`, `k_budding_s`, and apoptotic-body parameters cannot be estimated separately from particle size alone.

### Hold Fixed or Exclude From This Fit

Cargo potency, purity, isolation efficiency, subtype-specific biogenesis rates, apoptotic-body rates, and manufacturing objective weights should remain fixed or be excluded unless matching assays are supplied. RNA expression should not be used to set an assay recovery factor or a physical pulse parameter.

## Bridging RNA Sequencing to Constitutive Parameters

The bridge should use a small number of latent biological programs rather than thousands of gene-specific coefficients:

```text
RNA counts -> pathway activity z_m -> module capacity modifier -> model parameter
log theta_i,c = log theta_i,reference + sum_m beta_i,m * z_m,c
```

Pre-exposure expression from the same cell population can inform constitutive parameter modifiers. Post-exposure expression is initially an output or validation target. Using post-exposure expression to set parameters and then claiming that the model predicts the same exposure response would leak outcome information into the inputs.

A practical first mapping is:

| RNA program or gene | Model interpretation | Candidate parameter or extension | Use now |
|---|---|---|---|
| `SLC8A2` and calcium-export program | NCX capacity | `NCX_Vmax_uM_s` modifier | Prior sign and relative modifier after timing is known |
| `ATP2B1` and calcium-pump program | plasma-membrane calcium clearance | `PMCA_Vmax_uM_s` modifier | Gene-set score from full matrix |
| `ATP2A2` and ER calcium program | SERCA capacity | `SERCA_Vmax_uM_s` modifier | Gene-set score from full matrix |
| `ANXA1`, annexins, and repair genes | calcium-dependent membrane repair | annexin capacity state and `resealing_annexin_weight` | Latent repair score, not single-gene fitting |
| `CAPN2` and calpain program | cytoskeletal remodeling and shedding | calpain capacity and `actin_calpain_weight` | Prior sign; needs protein or activity validation |
| `TMEM16F` program | phosphatidylserine scrambling | scramblase capacity and `K_scramblase_uM` | Prefer protein or PS assay for magnitude |
| `HMOX1`, `MT2A`, `GCLC`, `GPX4`, and `SOD2` | inducible antioxidant response | new antioxidant-capacity state affecting ROS decay | Use a pathway score; do not raise ROS production simply because response genes rise |
| `RAB27`, `RAB11`, `RAB35`, `TSG101`, `PDCD6IP`, and `SMPD3` programs | MVB docking, ESCRT loading, and ceramide pathways | grouped MVB and ILV capacity modifiers | Only after subtype evidence or stronger priors |
| apoptosis and p53 programs | injury response | apoptosis propensity and viability observation model | Validation or weak prior unless matched viability exists |
| cholesterol homeostasis | membrane composition and vesiculation | slow membrane-composition modifier | Direction requires protein or lipid evidence |

The shared significant-gene table includes direct overlaps such as `KCNN4`, `SLC8A2`, `ANXA1`, `MT2A`, and `HMOX1`. `KCNN4` represents a calcium-activated potassium channel not explicitly included in the present equations. It is better represented as a new channel branch than folded into nonspecific pore-mediated potassium flux.

The EGTA series may be the most informative transcriptomic component. If treatment, EGTA timing, free calcium, donor pairing, and RNA collection time are known, interaction contrasts can estimate which pathway scores depend on extracellular calcium. Those scores can constrain the signs and relative strengths of calcium-to-repair, calcium-to-stress, and calcium-to-release couplings.

For multi-hour experiments, introduce slow transcriptional states only after the static bridge is validated:

```text
dz_m/dt = activation_m(Ca, ROS, damage) - decay_m * z_m
theta_i(t) = theta_i,baseline * exp(sum_m beta_i,m * z_m(t))
```

Fit at the pathway level with sign-constrained, regularized coefficients. A time-course RNA experiment is needed to identify activation and decay constants. One post-exposure RNA time point can test direction and consistency but cannot determine transcriptional kinetics.

## Overfit Avoidance and Identifiability

1. Treat each biological culture or donor as an experimental unit. Do not treat Exoid size bins, technical reads, or instrument passes as independent samples.
2. Use a prespecified primary endpoint such as log particle concentration per million viable cell hours. Treat weighted median diameter or broad size fractions as secondary endpoints.
3. Compare the mechanistic model against simple null models: sham-only baseline, condition means, and a small time-by-dose regression. A mechanistic fit should improve held-out prediction, not merely training error.
4. Split validation by donor or biological replicate. If donors are unavailable, leave out an entire exposure condition or harvest time. Random row-level splits would leak correlated measurements.
5. Use literature-informed bounds and shrinkage priors. Report posterior or profile-likelihood intervals and parameter correlations, not only an optimizer result.
6. Run multistart optimization and simulation-based recovery tests. Parameters that cannot be recovered from synthetic data generated under this design should be fixed or combined before fitting real data.
7. Use profile likelihood to identify finite confidence intervals and practically nonidentifiable directions ([Raue et al. 2009](https://pubmed.ncbi.nlm.nih.gov/19505944/)).
8. Bootstrap cultures or donors, not size bins. Preserve paired designs and batches during resampling.
9. Perform posterior predictive checks for concentration, size distributions, viability, and RNA pathway scores. Inspect systematic residuals by time, exposure, donor, and batch.
10. Lock the parameter list and priors before examining the held-out data. Use the holdout once for model selection, then collect a new experiment for confirmation.

## Recommended Sequence

1. Resolve the collaborator questions and construct a machine-readable sample manifest.
2. Confirm per-sample Exoid concentration and whether the time points are cumulative or interval collections.
3. Add a supernatant-particle observation layer with viable-cell normalization and loss or uptake.
4. Re-run sensitivity screening using measured pulse and exposure inputs.
5. Calibrate upstream electrical, calcium, and repair parameters to literature data or matching assays, not to EV concentration alone.
6. Fit a four-to-six-parameter reduced EV observation model with blocked cross-validation.
7. Build pathway scores from raw RNA counts and test the EGTA interaction before introducing transcript-to-parameter modifiers.
8. Validate predictions on a withheld donor, exposure combination, or new experiment.

The most valuable next experiment would cross pulse number and field amplitude rather than confounding them, include time-matched sham controls, and measure particle concentration, viable cell number, and at least one upstream anchor such as calcium or phosphatidylserine exposure at the same time points.
