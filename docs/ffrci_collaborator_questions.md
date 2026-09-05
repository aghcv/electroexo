# Experimental Metadata Needed for Electro Exocytosis Model Calibration

Thank you for sharing the Exoid particle measurements and RNA sequencing results. We can use the data to test and calibrate parts of the electro exocytosis model, but several experimental details are needed to connect the file labels to physical exposure conditions and to define the correct statistical units. The questions marked essential should be resolved before any parameter fitting. A sample manifest or protocol file is preferable to narrative answers where one already exists.

The current Exoid export shows a strong late decrease for `5p40kV`: total particle concentration falls by approximately 57% from 1 to 3 hours and by 58% from 30 minutes to 3 hours. That decrease appears in all three broad size ranges and all three p-labeled distributions. The `1p20kV` and `3p40kV` trajectories are nonmonotonic. We need the answers below to determine whether these patterns reflect extracellular particle kinetics, sample collection, or data processing.

## Immediate Questions From the Concentration Drop Audit

1. Is summing the concentration values across all reported size bins the correct way to recover the total concentration for each sample? If not, please provide the instrument-reported total and explain the units of each size-bin value.
2. Within every stimulated condition and harvest time, the summed total is numerically identical across the available `p1`, `p2`, and `p3` records, even though their size distributions differ. Were the distributions normalized to a shared total, averaged or rescaled during export, or derived from repeated scans of one pooled sample?
3. Are `p1`, `p2`, and `p3` independent biological cultures, technical runs, or processed distributions from the same sample? Does a given p label identify the same biological culture across 30 minutes, 1 hour, and 3 hours?
4. Could you share the uncombined raw Exoid export for every biological sample, including the instrument-reported total concentration, dilution factor, accepted event count, measurement volume, and quality-control flags?
5. Were the 30 minute, 1 hour, and 3 hour samples separate destructive harvests, or were they serially withdrawn from the same culture? Did the medium accumulate from exposure until each harvest, or was it replaced so that the samples represent collection intervals?
6. What volume of medium was present and collected at each time? If serial sampling was used, how much volume was removed and replaced, and were concentrations corrected for dilution, evaporation, or prior withdrawals?
7. Were sham cultures collected at 30 minutes, 1 hour, and 3 hours? The file contains `Sham2` and `sham media` but no sham time labels. Please map each control to its harvest time and treated comparison, or confirm that the controls are single-time measurements.
8. What is the experimental distinction between `Sham2` and `sham media`? Their exported all-size concentrations differ by approximately 1.85-fold, so they should not be combined until the protocol difference is known.
9. Were viable-cell count, viability, apoptosis or necrosis, LDH release, and debris measured at each harvest, especially for `5p40kV`? Please provide the measurements per culture rather than only a group mean if available.
10. Were samples measured in randomized order, and were storage time, freeze-thaw cycles, dilution, cartridge, and instrument settings balanced across condition and harvest time? Please identify any reruns, exclusions, out-of-range measurements, or background-subtraction steps.

## Essential Items Before Fitting

1. Does `20kV` or `40kV` in a sample label mean generator voltage, electrode voltage measured at the chamber, or electric field strength in kV per centimeter?
2. What pulse width, waveform, polarity, repetition rate, electrode gap, and exposure geometry were used for each condition?
3. Are `p1`, `p2`, and `p3` independent biological replicates, technical measurement replicates, separate donors, or repeated measurements from one pooled sample? In particular, does `p1` identify the same biological source across 30 minutes, 1 hour, and 3 hours, or do the replicate labels restart independently at each harvest?
4. Were the 30 minute, 1 hour, and 3 hour measurements collected from separate cultures or by serial sampling from the same culture? If they were separate destructive harvests, were cultures matched by donor, isolation batch, or another blocking variable? Did conditioned medium accumulate from exposure until harvest, or was the medium replaced between collection intervals?
5. Why is the summed particle concentration numerically indistinguishable across the available replicates within every condition and time, with a maximum relative spread of approximately `9e-11`, while the size distributions differ? Does the Exoid export contain a condition-level concentration applied to separate replicate distributions? Please share the uncombined, unrounded per-sample totals if they exist.
6. Which control should be compared with each treated sample, and what is the difference between `Sham2` and `sham media`?
7. What do `Treatment 1` and `Treatment 2` mean in the RNA sample names, and when was RNA collected relative to pulse exposure?
8. Are the reported RNA `count` columns raw integer read counts, estimated counts, averaged counts, or transformed values? The shared values include fractions.
9. Which RNA groups are `GroupA` and `GroupB` in the differential-expression table, and does positive `logFC` mean higher expression after 3p40kV treatment?

## Sample Identity and Experimental Units

1. Could you provide a sample manifest with one row per culture or sample and the following columns where applicable?
   - Unique sample identifier
   - Donor identifier
   - Cell type and source
   - Treatment condition
   - Pulse settings
   - Collection time
   - Control type
   - Biological replicate identifier
   - Technical replicate identifier
   - Culture batch and processing batch
   - Exoid file label
   - RNA sample name
2. Were the CD4 positive cells obtained from one donor, multiple donors, or pooled donors?
3. If multiple donors were used, were conditions paired within donor? Which samples came from each donor?
4. What do the replicate labels `p1`, `p2`, and `p3` represent, and can the same label be followed across harvest times as one matched sample set?
5. Does the missing `3p40kV 1h CD4+ p3` dataset reflect a failed measurement, excluded sample, or omitted export?
6. Do the Exoid samples and RNA samples come from the same cultures, matched parallel cultures, or separate experiments?
7. Please provide the CD4 isolation method, purity, activation state, passage or culture age, seeding density, cell concentration, culture volume, and time in culture before exposure.
8. We understand that five million cells were used per batch. Does this mean five million starting cells in every treated and sham culture? What exact conditioned-medium volume corresponds to each Exoid concentration, and did either cell number or volume differ across conditions or harvest times?

## Nanosecond Pulse Exposure

1. For each condition, please report the programmed and measured pulse amplitude and state whether it is voltage or field strength.
2. What was the electrode spacing, and how was field strength calculated or measured at the sample?
3. What were the pulse width at half maximum, rise time, fall time, waveform, polarity, and pulse-to-pulse interval?
4. What repetition rate was used? Were pulses delivered as one train or multiple bursts?
5. What exposure vessel and electrode geometry were used, including the sample volume and electrode material?
6. Was the waveform verified with an oscilloscope at the load? If available, could you share representative voltage and current traces?
7. What were the exposure-medium composition, conductivity, osmolarity, pH, extracellular calcium concentration, and starting temperature?
8. Was temperature measured during or immediately after exposure? If so, what was the maximum temperature rise for each condition?
9. Was the sham handled in the same vessel for the same duration with electrodes present and the pulse generator disabled?

## EV Collection Time Course

1. Please confirm that `30min`, `1h`, and `3h` denote time after the nsPEF exposure ended.
2. Were these destructive harvests from separate wells or serial withdrawals from the same well? If separate wells were used, which metadata identify wells that belong to the same donor or culture preparation across time?
3. Did EV-containing medium accumulate continuously from time zero to each harvest, or do the samples represent intervals such as 0 to 30 minutes, 30 to 60 minutes, and 1 to 3 hours?
4. Was medium replaced immediately before exposure, immediately after exposure, or between harvests?
5. What volume was collected at each time, and was that volume corrected for prior sampling or evaporation?
6. Were sham controls collected at each time? If yes, which sham label maps to each harvest time?
7. Were viable-cell count, total-cell count, viability, apoptosis, necrosis, or LDH release measured at each collection time?
8. Should particle release be normalized per milliliter, per starting cell number, per viable cell number at harvest, or per viable-cell hour?
9. Were cells allowed to reattach or recover after exposure, and were any cells removed with the collected medium?

## Exoid Particle Measurements

1. What is the exact Exoid instrument model, cartridge or pore type, software version, and analysis version?
2. What calibration beads, concentration standards, pressure, voltage, detection threshold, and accepted size range were used?
3. How many instrument runs were made per biological sample, and were their outputs averaged before export?
4. What sample dilution was used, and are concentrations already back-calculated to the original conditioned-medium volume?
5. Does each size-bin value represent concentration in that bin, concentration density per nanometer, or a smoothed distribution? Is summing the bins the correct way to recover total particles per milliliter?
6. Could you share the Exoid total concentration reported for each individual sample and, if available, the uncombined raw export?
7. Why do all available replicates within a condition and time have numerically indistinguishable summed concentrations? Are the reported decimals copied from a shared condition-level estimate, normalized to a shared total, or independently measured?
8. Were background, blank-medium, buffer, or cartridge counts subtracted? Please describe the subtraction rule.
9. What were the lower and upper limits of detection and quantification? How were zero or missing bins handled?
10. What preprocessing occurred before measurement, including centrifugation speeds and durations, filtration, concentration, size-exclusion chromatography, ultracentrifugation, precipitation, or other isolation steps?
11. How long and at what temperature were samples stored before measurement? How many freeze-thaw cycles occurred?
12. Were particle identity and co-isolates assessed using EV markers, negative markers, protein-to-particle ratio, microscopy, detergent lysis, or an orthogonal particle method?

## RNA Sequencing Design and Processing

1. Was RNA extracted from the producer cells, extracellular vesicles, conditioned medium, or another fraction?
2. What was the RNA collection time after exposure for each sample?
3. What exposure settings correspond to `Treatment_1` and `Treatment_2` for CD4 cells and HUVECs?
4. What are the units of `EGTA_0.125`, `0.25_EGTA`, and `0.5_EGTA`? When was EGTA added and removed relative to exposure?
5. What was the total extracellular calcium concentration in each EGTA condition, and was free calcium calculated or measured?
6. Are the three samples in each named group biological replicates? Do they correspond to donors or batches that should be included as blocking factors?
7. Which genome build and gene annotation release were used?
8. What library preparation, sequencing platform, read layout, read length, and target depth were used?
9. Which aligner or quantifier produced the gene-level values, and how were multimapping reads and low-expression genes handled?
10. Could you provide the raw integer gene-count matrix, sample design table, and quality-control report, including mapping rate, library size, duplication, RNA integrity, and sample-level PCA?
11. Which differential-expression package, version, normalization, dispersion method, design formula, covariates, contrasts, and filtering threshold were used?
12. Were any samples excluded? Several nonprimary groups have much lower reported count totals than their replicates, so the exclusion and QC decisions would help us avoid propagating outliers.
13. Which gene-set collection and version were used for enrichment, and what ranking statistic and permutation scheme were used?
14. Some `TAG_%` entries in the pathway file appear as dates such as `11-May` or `13-Oct`. Could you share the original enrichment output before spreadsheet date conversion?

## Measurements That Can Anchor Model Layers

1. Are calcium, ROS, ATP, mitochondrial potential, membrane permeability, phosphatidylserine exposure, or membrane-resealing measurements available for any of these exposures?
2. Are viability, apoptosis, necrosis, cell count, or debris measurements available at 30 minutes, 1 hour, 3 hours, or later?
3. Were EV-associated markers or size-fractionated measurements collected that could distinguish small EV enriched material from medium or large EVs and injury-associated particles?
4. Are total EV protein, particle-to-protein ratio, cargo measurements, or functional potency assays available?
5. Is there an unexposed baseline sample at time zero for particle release and RNA expression?

## Preferred Files

If available, the most useful next files are a sample manifest, the complete pulse and exposure protocol, uncombined Exoid exports with per-sample total concentrations, cell-count and viability measurements, the raw integer RNA count matrix, the RNA design table, and the RNA quality-control report. These files would let us distinguish measured inputs from parameters that require literature priors and from parameters that can be estimated using the experiment.
