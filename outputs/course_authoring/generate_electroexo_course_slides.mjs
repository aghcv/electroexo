import fs from "node:fs/promises";
import path from "node:path";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const COURSE_TITLE = "Computational Systems Biology: Electro-Exocytosis";
const COURSE_SUBTITLE =
  "From compartment models and PK/PD analogies to the electroexo multiscale framework";

const SLIDE_SIZE = { width: 1280, height: 720 };
const FRAME = { left: 52, top: 46, width: 1176, height: 628 };
const COLORS = {
  bg: "#FFFFFF",
  text: "#111827",
  muted: "#4B5563",
  faint: "#6B7280",
  panel: "#F3F4F6",
  panel2: "#E5E7EB",
  accent: "#0F4C81",
  accentSoft: "#DCE9F5",
  accentDeep: "#0B3558",
  success: "#0F766E",
  warning: "#9A6700",
};

const FIG_DIR = "/Users/aghorban/code/electro-exocytosis/figs";

const IMAGE_PATHS = {
  framework: path.join(FIG_DIR, "fig01.png"),
  dosimetrySchematic: path.join(FIG_DIR, "figA1_dosimetry_schematic.png"),
  temperatureProfiles: path.join(FIG_DIR, "temperature_profiles.png"),
  membraneVoltage: path.join(FIG_DIR, "membrane_voltage_by_pulse_width.png"),
  layer3: path.join(FIG_DIR, "layer3_timeseries_comparison.png"),
  remodeling: path.join(FIG_DIR, "remodeling_repair_timeseries.png"),
  evBiogenesis: path.join(FIG_DIR, "ev_biogenesis_timeseries.png"),
  storyboard: path.join(FIG_DIR, "multilayer_storyboard.png"),
};

const LAYER_LABELS = [
  "1 Pulse and dosimetry",
  "2 Electrodynamics",
  "3 Ion transport and stress",
  "4 Remodeling and repair",
  "5 EV biogenesis and release",
  "6 Cargo and potency",
  "7 Injury and quality gate",
  "8 Manufacturing and QC",
];

const DECKS = [
  {
    slug: "00_Course_Overview_and_Syllabus",
    title: "Course Overview and Syllabus",
    subtitle: COURSE_SUBTITLE,
    sessionLabel: "Overview",
    slides: [
      {
        kind: "title",
        title: COURSE_TITLE,
        subtitle:
          "A senior undergraduate and graduate module course built around the electroexo codebase and the electro-exocytosis manuscript.",
        kicker: "Overview deck",
      },
      {
        kind: "bullets",
        title: "Audience, premise, and outcomes",
        bullets: [
          "Students start with basic physiology, transport, and differential equations, then climb toward a multiscale research framework.",
          "The course treats electro-exocytosis as a systems-biology problem rather than a single-pathway story.",
          "By the end, students should be able to read, critique, and extend the layered electroexo model.",
          "Every session ties biological abstractions to concrete computational modules, states, and assumptions.",
        ],
        calloutTitle: "Main repositories",
        calloutBody:
          "electroexo package repository\nelectro-exocytosis manuscript repository",
      },
      {
        kind: "table",
        title: "Ten-session course arc",
        headers: ["Session", "Focus", "Framework destination"],
        rows: [
          ["1", "Systems-biology framing", "Why electro-exocytosis needs a layered model"],
          ["2", "Compartment models and PK/PD", "Mass balances, states, and effect mappings"],
          ["3", "Cell biology foundations", "Organelles, calcium, exocytosis, EV routes"],
          ["4", "nsPEF biophysics and dosimetry", "Pulse descriptors and exposure metrics"],
          ["5", "ODEs and nonlinear systems", "Coupled dynamics and multiscale numerics"],
          ["6", "Layers 1-2 in electroexo", "Pulse, dosimetry, electrodynamics"],
          ["7", "Layer 3 in electroexo", "Ion transport, ROS, ATP, mitochondria"],
          ["8", "Layers 4-5 in electroexo", "Repair logic and EV release"],
          ["9", "Layers 6-8 and modifiers", "Cargo, quality, manufacturing, phenotype"],
          ["10", "Calibration and extension", "Evidence workflow and capstone ideas"],
        ],
        note:
          "The progression is intentionally cumulative: every prerequisite topic earns its place by making a later electroexo submodule readable.",
      },
      {
        kind: "module-map",
        title: "The framework we are teaching toward",
        activeLayers: [1, 2, 3, 4, 5, 6, 7, 8],
        body:
          "The package resolves nanosecond pulse descriptors first, then integrates a reduced ODE system over longer biological timescales. Cell state acts across every layer instead of appearing as one downstream box.",
      },
      {
        kind: "bullets",
        title: "How the course uses the project assets",
        bullets: [
          "The manuscript defines the biological rationale, module map, and preliminary figures.",
          "The electroexo package shows how those ideas were operationalized into code, parameters, examples, tests, and outputs.",
          "The evidence workbook and calibration targets teach students how literature becomes model structure and parameter constraints.",
          "Assignments can ask students to inspect assumptions, run scenarios, or propose better submodels without pretending the framework is already validated.",
        ],
        calloutTitle: "Important framing",
        calloutBody:
          "Version 0.1.0 is structurally complete but scientifically provisional. That makes it ideal for teaching model design, critique, and extension.",
      },
      {
        kind: "readings",
        title: "Core reading spine",
        primaryReadings: [
          "electro-exocytosis manuscript: framework overview, methods, module tables",
          "README and docs/model_assumptions.md in electroexo",
          "Kotnik et al. 2019 on membrane electroporation mechanisms and models",
        ],
        secondaryReadings: [
          "Hodgkin and Huxley 1952 for state-based biophysical modeling",
          "De Young and Keizer 1992 for compartmental calcium dynamics",
          "Hucka et al. 2003 for model representation and systems-biology context",
        ],
        activity:
          "Opening assignment: draw a one-page concept map from pulse protocol to EV engineering output, then annotate which links already have code support in electroexo.",
      },
    ],
  },
  {
    slug: "01_Systems_Biology_Framing_and_Electro_Exocytosis",
    title: "Session 1: Systems-Biology Framing and Electro-Exocytosis",
    subtitle: "Why this is a coupled multiscale problem",
    sessionLabel: "Session 1",
    slides: [
      {
        kind: "title",
        title: "Systems-Biology Framing and Electro-Exocytosis",
        subtitle:
          "We begin by defining the biological question, the engineering objective, and the modeling scope.",
        kicker: "Session 1",
      },
      {
        kind: "bullets",
        title: "Why electro-exocytosis is not a one-equation problem",
        bullets: [
          "An nsPEF protocol perturbs membranes on nanosecond scales but EV release unfolds over seconds to hours.",
          "The outputs that matter are multidimensional: yield, subtype balance, cargo, potency, purity, and viability.",
          "The same pulse can support productive remodeling in one regime and injury-dominant shedding in another.",
          "A useful model must connect electrical forcing to intracellular signaling, repair, trafficking, and manufacturing-facing outputs.",
        ],
        calloutTitle: "Teaching goal",
        calloutBody:
          "Students should learn to ask: what level of abstraction preserves mechanism without demanding an impossible whole-cell simulation?",
      },
      {
        kind: "bullets",
        title: "The central systems-biology idea",
        bullets: [
          "Represent the pulse as an external forcing function rather than simulating the power supply in full detail.",
          "Represent organelles and pathways as compartments, states, and fluxes with interpretable interfaces.",
          "Use state variables to carry memory across scales: permeability, calcium, ROS, ATP, repair state, EV pools.",
          "Use outputs that match experiments and engineering decisions, not only hidden internal variables.",
        ],
        calloutTitle: "Framework language",
        calloutBody:
          "Inputs -> states -> couplings -> outputs -> evidence checks",
      },
      {
        kind: "module-map",
        title: "Framework map introduced on day one",
        activeLayers: [1, 2, 3, 4, 5, 7, 8],
        body:
          "Students do not need every submodule on day one, but they should see the full map early so later prerequisites feel purposeful instead of disconnected.",
      },
      {
        kind: "figure",
        title: "Integrated mechanism and engineering view",
        imagePath: IMAGE_PATHS.framework,
        imageFit: "contain",
        notes: [
          "The manuscript figure is ideal for orienting students to the pathway logic before any equations appear.",
          "Use it to distinguish membrane perturbation, calcium redistribution, repair, EV routes, and final engineering outputs.",
          "Emphasize that the model aims to separate useful EV production from nonspecific stress debris.",
        ],
        sourceLabel: "Source figure: electro-exocytosis manuscript / figs/fig01.png",
      },
      {
        kind: "readings",
        title: "Readings and teaching use",
        primaryReadings: [
          "electro-exocytosis manuscript: Introduction and Results framework table",
          "README.md in electroexo",
          "Kotnik et al. 2019, Annual Review of Biophysics",
        ],
        secondaryReadings: [
          "Beebe et al. 2003 on nanosecond pulsed fields",
          "Manoochehri et al. 2025 on EV engineering context",
        ],
        activity:
          "In class: students identify three places where experimental observables naturally attach to the model and three places where placeholder assumptions are still obvious.",
      },
    ],
  },
  {
    slug: "02_Compartment_Models_and_PK_PD_Analogies",
    title: "Session 2: Compartment Models and PK/PD Analogies",
    subtitle: "Mass balances, states, fluxes, thresholds, and effect mappings",
    sessionLabel: "Session 2",
    slides: [
      {
        kind: "title",
        title: "Compartment Models and PK/PD Analogies",
        subtitle:
          "We use familiar exposure-response ideas to prepare students for reduced cell-state modeling.",
        kicker: "Session 2",
      },
      {
        kind: "bullets",
        title: "Compartmental thinking from first principles",
        bullets: [
          "A compartment is an abstraction that stores a quantity, exchanges it with other compartments, and evolves over time.",
          "The core modeling grammar is mass balance: rate of change = inflows - outflows + sources - sinks.",
          "Compartments can be spatial, biochemical, or phenomenological, as long as their interfaces are explicit.",
          "This is exactly the mindset needed for cytosol, ER, mitochondria, MVB pools, and damage states.",
        ],
        calloutTitle: "Bridge to electroexo",
        calloutBody:
          "Layer 3 and Layer 5 are both compartment models; they simply track different biological objects.",
      },
      {
        kind: "table",
        title: "PK/PD analogies that help students enter the framework",
        headers: ["PK/PD idea", "Cell-systems analogy", "Where it appears in electroexo"],
        rows: [
          ["Dose", "Pulse train and energy deposition", "pulse.py, dosimetry.py"],
          ["Exposure compartment", "Membrane or organelle state", "electrodynamics.py"],
          ["Concentration-time profile", "Ca2+, ROS, ATP trajectories", "ion_transport.py"],
          ["Effect model", "Repair activation or EV release gate", "remodeling_repair.py, ev_release.py"],
          ["Therapeutic window", "Productive secretory window", "storyboard and scenario comparisons"],
        ],
        note:
          "The analogy is pedagogical, not literal. It helps students see why reduced state models can still be mechanistically useful.",
      },
      {
        kind: "bullets",
        title: "Ordinary differential equations as the default language",
        bullets: [
          "Once states and fluxes are defined, each compartment becomes an ODE for its time evolution.",
          "Thresholds and saturation enter through Hill functions, logistic gates, or bounded rate laws.",
          "The same structure can represent calcium clearance, organelle recovery, docking turnover, or damage accumulation.",
          "Students should become comfortable reading biology as a collection of interacting dynamical balances.",
        ],
        calloutTitle: "Canonical template",
        calloutBody:
          "dx/dt = inflow(x,u) - outflow(x,u) + coupling(x,z,u)",
      },
      {
        kind: "module-map",
        title: "Where compartment modeling lives in this framework",
        activeLayers: [3, 5, 6, 7, 8],
        body:
          "The reduced cell model is not one giant equation. It is a linked family of smaller compartment models whose outputs become each other's inputs.",
      },
      {
        kind: "readings",
        title: "Readings and student exercise",
        primaryReadings: [
          "Hodgkin and Huxley 1952",
          "De Young and Keizer 1992",
          "Slepchenko et al. 2003 on quantitative cell biology with Virtual Cell",
        ],
        secondaryReadings: [
          "Smith, Wagner, and Keizer 1996 on calcium buffering approximations",
          "Mosley et al. 2013 and Ying et al. 2013 on vesicle modeling analogies",
        ],
        activity:
          "Homework: write a two-compartment calcium model and explain how it would need to change before it could serve as a precursor to electroexo Layer 3.",
      },
    ],
  },
  {
    slug: "03_Cell_Biology_Foundations_for_the_Framework",
    title: "Session 3: Cell Biology Foundations for the Framework",
    subtitle: "Membranes, organelles, calcium microdomains, exocytosis, and EV routes",
    sessionLabel: "Session 3",
    slides: [
      {
        kind: "title",
        title: "Cell Biology Foundations for the Framework",
        subtitle:
          "Students need a compact but rigorous map of the biology before the computational modules feel intuitive.",
        kicker: "Session 3",
      },
      {
        kind: "bullets",
        title: "Membrane and organelle essentials",
        bullets: [
          "The plasma membrane is only one electrically responsive interface; nsPEF can also perturb ER, mitochondria, and endosomal structures.",
          "Each organelle matters because it stores gradients, constrains traffic, or participates in recovery and stress.",
          "In the framework, organelles are represented by the states they contribute, not by geometric perfection.",
          "Students should be able to explain why ER and mitochondria cannot be collapsed into generic 'cell stress' too early.",
        ],
        calloutTitle: "Key organelles",
        calloutBody:
          "ER: calcium store\nMitochondria: ATP and ROS\nMVB/endosomes: exosome route\nPlasma membrane: poration and repair interface",
      },
      {
        kind: "bullets",
        title: "Calcium microdomains and buffering matter",
        bullets: [
          "A small average cytosolic calcium rise can hide much larger local calcium near pores, channels, or fusion machinery.",
          "Buffers, pumps, exchangers, and organelle uptake shape both the amplitude and the persistence of the signal.",
          "This is why Layer 4 uses a submembrane calcium proxy instead of relying only on bulk Ca2+.",
          "Students should distinguish calcium as messenger, stressor, and repair trigger.",
        ],
        calloutTitle: "Relevant references",
        calloutBody:
          "Neher and Augustine 1992\nFaas et al. 2011\nSmith et al. 1996",
      },
      {
        kind: "bullets",
        title: "EV biogenesis routes to keep straight",
        bullets: [
          "Small EV output is linked to endosomal maturation, ILV loading, docking, and fusion of MVB-related compartments.",
          "Medium/large EV output is more closely tied to budding, actomyosin state, and membrane remodeling.",
          "Apoptotic bodies belong to a distinct injury-associated regime and should not be interpreted as productive secretion.",
          "The framework is useful because it predicts route competition, not just total particle count.",
        ],
        calloutTitle: "Vocabulary to stabilize early",
        calloutBody:
          "MVB, ILV, budding precursor, apoptotic commitment, secretory bias, lysosomal routing",
      },
      {
        kind: "module-map",
        title: "Biology-heavy regions of the framework",
        activeLayers: [3, 4, 5, 6],
        body:
          "This session supplies the biological meaning for the state variables students will later see in ion_transport.py, remodeling_repair.py, and ev_release.py.",
      },
      {
        kind: "readings",
        title: "Readings and discussion prompt",
        primaryReadings: [
          "Futter et al. 1996 on MVB maturation and lysosomal fusion",
          "Neher and Augustine 1992",
          "Faas et al. 2011",
        ],
        secondaryReadings: [
          "Muratori et al. 2017 and 2021 on repair and phosphatidylserine exposure",
          "Williams et al. 2023 and 2025 on repair-coupled vesicle release",
        ],
        activity:
          "In class: students classify which biological processes should be states, which should be parameters, and which should remain external inputs in a reduced model.",
      },
    ],
  },
  {
    slug: "04_nsPEF_Biophysics_and_Dosimetry",
    title: "Session 4: nsPEF Biophysics and Dosimetry",
    subtitle: "Pulse descriptors, exposure geometry, heating, and field effects",
    sessionLabel: "Session 4",
    slides: [
      {
        kind: "title",
        title: "nsPEF Biophysics and Dosimetry",
        subtitle:
          "Before students can understand Layer 1 and Layer 2, they need a compact language for pulse delivery and exposure metrics.",
        kicker: "Session 4",
      },
      {
        kind: "bullets",
        title: "Pulse protocol vocabulary to normalize",
        bullets: [
          "Amplitude, pulse width, pulse number, repetition rate, waveform, and train duration all matter.",
          "Equivalent 'dose' claims can hide very different electrical and thermal histories.",
          "Pulse width is not just a label; in this framework it changes charging dynamics and downstream access to organelles.",
          "Geometry and conductivity determine how nominal protocol settings translate into effective exposure.",
        ],
        calloutTitle: "Code anchors",
        calloutBody:
          "electro_exocytosis/models/pulse.py\nelectro_exocytosis/models/dosimetry.py",
      },
      {
        kind: "bullets",
        title: "The reduced Layer 1 logic",
        bullets: [
          "pulse.py converts protocol inputs into peak field, pulse duration, train duration, energy density, and a dose index.",
          "dosimetry.py then applies geometry factors, Joule heat density, and a thermal-retention interpretation.",
          "This is intentionally simpler than a full finite-element solver, but it preserves exposure bookkeeping.",
          "Students should learn to ask what data are needed before replacing these placeholders with exposure-specific models.",
        ],
        calloutTitle: "Teaching formula",
        calloutBody:
          "u ~ sigma * E^2 * pulse_width * pulse_number * waveform_factor",
      },
      {
        kind: "figure",
        title: "Dosimetry schematic from the manuscript",
        imagePath: IMAGE_PATHS.dosimetrySchematic,
        imageFit: "contain",
        notes: [
          "Use this to distinguish pulse descriptor inputs from thermal and exposure outputs.",
          "Students should see how geometry, conductivity, and repetition rate alter interpretation even before biology starts.",
        ],
        sourceLabel: "Source figure: electro-exocytosis manuscript / figs/figA1_dosimetry_schematic.png",
      },
      {
        kind: "figure",
        title: "Preliminary temperature-profile behavior",
        imagePath: IMAGE_PATHS.temperatureProfiles,
        imageFit: "contain",
        notes: [
          "This figure supports the idea that retained heating is not identical to absorbed electrical dose.",
          "It also creates a natural transition into why Layer 1 exposes multiple dosimetry interpretations rather than one scalar.",
        ],
        sourceLabel: "Source figure: electro-exocytosis manuscript / figs/temperature_profiles.png",
      },
      {
        kind: "readings",
        title: "Readings and classroom use",
        primaryReadings: [
          "Beebe et al. 2003",
          "Kotnik and Miklavcic 2000",
          "Kotnik et al. 2019",
        ],
        secondaryReadings: [
          "Orlacchio et al. 2023 for absorbed energy and temperature targets",
          "Krassowska and Filev 2007 for electroporation modeling context",
        ],
        activity:
          "Mini-lab: students compare two protocols with similar nominal amplitude but different pulse widths and explain why the downstream model should not treat them as equivalent.",
      },
    ],
  },
  {
    slug: "05_ODEs_Nonlinearity_and_Multiscale_Numerics",
    title: "Session 5: ODEs, Nonlinearity, and Multiscale Numerics",
    subtitle: "How the framework moves from biology to a solvable coupled system",
    sessionLabel: "Session 5",
    slides: [
      {
        kind: "title",
        title: "ODEs, Nonlinearity, and Multiscale Numerics",
        subtitle:
          "This session equips students to read the simulation pipeline rather than treating it as a black box.",
        kicker: "Session 5",
      },
      {
        kind: "bullets",
        title: "Why ODE systems are the workhorse here",
        bullets: [
          "Each biological layer contributes states, rates, and couplings that can be written as first-order time-evolution rules.",
          "ODEs are expressive enough for thresholds, saturation, recovery, accumulation, and bounded state variables.",
          "They are also transparent enough to support teaching, debugging, and staged calibration.",
          "The point is not to claim cellular completeness; it is to preserve interpretable mechanism across timescales.",
        ],
        calloutTitle: "Numerical teaching target",
        calloutBody:
          "Students should be able to trace where every derivative term comes from and what biological story it encodes.",
      },
      {
        kind: "bullets",
        title: "Nonlinear building blocks students must recognize",
        bullets: [
          "Logistic gates encode threshold-like activation such as membrane permeability or organelle response.",
          "Hill functions encode cooperative or saturating responses such as calcium sensing and apoptotic commitment.",
          "Bounded recovery terms keep repair, ROS, ATP, and viability dynamics interpretable.",
          "Composed nonlinearities create regime behavior, which is why sensitivity and identifiability become important later.",
        ],
        calloutTitle: "Common motifs",
        calloutBody:
          "logistic(deltaV)\nHill(Ca)\nfirst-order recovery\nstress-dependent penalties",
      },
      {
        kind: "bullets",
        title: "How electroexo assembles the full simulation",
        bullets: [
          "simulation.py builds pulse descriptors, dosimetry, and electrodynamic states before constructing the longer-timescale ODE right-hand side.",
          "Ion states and EV states are integrated together, while remodeling observables are recomputed as couplers.",
          "A multiscale scheduler chooses the output time grid so nanosecond inputs can still drive minute-scale trajectories.",
          "The resulting outputs are stored as standard-language tables, summaries, plots, and metadata bundles.",
        ],
        calloutTitle: "Code anchors",
        calloutBody:
          "electro_exocytosis/simulation.py\nelectro_exocytosis/numerics/multiscale.py\nelectro_exocytosis/numerics/solvers.py",
      },
      {
        kind: "module-map",
        title: "Where coupling and numerics sit in the framework",
        activeLayers: [1, 2, 3, 4, 5],
        body:
          "The numerics layer is not biologically separate. It is the glue that keeps layered abstractions coherent when they operate on very different timescales.",
      },
      {
        kind: "readings",
        title: "Readings and coding exercise",
        primaryReadings: [
          "Hodgkin and Huxley 1952",
          "Hucka et al. 2003",
          "electroexo simulation.py and tests/test_solver_smoke.py",
        ],
        secondaryReadings: [
          "Krassowska and Filev 2007",
          "Kotnik et al. 2019",
        ],
        activity:
          "Homework: identify three nonlinearities in the current framework, explain why each was chosen, and describe what experimental evidence would justify replacing it.",
      },
    ],
  },
  {
    slug: "06_Layers_1_and_2_Pulse_to_Electrodynamics",
    title: "Session 6: Layers 1-2, Pulse to Electrodynamics",
    subtitle: "From pulse descriptors to membrane and organelle perturbation states",
    sessionLabel: "Session 6",
    slides: [
      {
        kind: "title",
        title: "Layers 1-2: Pulse to Electrodynamics",
        subtitle:
          "Students now read the actual implementation of the first two layers and connect it to the manuscript logic.",
        kicker: "Session 6",
      },
      {
        kind: "bullets",
        title: "Implementation walkthrough for Layers 1 and 2",
        bullets: [
          "pulse.py computes peak field, duration, dose index, train duration, and absorbed-energy proxies.",
          "dosimetry.py adds geometry and thermal-retention interpretations, keeping the exposure model selectable.",
          "electrodynamics.py converts effective field into a Schwan-style membrane voltage with nanosecond charging effects.",
          "Reduced permeability and pore-density proxies become the bridge into the ion-transport layer.",
        ],
        calloutTitle: "Key classes and outputs",
        calloutBody:
          "PulseDescriptors\nDosimetryResult\nElectrodynamicsState",
      },
      {
        kind: "double-figure",
        title: "Two figures students should learn to interpret together",
        leftImagePath: IMAGE_PATHS.temperatureProfiles,
        rightImagePath: IMAGE_PATHS.membraneVoltage,
        leftLabel: "Layer 1: retained heating and dose interpretation",
        rightLabel: "Layer 2: charging and voltage-threshold effects",
        notes: [
          "These figures teach that pulse width and conductivity are not bookkeeping details; they materially shape the upstream state seen by later layers.",
          "They also justify separating pulse, dosimetry, and electrodynamics into distinct modules instead of collapsing them into one black box.",
        ],
      },
      {
        kind: "table",
        title: "Calibration targets that matter first for Layers 1 and 2",
        headers: ["Target paper", "Useful data", "Model parameters affected"],
        rows: [
          ["Orlacchio 2023", "absorbed energy density, temperature, viability", "waveform and thermal retention factors"],
          ["Pakhomov 2007", "conductance increase and recovery", "pore reseal and permeability scales"],
          ["Nesin 2012", "leak current and channel inhibition", "leak and VG-channel modifiers"],
          ["Baker 2024", "organelle electroporation geometry", "organelle voltage fractions and geometry factors"],
        ],
        note:
          "This is a good teaching point: students should see that 'more literature' only helps when the paper constrains a real model quantity.",
      },
      {
        kind: "module-map",
        title: "Framework focus for this session",
        activeLayers: [1, 2],
        body:
          "Layer 1 defines the controllable engineering input. Layer 2 defines the immediate cell-level perturbation state that everything downstream must interpret.",
      },
      {
        kind: "readings",
        title: "Readings and student task",
        primaryReadings: [
          "pulse.py, dosimetry.py, electrodynamics.py",
          "Kotnik et al. 2019",
          "Krassowska and Filev 2007",
        ],
        secondaryReadings: [
          "Pakhomov et al. 2007",
          "Nesin et al. 2012",
          "Orlacchio et al. 2023",
        ],
        activity:
          "Code reading task: students annotate which parts of Layers 1-2 are structural, which are placeholder parameterizations, and which would likely need cell-type-specific adaptation.",
      },
    ],
  },
  {
    slug: "07_Layer_3_Ion_Transport_and_Bioenergetics",
    title: "Session 7: Layer 3, Ion Transport and Bioenergetics",
    subtitle: "Calcium, ER release, mitochondria, ROS, ATP, and osmotic stress",
    sessionLabel: "Session 7",
    slides: [
      {
        kind: "title",
        title: "Layer 3: Ion Transport and Bioenergetics",
        subtitle:
          "This is the first major ODE layer and the main bridge from electroporation to longer-lived biology.",
        kicker: "Session 7",
      },
      {
        kind: "bullets",
        title: "State vector and flux logic in Layer 3",
        bullets: [
          "The current state vector tracks cytosolic, ER, and mitochondrial calcium, membrane potential, ROS, ATP, Na, K, Cl, and osmotic stress.",
          "Fluxes include pore-mediated calcium entry, ER release, SERCA refilling, PMCA and NCX clearance, and mitochondrial uptake/release.",
          "Additional terms convert ionic perturbation into ROS generation, ATP depletion, depolarization, and osmotic burden.",
          "Students should notice that this layer is already a compact systems-biology module on its own.",
        ],
        calloutTitle: "Code anchors",
        calloutBody:
          "electro_exocytosis/models/ion_transport.py\nexamples/compare_ion_transport_bioenergetics.py",
      },
      {
        kind: "figure",
        title: "Preliminary Layer 3 behavior from the manuscript figures",
        imagePath: IMAGE_PATHS.layer3,
        imageFit: "contain",
        notes: [
          "The comparison plot is useful because it shows how different assumptions generate distinct calcium, ROS, mitochondrial, and osmotic histories.",
          "Students should relate these trajectories back to named fluxes and timescales rather than reading them only as curves.",
        ],
        sourceLabel: "Source figure: electro-exocytosis manuscript / figs/layer3_timeseries_comparison.png",
      },
      {
        kind: "table",
        title: "High-priority calibration targets for Layer 3",
        headers: ["Source", "Observable", "Why it matters"],
        rows: [
          ["Semenov 2013", "Ca2+ thresholds and slopes", "sets pore and ER-release gains"],
          ["Bagalkot 2018 / Yun 2024", "pulse-duration-dependent Ca2+ responses", "separates pore and VGCC contributions"],
          ["Nuccitelli 2013", "pulse-number ROS dependence", "constrains ROS source terms"],
          ["Radzeviciute-Valciuke 2024", "ATP depletion under calcium electroporation", "bounds stress and damage coupling"],
        ],
        note:
          "Layer 3 is where experimental observables begin to map very naturally onto state trajectories and rate constants.",
      },
      {
        kind: "module-map",
        title: "Framework focus for this session",
        activeLayers: [3],
        body:
          "Layer 3 translates a brief permeability event into longer-lived signaling and metabolic consequences, which is why so many downstream modules depend on it.",
      },
      {
        kind: "readings",
        title: "Readings and computational exercise",
        primaryReadings: [
          "ion_transport.py",
          "De Young and Keizer 1992",
          "Semenov et al. 2013",
        ],
        secondaryReadings: [
          "Neher and Augustine 1992",
          "Bagalkot et al. 2018",
          "Nuccitelli et al. 2013",
        ],
        activity:
          "Homework: choose one Layer 3 flux and propose a more evidence-backed alternative, including what experimental data would be needed to fit it.",
      },
    ],
  },
  {
    slug: "08_Layers_4_and_5_Repair_and_EV_Release",
    title: "Session 8: Layers 4-5, Repair and EV Release",
    subtitle: "From calcium-dependent remodeling to subtype-resolved vesicle output",
    sessionLabel: "Session 8",
    slides: [
      {
        kind: "title",
        title: "Layers 4-5: Repair and EV Release",
        subtitle:
          "Students now study how the framework converts upstream stress into repair-compatible or injury-compatible vesicle behavior.",
        kicker: "Session 8",
      },
      {
        kind: "bullets",
        title: "What Layer 4 adds beyond calcium alone",
        bullets: [
          "Layer 4 converts bulk and local calcium into PS exposure, calpain activity, annexin recruitment, lysosomal repair, actomyosin tension, and repair state.",
          "This layer decides whether the cell interprets the perturbation as recoverable remodeling or escalating damage.",
          "The submembrane calcium proxy is pedagogically important because it links microdomain biology to reduced modeling.",
          "Repair-associated shedding becomes a mechanistic input to the EV module rather than an afterthought.",
        ],
        calloutTitle: "Code anchor",
        calloutBody:
          "electro_exocytosis/models/remodeling_repair.py",
      },
      {
        kind: "bullets",
        title: "What Layer 5 adds beyond 'EV count'",
        bullets: [
          "The model tracks MVB pool, ILV load, docked MVB pool, budding pool, and apoptotic commitment.",
          "It separates secretory bias from lysosomal routing and distinguishes small EV, medium/large EV, and apoptotic-body outputs.",
          "This creates a route-competition picture instead of a single release-rate scalar.",
          "Students should see why this is essential for an EV-engineering framing.",
        ],
        calloutTitle: "Code anchor",
        calloutBody:
          "electro_exocytosis/models/ev_release.py\nexamples/compare_ev_biogenesis_release.py",
      },
      {
        kind: "double-figure",
        title: "Repair-state and EV-biogenesis figures for class discussion",
        leftImagePath: IMAGE_PATHS.remodeling,
        rightImagePath: IMAGE_PATHS.evBiogenesis,
        leftLabel: "Layer 4: repair and remodeling trajectories",
        rightLabel: "Layer 5: EV pool and release trajectories",
        notes: [
          "These two figures should be taught together so students can watch how upstream repair logic reshapes downstream route selection.",
          "A productive secretory window is not just 'more EV'; it is a pattern of repair support, route bias, and controlled stress.",
        ],
      },
      {
        kind: "module-map",
        title: "Framework focus for this session",
        activeLayers: [4, 5],
        body:
          "This is the pivot from cellular response to vesicle engineering. Layer 4 decides the membrane-repair interpretation; Layer 5 decides how that interpretation becomes release behavior.",
      },
      {
        kind: "readings",
        title: "Readings and student exercise",
        primaryReadings: [
          "remodeling_repair.py and ev_release.py",
          "Muratori et al. 2017 and 2021",
          "Williams et al. 2023 and 2025",
        ],
        secondaryReadings: [
          "Bhattacharya et al. 2022",
          "Hellwich et al. 2026",
        ],
        activity:
          "In class: students debate which Layer 5 states should remain explicit in a minimal model and which could be collapsed without losing the engineering interpretation.",
      },
    ],
  },
  {
    slug: "09_Layers_6_to_8_and_Cell_State_Modifiers",
    title: "Session 9: Layers 6-8 and Cell-State Modifiers",
    subtitle: "Cargo, potency, injury, QC, manufacturing, and phenotype adaptation",
    sessionLabel: "Session 9",
    slides: [
      {
        kind: "title",
        title: "Layers 6-8 and Cell-State Modifiers",
        subtitle:
          "The framework now expands from vesicle release to engineering-relevant product interpretation and phenotype dependence.",
        kicker: "Session 9",
      },
      {
        kind: "bullets",
        title: "Why the framework does not stop at EV release",
        bullets: [
          "A translational EV model must say more than 'release increased'. It must address cargo, potency, purity, and producer-cell viability.",
          "Layer 6 currently uses placeholder cargo-enrichment proxies linked to Ca2+, ROS, and ATP.",
          "Layer 7 translates accumulated damage into apoptosis, necrosis, purity, and a pass/fail quality gate.",
          "Layer 8 translates biological output into isolation efficiency, purity factor, batch consistency, and scalability.",
        ],
        calloutTitle: "Code anchors",
        calloutBody:
          "cargo_potency.py\ninjury_quality.py\nmanufacturing_qc.py",
      },
      {
        kind: "bullets",
        title: "Cell state as a cross-cutting modifier layer",
        bullets: [
          "electroexo treats phenotype adaptation as modifiers on calcium handling, baseline EV release, and stress sensitivity.",
          "This is a teaching-friendly compromise between one universal model and a separate model for every cell type.",
          "Students should recognize both the usefulness and the limits of scalar modifier strategies.",
          "This is also the right place to discuss disease state, activation state, and heterogeneity as modeling challenges.",
        ],
        calloutTitle: "Code anchor",
        calloutBody:
          "electro_exocytosis/models/cell_state.py",
      },
      {
        kind: "figure",
        title: "Integrated storyboard across productive and injurious regimes",
        imagePath: IMAGE_PATHS.storyboard,
        imageFit: "contain",
        notes: [
          "This figure is ideal for showing how the same layered framework supports regime interpretation from mild reversible to injury-dominant windows.",
          "Use it to connect yield, subtype balance, repair state, and viability in one integrated classroom discussion.",
        ],
        sourceLabel: "Source figure: electro-exocytosis manuscript / figs/multilayer_storyboard.png",
      },
      {
        kind: "module-map",
        title: "Framework focus for this session",
        activeLayers: [6, 7, 8],
        body:
          "These layers convert mechanism into engineering meaning. They also force students to confront how much of the framework is still placeholder and where future evidence is most needed.",
      },
      {
        kind: "readings",
        title: "Readings and project prompt",
        primaryReadings: [
          "cargo_potency.py, injury_quality.py, manufacturing_qc.py, cell_state.py",
          "README placeholder-status section",
          "docs/model_assumptions.md",
        ],
        secondaryReadings: [
          "Manoochehri et al. 2025 for application context",
          "Relevant MISEV/EV nomenclature background from the manuscript bibliography",
        ],
        activity:
          "Project prompt: propose a more defensible Layer 6 or Layer 8 submodel, explain what evidence exists already, and identify what data are still missing.",
      },
    ],
  },
  {
    slug: "10_Calibration_Evidence_and_Extending_the_Framework",
    title: "Session 10: Calibration, Evidence, and Extending the Framework",
    subtitle: "How students move from reading the framework to improving it",
    sessionLabel: "Session 10",
    slides: [
      {
        kind: "title",
        title: "Calibration, Evidence, and Extending the Framework",
        subtitle:
          "The last session teaches how computational systems biology remains honest: literature curation, staged calibration, and uncertainty-aware extension.",
        kicker: "Session 10",
      },
      {
        kind: "bullets",
        title: "The evidence workflow built into the repository",
        bullets: [
          "The evidence workbook and calibration_targets.csv are not decorative; they are the bridge between literature and parameter replacement.",
          "EvidenceLoader exposes module maps, literature tracking, and calibration targets in a form students can query and critique.",
          "A good calibration target names a measurable observable and a concrete model quantity it constrains.",
          "This keeps students from confusing narrative support with fit-ready evidence.",
        ],
        calloutTitle: "Code anchors",
        calloutBody:
          "electro_exocytosis/evidence/evidence_loader.py\nelectro_exocytosis/evidence/calibration_targets.csv\ndocs/fulltext_calibration_opportunities.md",
      },
      {
        kind: "table",
        title: "A staged calibration hierarchy for the course",
        headers: ["Stage", "What to fit first", "Why that order matters"],
        rows: [
          ["1", "Pulse, dosimetry, permeability, resealing", "stabilizes the exposure-to-perturbation map"],
          ["2", "Ca2+, ROS, ATP, osmotic stress", "constrains the main signaling layer before downstream tuning"],
          ["3", "Repair and shedding observables", "separates recoverable remodeling from injury"],
          ["4", "Subtype-resolved EV release", "fits productive versus injurious output regimes"],
          ["5", "Cargo, potency, QC, manufacturing", "adds translational interpretation after upstream states are credible"],
        ],
        note:
          "This staged order is one of the most important conceptual lessons in the course: do not fit downstream outputs before upstream state logic is constrained.",
      },
      {
        kind: "bullets",
        title: "Software workflow students should leave able to use",
        bullets: [
          "Define scenarios in YAML and understand what belongs in the scenario versus parameter overrides.",
          "Run simulations from the CLI and interpret summaries, time-series tables, and plot bundles.",
          "Trace assumptions from docs/model_assumptions.md back into specific code modules and tests.",
          "Use examples and storyboard figures to compare regimes instead of staring only at raw numbers.",
        ],
        calloutTitle: "Concrete files",
        calloutBody:
          "docs/input_schema.md\nelectro_exocytosis/cli.py\nexamples/*.py\ntests/*.py",
      },
      {
        kind: "bullets",
        title: "Capstone directions and open problems",
        bullets: [
          "Replace one placeholder rate law with a literature-backed alternative and justify the required parameters.",
          "Add cell-type-specific modifier logic for a chosen producer cell class.",
          "Design a calibration study that fits Layer 3 or Layer 4 against published assays before touching EV output.",
          "Propose a direct-EV-engineering branch that shares pulse physics but diverges in state variables and QC outputs.",
        ],
        calloutTitle: "Final teaching message",
        calloutBody:
          "Good systems-biology modeling is not only about solving equations. It is about choosing abstractions that can grow honestly with evidence.",
      },
      {
        kind: "readings",
        title: "Closing readings and final assignment",
        primaryReadings: [
          "docs/fulltext_calibration_opportunities.md",
          "evidence_loader.py and calibration_targets.csv",
          "electro-exocytosis manuscript appendix framework tables",
        ],
        secondaryReadings: [
          "README disclaimer and placeholder sections",
          "Any session-specific primary paper chosen for the final project",
        ],
        activity:
          "Final assignment: submit a proposed framework extension with one biological motivation, one computational form, one calibration plan, and one validation limitation.",
      },
    ],
  },
];

async function ensureDir(dir) {
  await fs.mkdir(dir, { recursive: true });
}

function addTextBox(slide, opts) {
  const shape = slide.shapes.add({
    geometry: "textbox",
    position: {
      left: opts.left,
      top: opts.top,
      width: opts.width,
      height: opts.height,
    },
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  shape.text = opts.text;
  shape.text.style = {
    fontSize: opts.fontSize ?? 18,
    bold: opts.bold ?? false,
    color: opts.color ?? COLORS.text,
    alignment: opts.alignment ?? "left",
  };
  return shape;
}

function addPanel(slide, { left, top, width, height, fill = COLORS.panel, lineFill = COLORS.panel2, radius = "rounded-xl" }) {
  return slide.shapes.add({
    geometry: "roundRect",
    position: { left, top, width, height },
    fill,
    line: { style: "solid", fill: lineFill, width: 1 },
    borderRadius: radius,
  });
}

function addHeader(slide, deck, slideTitle, slideIndex, slideCount) {
  slide.background.fill = COLORS.bg;
  slide.shapes.add({
    geometry: "rect",
    position: { left: FRAME.left, top: 20, width: FRAME.width, height: 6 },
    fill: COLORS.accent,
    line: { style: "solid", fill: COLORS.accent, width: 0 },
  });
  addTextBox(slide, {
    text: COURSE_TITLE,
    left: FRAME.left,
    top: 28,
    width: 520,
    height: 22,
    fontSize: 13,
    bold: true,
    color: COLORS.accentDeep,
  });
  addTextBox(slide, {
    text: deck.sessionLabel,
    left: FRAME.left + FRAME.width - 220,
    top: 28,
    width: 220,
    height: 22,
    fontSize: 13,
    bold: true,
    color: COLORS.faint,
    alignment: "right",
  });
  addTextBox(slide, {
    text: slideTitle,
    left: FRAME.left,
    top: 58,
    width: 900,
    height: 44,
    fontSize: 34,
    bold: true,
    color: COLORS.text,
  });
  addTextBox(slide, {
    text: `${slideIndex + 1} / ${slideCount}`,
    left: FRAME.left + FRAME.width - 90,
    top: 684,
    width: 90,
    height: 20,
    fontSize: 12,
    color: COLORS.faint,
    alignment: "right",
  });
}

function addBulletList(slide, bullets, { left, top, width, height, fontSize = 21, color = COLORS.text }) {
  const text = bullets.map((bullet) => `• ${bullet}`).join("\n\n");
  addTextBox(slide, {
    text,
    left,
    top,
    width,
    height,
    fontSize,
    color,
  });
}

function addCallout(slide, { title, body, left, top, width, height }) {
  addPanel(slide, { left, top, width, height, fill: COLORS.accentSoft, lineFill: COLORS.accentSoft });
  addTextBox(slide, {
    text: title,
    left: left + 18,
    top: top + 16,
    width: width - 36,
    height: 24,
    fontSize: 21,
    bold: true,
    color: COLORS.accentDeep,
  });
  addTextBox(slide, {
    text: body,
    left: left + 18,
    top: top + 52,
    width: width - 36,
    height: height - 64,
    fontSize: 18,
    color: COLORS.accentDeep,
  });
}

function addModuleBoxes(slide, activeLayers, left, top, width) {
  const gap = 14;
  const cols = 4;
  const boxWidth = (width - gap * (cols - 1)) / cols;
  const boxHeight = 86;
  LAYER_LABELS.forEach((label, index) => {
    const row = Math.floor(index / cols);
    const col = index % cols;
    const x = left + col * (boxWidth + gap);
    const y = top + row * (boxHeight + gap);
    const active = activeLayers.includes(index + 1);
    addPanel(slide, {
      left: x,
      top: y,
      width: boxWidth,
      height: boxHeight,
      fill: active ? COLORS.accent : COLORS.panel,
      lineFill: active ? COLORS.accent : COLORS.panel2,
      radius: "rounded-lg",
    });
    addTextBox(slide, {
      text: label,
      left: x + 12,
      top: y + 14,
      width: boxWidth - 24,
      height: boxHeight - 28,
      fontSize: 18,
      bold: true,
      color: active ? "#FFFFFF" : COLORS.text,
    });
  });
  addPanel(slide, {
    left,
    top: top + 2 * (boxHeight + gap) + 8,
    width,
    height: 52,
    fill: COLORS.accentSoft,
    lineFill: COLORS.accentSoft,
    radius: "rounded-lg",
  });
  addTextBox(slide, {
    text:
      "Cross-cutting modifier layer: cell type, disease state, calcium handling, baseline EV release, and stress sensitivity tune every mechanistic block.",
    left: left + 16,
    top: top + 2 * (boxHeight + gap) + 20,
    width: width - 32,
    height: 28,
    fontSize: 17,
    color: COLORS.accentDeep,
  });
}

async function addImage(slide, imagePath, position, fit = "contain") {
  const bytes = await fs.readFile(imagePath);
  slide.images.add({
    blob: bytes,
    contentType: "image/png",
    alt: path.basename(imagePath),
    fit,
    position,
  });
}

function styleTable(table, rows, columns) {
  table.styleOptions = { headerRow: true, bandedRows: true };
  table.borders.assign({ style: "solid", fill: COLORS.panel2, width: 1 });
  for (let col = 0; col < columns; col += 1) {
    const cell = table.getCell(0, col);
    cell.fill = COLORS.accent;
    cell.text.style = { fontSize: 15, bold: true, color: "#FFFFFF" };
  }
  for (let row = 1; row < rows; row += 1) {
    for (let col = 0; col < columns; col += 1) {
      table.getCell(row, col).text.style = { fontSize: 15, color: COLORS.text };
    }
  }
}

async function renderSlide(slide, deck, slideSpec, slideIndex, slideCount) {
  switch (slideSpec.kind) {
    case "title": {
      slide.background.fill = COLORS.bg;
      slide.shapes.add({
        geometry: "rect",
        position: { left: 0, top: 0, width: 1280, height: 720 },
        fill: COLORS.bg,
        line: { style: "solid", fill: COLORS.bg, width: 0 },
      });
      slide.shapes.add({
        geometry: "rect",
        position: { left: 844, top: 0, width: 436, height: 720 },
        fill: COLORS.accentSoft,
        line: { style: "solid", fill: COLORS.accentSoft, width: 0 },
      });
      slide.shapes.add({
        geometry: "rect",
        position: { left: 844, top: 0, width: 18, height: 720 },
        fill: COLORS.accent,
        line: { style: "solid", fill: COLORS.accent, width: 0 },
      });
      addTextBox(slide, {
        text: slideSpec.kicker ?? deck.sessionLabel,
        left: 72,
        top: 76,
        width: 260,
        height: 24,
        fontSize: 18,
        bold: true,
        color: COLORS.accent,
      });
      addTextBox(slide, {
        text: slideSpec.title,
        left: 72,
        top: 128,
        width: 680,
        height: 180,
        fontSize: 50,
        bold: true,
        color: COLORS.text,
      });
      addTextBox(slide, {
        text: slideSpec.subtitle,
        left: 72,
        top: 332,
        width: 700,
        height: 150,
        fontSize: 24,
        color: COLORS.muted,
      });
      addTextBox(slide, {
        text: `${COURSE_TITLE}\n${deck.subtitle}`,
        left: 888,
        top: 112,
        width: 310,
        height: 180,
        fontSize: 24,
        bold: true,
        color: COLORS.accentDeep,
      });
      addTextBox(slide, {
        text:
          "Built from the live electroexo package, the electro-exocytosis manuscript, and the current evidence-driven module structure.",
        left: 888,
        top: 340,
        width: 300,
        height: 200,
        fontSize: 18,
        color: COLORS.accentDeep,
      });
      addTextBox(slide, {
        text: `${slideIndex + 1} / ${slideCount}`,
        left: 1120,
        top: 676,
        width: 92,
        height: 20,
        fontSize: 12,
        color: COLORS.faint,
        alignment: "right",
      });
      break;
    }
    case "bullets": {
      addHeader(slide, deck, slideSpec.title, slideIndex, slideCount);
      addBulletList(slide, slideSpec.bullets, {
        left: FRAME.left,
        top: 132,
        width: 716,
        height: 470,
        fontSize: 20,
      });
      addCallout(slide, {
        title: slideSpec.calloutTitle,
        body: slideSpec.calloutBody,
        left: 812,
        top: 146,
        width: 380,
        height: 360,
      });
      break;
    }
    case "module-map": {
      addHeader(slide, deck, slideSpec.title, slideIndex, slideCount);
      addModuleBoxes(slide, slideSpec.activeLayers, FRAME.left, 150, FRAME.width);
      addTextBox(slide, {
        text: slideSpec.body,
        left: FRAME.left,
        top: 520,
        width: FRAME.width,
        height: 90,
        fontSize: 20,
        color: COLORS.muted,
      });
      break;
    }
    case "figure": {
      addHeader(slide, deck, slideSpec.title, slideIndex, slideCount);
      addPanel(slide, { left: FRAME.left, top: 138, width: 780, height: 448, fill: "#FAFAFA" });
      await addImage(slide, slideSpec.imagePath, { left: FRAME.left + 18, top: 154, width: 744, height: 412 }, slideSpec.imageFit);
      addCallout(slide, {
        title: "What to notice",
        body: slideSpec.notes.map((note) => `• ${note}`).join("\n\n"),
        left: 858,
        top: 148,
        width: 334,
        height: 336,
      });
      if (slideSpec.sourceLabel) {
        addTextBox(slide, {
          text: slideSpec.sourceLabel,
          left: FRAME.left,
          top: 602,
          width: 820,
          height: 18,
          fontSize: 12,
          color: COLORS.faint,
        });
      }
      break;
    }
    case "double-figure": {
      addHeader(slide, deck, slideSpec.title, slideIndex, slideCount);
      addPanel(slide, { left: FRAME.left, top: 138, width: 560, height: 284, fill: "#FAFAFA" });
      addPanel(slide, { left: 668, top: 138, width: 560, height: 284, fill: "#FAFAFA" });
      await addImage(slide, slideSpec.leftImagePath, { left: FRAME.left + 14, top: 150, width: 532, height: 252 }, "contain");
      await addImage(slide, slideSpec.rightImagePath, { left: 682, top: 150, width: 532, height: 252 }, "contain");
      addTextBox(slide, {
        text: slideSpec.leftLabel,
        left: FRAME.left,
        top: 432,
        width: 560,
        height: 24,
        fontSize: 18,
        bold: true,
        color: COLORS.text,
      });
      addTextBox(slide, {
        text: slideSpec.rightLabel,
        left: 668,
        top: 432,
        width: 560,
        height: 24,
        fontSize: 18,
        bold: true,
        color: COLORS.text,
      });
      addPanel(slide, { left: FRAME.left, top: 478, width: FRAME.width, height: 120, fill: COLORS.accentSoft, lineFill: COLORS.accentSoft });
      addTextBox(slide, {
        text: slideSpec.notes.map((note) => `• ${note}`).join("\n\n"),
        left: FRAME.left + 18,
        top: 496,
        width: FRAME.width - 36,
        height: 90,
        fontSize: 18,
        color: COLORS.accentDeep,
      });
      break;
    }
    case "table": {
      addHeader(slide, deck, slideSpec.title, slideIndex, slideCount);
      const values = [slideSpec.headers, ...slideSpec.rows];
      const table = slide.tables.add({
        rows: values.length,
        columns: slideSpec.headers.length,
        left: FRAME.left,
        top: 142,
        width: FRAME.width,
        height: 404,
        values,
      });
      styleTable(table, values.length, slideSpec.headers.length);
      addPanel(slide, {
        left: FRAME.left,
        top: 566,
        width: FRAME.width,
        height: 48,
        fill: COLORS.accentSoft,
        lineFill: COLORS.accentSoft,
      });
      addTextBox(slide, {
        text: slideSpec.note,
        left: FRAME.left + 14,
        top: 580,
        width: FRAME.width - 28,
        height: 22,
        fontSize: 16,
        color: COLORS.accentDeep,
      });
      break;
    }
    case "readings": {
      addHeader(slide, deck, slideSpec.title, slideIndex, slideCount);
      addPanel(slide, { left: FRAME.left, top: 144, width: 550, height: 318, fill: COLORS.panel });
      addPanel(slide, { left: 678, top: 144, width: 550, height: 318, fill: COLORS.panel });
      addTextBox(slide, {
        text: "Primary readings",
        left: FRAME.left + 16,
        top: 160,
        width: 250,
        height: 24,
        fontSize: 22,
        bold: true,
        color: COLORS.text,
      });
      addBulletList(slide, slideSpec.primaryReadings, {
        left: FRAME.left + 16,
        top: 198,
        width: 516,
        height: 246,
        fontSize: 18,
      });
      addTextBox(slide, {
        text: "Secondary readings",
        left: 694,
        top: 160,
        width: 250,
        height: 24,
        fontSize: 22,
        bold: true,
        color: COLORS.text,
      });
      addBulletList(slide, slideSpec.secondaryReadings, {
        left: 694,
        top: 198,
        width: 516,
        height: 246,
        fontSize: 18,
      });
      addCallout(slide, {
        title: "Activity or assignment",
        body: slideSpec.activity,
        left: FRAME.left,
        top: 490,
        width: FRAME.width,
        height: 108,
      });
      break;
    }
    default:
      throw new Error(`Unsupported slide kind: ${slideSpec.kind}`);
  }
}

async function exportDeck(deck, outdir, previewDir) {
  const presentation = Presentation.create({ slideSize: SLIDE_SIZE });
  for (let index = 0; index < deck.slides.length; index += 1) {
    const slide = presentation.slides.add();
    await renderSlide(slide, deck, deck.slides[index], index, deck.slides.length);
  }

  const deckPath = path.join(outdir, `${deck.slug}.pptx`);
  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(deckPath);

  const montageBlob = await presentation.export({ format: "webp", montage: true, scale: 1 });
  await fs.writeFile(path.join(previewDir, `${deck.slug}-montage.webp`), new Uint8Array(await montageBlob.arrayBuffer()));

  const titlePreview = await presentation.export({
    slide: presentation.slides.items[0],
    format: "png",
    scale: 1,
  });
  await fs.writeFile(path.join(previewDir, `${deck.slug}-slide01.png`), new Uint8Array(await titlePreview.arrayBuffer()));

  return deckPath;
}

async function writeOutline(outdir, deckPaths) {
  const lines = [
    COURSE_TITLE,
    COURSE_SUBTITLE,
    "",
    "Generated decks:",
    ...deckPaths.map((deckPath, index) => `${index + 1}. ${path.basename(deckPath)}`),
    "",
    "Session list:",
    ...DECKS.map((deck, index) => {
      const label = index === 0 ? "Overview" : deck.title;
      return `${index + 1}. ${label}`;
    }),
    "",
    "Primary source anchors:",
    "- /Users/aghorban/code/electroexo/README.md",
    "- /Users/aghorban/code/electroexo/docs/model_assumptions.md",
    "- /Users/aghorban/code/electroexo/docs/fulltext_calibration_opportunities.md",
    "- /Users/aghorban/code/electro-exocytosis/main.tex",
    "- /Users/aghorban/code/electro-exocytosis/source.bib",
  ];
  await fs.writeFile(path.join(outdir, "course_outline.txt"), `${lines.join("\n")}\n`, "utf8");
}

async function main() {
  const outdir = process.argv[2];
  const previewDir = process.argv[3];
  if (!outdir || !previewDir) {
    throw new Error("Usage: node generate_electroexo_course_slides.mjs <outdir> <previewDir>");
  }
  await ensureDir(outdir);
  await ensureDir(previewDir);

  const deckPaths = [];
  for (const deck of DECKS) {
    deckPaths.push(await exportDeck(deck, outdir, previewDir));
  }
  await writeOutline(outdir, deckPaths);
  console.log(`Generated ${deckPaths.length} course decks in ${outdir}`);
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
