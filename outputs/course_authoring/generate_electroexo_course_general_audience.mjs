import fs from "node:fs/promises";
import path from "node:path";
import { execFileSync } from "node:child_process";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const SLIDE_SIZE = { width: 1280, height: 720 };
const FRAME = { left: 74, top: 94, width: 1132, height: 560 };
const COLORS = {
  bg: "#FFFFFF",
  teal: "#3B8599",
  tealDark: "#225C6A",
  text: "#111827",
  muted: "#4B5563",
  faint: "#6B7280",
  soft: "#E8F2F5",
  soft2: "#F4F7F8",
  accent: "#1D4ED8",
  red: "#B91C1C",
};

const ZIP_6931 =
  "/Users/aghorban/Library/CloudStorage/OneDrive-OldDominionUniversity/Papers/preparation/electro-exocytosis-comp-framework/presentation-material/BIEN6931-TopicsInBiomedicalEngineering/end of the semester back up/BIEN 6931 701 Topics in Biomedical Engr - 5132018 - 1103 PM.zip";
const PDF_5700_HOMEO =
  "/Users/aghorban/Library/CloudStorage/OneDrive-OldDominionUniversity/Papers/preparation/electro-exocytosis-comp-framework/presentation-material/BIEN5700-SystemPhysiology/01homeostasis.pdf";
const PDF_5700_SYLLABUS =
  "/Users/aghorban/Library/CloudStorage/OneDrive-OldDominionUniversity/Papers/preparation/electro-exocytosis-comp-framework/presentation-material/BIEN5700-SystemPhysiology/Syllabus sp18.pdf";

const FIG_ROOT = "/Users/aghorban/code/electro-exocytosis/figs";
const REPO_RESULTS = "/Users/aghorban/code/electroexo/results";

const FIGURES = {
  framework: path.join(FIG_ROOT, "fig01.png"),
  dosimetrySchematic: path.join(FIG_ROOT, "figA1_dosimetry_schematic.png"),
  temperatureProfiles: path.join(FIG_ROOT, "temperature_profiles.png"),
  membraneVoltage: path.join(FIG_ROOT, "membrane_voltage_by_pulse_width.png"),
  layer3: path.join(FIG_ROOT, "layer3_timeseries_comparison.png"),
  remodeling: path.join(FIG_ROOT, "remodeling_repair_timeseries.png"),
  evBiogenesis: path.join(FIG_ROOT, "ev_biogenesis_timeseries.png"),
  storyboard: path.join(FIG_ROOT, "multilayer_storyboard.png"),
  evRegimeMap: path.join(REPO_RESULTS, "ev_release_regime_map/ev_release_regime_map.png"),
  humanCellPulse: path.join(REPO_RESULTS, "human_cell_3d/storyboard_frames/01_01_nsPEF_exposure.png"),
  humanCellElectro: path.join(REPO_RESULTS, "human_cell_3d/storyboard_frames/02_02_membrane_electrodynamics.png"),
  humanCellCalcium: path.join(REPO_RESULTS, "human_cell_3d/storyboard_frames/03_03_Ca_influx_ER_release.png"),
  humanCellStress: path.join(REPO_RESULTS, "human_cell_3d/storyboard_frames/04_04_ionic_mito_ROS_ATP.png"),
  humanCellRepair: path.join(REPO_RESULTS, "human_cell_3d/storyboard_frames/05_05_remodeling_repair.png"),
  humanCellEV: path.join(REPO_RESULTS, "human_cell_3d/storyboard_frames/06_06_EV_release.png"),
};

const DECKS = [
  {
    slug: "00_Course_Overview_and_Roadmap",
    label: "Overview",
    title: "Course Overview and Roadmap",
    subtitle: "A general-audience route into computational systems biology and electro-exocytosis",
    slides: [
      {
        kind: "title",
        title: "Computational Systems Biology:\nElectro-Exocytosis",
        subtitle:
          "Designed for senior undergraduates who may know basic biology and engineering, but may know nothing about this project yet.",
      },
      {
        kind: "big-idea",
        title: "What this course is trying to do",
        statement:
          "Teach the background students need before they ever see the electroexo framework.",
        supporting:
          "We begin with homeostasis, transport, compartments, ODEs, PK/PD, electrophysiology, calcium biology, and exocytosis. Only then do we assemble the project framework.",
      },
      {
        kind: "bullets",
        title: "Audience and assumptions",
        bullets: [
          "Senior undergraduates in biomedical engineering or related areas.",
          "Comfort with calculus is helpful, but the biological meaning always comes first.",
          "No prior knowledge of extracellular vesicles, nsPEF, or the electroexo codebase is assumed.",
        ],
      },
      {
        kind: "table",
        title: "How the course is organized",
        headers: ["Block", "Main purpose", "Source inspiration"],
        rows: [
          ["Physiology block", "Build biological intuition", "BIEN5700 Systems Physiology"],
          ["Modeling block", "Build mathematical and computational intuition", "BIEN6931 Topics in Biomedical Engineering"],
          ["Electroexo block", "Integrate both into the project framework", "electroexo + electro-exocytosis"],
        ],
        note: "This separation keeps the course student-centered instead of project-centered.",
      },
      {
        kind: "table",
        title: "Session roadmap",
        headers: ["Session", "Main idea", "Why it matters later"],
        rows: [
          ["1", "Homeostasis and systems thinking", "Gives students a physiological lens"],
          ["2", "Transport, compartments, and mass balance", "Introduces the grammar of reduced models"],
          ["3", "ODEs and simulation logic", "Prepares students to read model equations"],
          ["4", "PK, PD, and PBPK", "Shows useful compartment-model patterns"],
          ["5", "Bioelectricity and electrophysiology", "Builds toward pulsed-field thinking"],
          ["6", "Calcium, organelles, and bioenergetics", "Builds the Layer 3 biological background"],
          ["7", "Exocytosis, repair, and extracellular vesicles", "Builds the Layer 4-5 biological background"],
          ["8", "nsPEF, electroporation, and dose", "Introduces the pulse perturbation itself"],
          ["9", "The electroexo framework", "Assembles the full layered model"],
          ["10", "Evidence, calibration, and projects", "Teaches honest model extension"],
        ],
        note: "The early sessions are prerequisites, not warm-up filler.",
      },
      {
        kind: "bullets",
        title: "Project-style expectations adapted from BIEN6931",
        bullets: [
          "Students should learn to state assumptions clearly.",
          "Students should connect equations to physical laws and biological meaning.",
          "Students should present figures, simulations, and limitations transparently.",
        ],
      },
    ],
  },
  {
    slug: "01_Homeostasis_and_Systems_Thinking",
    label: "Session 1",
    title: "Homeostasis and Systems Thinking",
    subtitle: "Start from physiology before modeling",
    slides: [
      {
        kind: "title",
        title: "Session 1\nHomeostasis and Systems Thinking",
        subtitle:
          "Borrow the organism-first teaching logic of BIEN5700 before we shift to cell-scale computational modeling.",
      },
      {
        kind: "figure",
        title: "The physiology course model we are borrowing from",
        imageKey: "bien5700_intro",
        takeaway:
          "BIEN5700 organizes learning by lesson, one concept at a time. That is the teaching stance we use here as well.",
      },
      {
        kind: "big-idea",
        title: "What is physiology in this course?",
        statement:
          "Physiology is the mechanistic study of how living systems function and stay organized.",
        supporting:
          "That definition matters because modeling is only useful when students can say what the system is, what it does, and what keeps it stable or unstable.",
      },
      {
        kind: "bullets",
        title: "What is homeostasis?",
        bullets: [
          "A living system maintains internal variables within useful ranges.",
          "Control is achieved through sensing, comparison, response, and feedback.",
          "Failure of control is often the starting point for disease and for intervention.",
        ],
      },
      {
        kind: "bullets",
        title: "Levels of biological organization",
        bullets: [
          "Organelle",
          "Cell",
          "Tissue",
          "Organ",
          "Whole organism",
        ],
      },
      {
        kind: "figure",
        title: "Feedback is the language of physiological control",
        imageKey: "bien5700_feedback",
        takeaway:
          "Students should be able to distinguish negative feedback, positive feedback, and situations where feedback alone is not enough to explain behavior.",
      },
      {
        kind: "big-idea",
        title: "Why this matters for electro-exocytosis",
        statement:
          "Electro-exocytosis is not just a pulse-response problem. It is a homeostasis-and-perturbation problem.",
        supporting:
          "Pulses disturb membranes, ions, organelles, and repair processes. The cell responds through feedback, recovery, adaptation, or failure.",
      },
    ],
  },
  {
    slug: "02_Transport_Compartments_and_Mass_Balance",
    label: "Session 2",
    title: "Transport, Compartments, and Mass Balance",
    subtitle: "The first real modeling language students need",
    slides: [
      {
        kind: "title",
        title: "Session 2\nTransport, Compartments, and Mass Balance",
        subtitle:
          "This session bridges physiology and mathematics by turning transport ideas into simple compartments and rates.",
      },
      {
        kind: "figure",
        title: "A physiology-first transport slide worth imitating",
        imageKey: "bien5700_transport",
        takeaway:
          "Even in a content-heavy course, BIEN5700 keeps this slide focused on one comparison only: physical versus physiological transport.",
      },
      {
        kind: "bullets",
        title: "Why transport comes before equations",
        bullets: [
          "Concentrations change because material moves, reacts, or is stored.",
          "If students understand transport, they can understand why state variables rise or fall.",
          "If they do not, differential equations feel arbitrary.",
        ],
      },
      {
        kind: "figure",
        title: "A compartment model is a well-mixed box with flows",
        imageKey: "bien6931_compartment_intro",
        takeaway:
          "The main simplifying assumption is that the compartment is well mixed over the timescale of interest.",
      },
      {
        kind: "big-idea",
        title: "The mass balance idea",
        statement: "Rate of change = inflow - outflow + production - consumption",
        supporting:
          "This one sentence is the backbone of compartmental modeling, PK/PD, calcium models, and later electroexo state equations.",
      },
      {
        kind: "bullets",
        title: "Why compartments are so useful",
        bullets: [
          "They keep the model interpretable.",
          "They let us connect variables to measurements.",
          "They are flexible enough for plasma, cytosol, ER, mitochondria, or EV pools.",
        ],
      },
      {
        kind: "big-idea",
        title: "Electroexo will reuse this pattern repeatedly",
        statement:
          "By the time students reach electroexo, they should already recognize compartments as a natural modeling choice.",
        supporting:
          "Calcium stores, mitochondrial states, EV pools, and damage states all inherit this logic.",
      },
    ],
  },
  {
    slug: "03_ODEs_and_Simulation_Basics",
    label: "Session 3",
    title: "ODEs and Simulation Basics",
    subtitle: "How state variables become time-evolving models",
    slides: [
      {
        kind: "title",
        title: "Session 3\nODEs and Simulation Basics",
        subtitle:
          "Students do not need abstract mathematics first. They need to know what a changing biological state looks like mathematically.",
      },
      {
        kind: "figure",
        title: "BIEN6931 begins ODE review with one core definition",
        imageKey: "bien6931_ode_intro",
        takeaway:
          "The teaching move to copy here is good: define a state variable, define its rate of change, and connect both to time.",
      },
      {
        kind: "bullets",
        title: "Words students must own",
        bullets: [
          "State variable",
          "Parameter",
          "Input",
          "Initial condition",
          "Output",
        ],
      },
      {
        kind: "big-idea",
        title: "Why ODEs appear so often in biology",
        statement:
          "Many biological questions ask how a well-mixed quantity changes over time.",
        supporting:
          "Concentrations, voltages, pool sizes, repair states, damage levels, and signaling activities are all natural ODE states.",
      },
      {
        kind: "bullets",
        title: "Linear and nonlinear behavior",
        bullets: [
          "Linear systems are easier to analyze and teach.",
          "Biological systems are often nonlinear because of thresholds, saturation, feedback, and coupling.",
          "Nonlinearity is not a detail. It is often the point.",
        ],
      },
      {
        kind: "bullets",
        title: "Why numerical simulation is necessary",
        bullets: [
          "Many realistic models do not have clean closed-form solutions.",
          "Simulation lets us study parameter changes, regimes, and sensitivity.",
          "It also forces us to be explicit about assumptions.",
        ],
      },
      {
        kind: "big-idea",
        title: "The key learning outcome of this session",
        statement:
          "Students should be able to read a simple ODE and explain the biological meaning of each term.",
        supporting:
          "That skill is more important here than advanced analytical solution techniques.",
      },
    ],
  },
  {
    slug: "04_PK_PD_and_PBPK_as_Modeling_Templates",
    label: "Session 4",
    title: "PK, PD, and PBPK as Modeling Templates",
    subtitle: "A reusable way to teach compartments, coupling, and interpretation",
    slides: [
      {
        kind: "title",
        title: "Session 4\nPK, PD, and PBPK as Modeling Templates",
        subtitle:
          "PK/PD is valuable here not because electroexo is a drug-delivery project, but because PK/PD teaches transferable modeling patterns.",
      },
      {
        kind: "figure",
        title: "Three phases between dose and effect",
        imageKey: "bien6931_pkpd_phases",
        takeaway:
          "This figure is excellent for teaching students how a single intervention can be split into linked stages with different meanings.",
      },
      {
        kind: "figure",
        title: "Pharmacokinetics is about concentration at the target site",
        imageKey: "bien6931_pk_vs_pd",
        takeaway:
          "Students should learn the simple slogan: PK is what the body does to the input; PD is what the input does to the system.",
      },
      {
        kind: "bullets",
        title: "Why PK/PD is still relevant to this course",
        bullets: [
          "It teaches exposure-response thinking.",
          "It shows how compartments can be linked in sequence.",
          "It separates upstream delivery from downstream effect.",
        ],
      },
      {
        kind: "figure",
        title: "PBPK adds anatomy and physiology to the compartments",
        imageKey: "bien6931_pbpk",
        takeaway:
          "PBPK demonstrates how a compartment model becomes more mechanistic when organ identity, flow, and physiology are made explicit.",
      },
      {
        kind: "big-idea",
        title: "The transferable lesson",
        statement:
          "Electroexo also separates an external input from internal transport, signaling, and effect layers.",
        supporting:
          "The vocabulary changes, but the modeling logic is familiar: perturbation, distribution, interaction, and outcome.",
      },
      {
        kind: "bullets",
        title: "What students should carry forward",
        bullets: [
          "Inputs can be staged.",
          "Compartments can be physiological.",
          "Outputs can be delayed relative to the intervention.",
        ],
      },
    ],
  },
  {
    slug: "05_Membranes_Gradients_and_Cellular_Electrophysiology",
    label: "Session 5",
    title: "Membranes, Gradients, and Cellular Electrophysiology",
    subtitle: "The bioelectric background needed before electroporation",
    slides: [
      {
        kind: "title",
        title: "Session 5\nMembranes, Gradients, and Cellular Electrophysiology",
        subtitle:
          "Students need a clean picture of ions, gradients, and membrane behavior before pulsed electric fields make sense.",
      },
      {
        kind: "figure",
        title: "Dynamic steady state depends on unequal ion distributions",
        imageKey: "bien5700_dynamic_state",
        takeaway:
          "This kind of slide is useful because it makes the membrane problem quantitative without overwhelming students with derivations yet.",
      },
      {
        kind: "bullets",
        title: "Three ideas to stabilize early",
        bullets: [
          "Cells maintain unequal ion concentrations across membranes.",
          "Membranes are selectively permeable, not perfectly open or perfectly closed.",
          "Voltage emerges from charge separation and permeability differences.",
        ],
      },
      {
        kind: "bullets",
        title: "Why electrophysiology belongs in this course",
        bullets: [
          "Pulsed electric fields interact with charged species and membranes.",
          "Membrane voltage is one of the first interpretable perturbation states after a pulse.",
          "Electrophysiology is therefore a prerequisite, not an optional side topic.",
        ],
      },
      {
        kind: "big-idea",
        title: "A useful teaching simplification",
        statement:
          "Students do not need the full Hodgkin-Huxley model before electroexo.",
        supporting:
          "They do need to understand what voltage means, why channels matter, and why short pulses can change membrane behavior.",
      },
      {
        kind: "figure",
        title: "Human-cell electrodynamics can be visualized spatially",
        imagePath: FIGURES.humanCellElectro,
        takeaway:
          "The 3D storyboard assets in the electroexo repository can help students map abstract electrical states onto cell structures.",
      },
      {
        kind: "big-idea",
        title: "Bridge to the next session",
        statement:
          "Once membranes and gradients are perturbed, calcium and organelle biology become unavoidable.",
        supporting:
          "That is why electrophysiology and calcium signaling must be taught as consecutive ideas.",
      },
    ],
  },
  {
    slug: "06_Calcium_Organelles_and_Bioenergetics",
    label: "Session 6",
    title: "Calcium, Organelles, and Bioenergetics",
    subtitle: "The biological heart of the middle layers",
    slides: [
      {
        kind: "title",
        title: "Session 6\nCalcium, Organelles, and Bioenergetics",
        subtitle:
          "This session prepares students for the states that dominate electroexo Layer 3 and the biology that feeds Layers 4 and 5.",
      },
      {
        kind: "bullets",
        title: "Why calcium deserves its own session",
        bullets: [
          "Calcium is a messenger.",
          "Calcium is a trigger.",
          "Calcium is also a stressor when control fails.",
        ],
      },
      {
        kind: "bullets",
        title: "The organelles students must remember",
        bullets: [
          "Endoplasmic reticulum: a major calcium store",
          "Mitochondria: ATP production and ROS coupling",
          "Endosomes and MVB-related compartments: later EV relevance",
        ],
      },
      {
        kind: "figure",
        title: "A spatial view of calcium influx and ER release",
        imagePath: FIGURES.humanCellCalcium,
        takeaway:
          "The point is not geometry for its own sake. The point is that calcium redistribution is tied to membranes and organelles, not just to a single bulk concentration.",
      },
      {
        kind: "figure",
        title: "Ion, mitochondrial, ROS, and ATP stress can also be visualized",
        imagePath: FIGURES.humanCellStress,
        takeaway:
          "This repository asset helps students see why calcium, mitochondrial state, ROS, and ATP are modeled together instead of as unrelated variables.",
      },
      {
        kind: "figure",
        title: "The current electroexo Layer 3 behavior is already interpretable",
        imagePath: FIGURES.layer3,
        takeaway:
          "Later, when students read the code, this figure will help them connect fluxes and state variables to recognizable biological trajectories.",
      },
      {
        kind: "big-idea",
        title: "Key takeaway",
        statement:
          "A brief electrical perturbation can create longer-lived calcium, metabolic, and stress trajectories.",
        supporting:
          "That time-scale expansion is one of the main reasons a layered computational framework is useful.",
      },
    ],
  },
  {
    slug: "07_Exocytosis_Membrane_Repair_and_Extracellular_Vesicles",
    label: "Session 7",
    title: "Exocytosis, Membrane Repair, and Extracellular Vesicles",
    subtitle: "The release biology students need before EV engineering",
    slides: [
      {
        kind: "title",
        title: "Session 7\nExocytosis, Membrane Repair, and Extracellular Vesicles",
        subtitle:
          "Students should understand the biology of release and repair before they encounter any electroexo-specific abstractions.",
      },
      {
        kind: "bullets",
        title: "Exocytosis in one sentence",
        bullets: [
          "Membrane-bound cargo is delivered to the cell boundary and released through regulated fusion processes.",
        ],
      },
      {
        kind: "bullets",
        title: "Membrane repair in one sentence",
        bullets: [
          "Cells actively detect membrane injury and recruit calcium-dependent processes to reseal, remodel, or shed damaged regions.",
        ],
      },
      {
        kind: "bullets",
        title: "Extracellular vesicles are not one thing",
        bullets: [
          "Small EVs",
          "Medium/large EVs",
          "Apoptotic bodies",
        ],
      },
      {
        kind: "figure",
        title: "The manuscript framework makes the release routes visible",
        imagePath: FIGURES.framework,
        takeaway:
          "This is one of the best visuals to teach that vesicle output can arise from multiple biological routes with different meanings.",
      },
      {
        kind: "figure",
        title: "Repair and release are linked in the current framework",
        imagePath: FIGURES.remodeling,
        takeaway:
          "Students should see that repair is not separate from release logic. Repair states help shape what kind of vesicles appear and when.",
      },
      {
        kind: "big-idea",
        title: "Key takeaway",
        statement:
          "Electro-exocytosis is meaningful only if we can separate productive release from injury-associated debris.",
        supporting:
          "That distinction becomes a major design principle for the later layers of the framework.",
      },
    ],
  },
  {
    slug: "08_nsPEF_Electroporation_and_Dose",
    label: "Session 8",
    title: "nsPEF, Electroporation, and Dose",
    subtitle: "Now the external perturbation can be introduced",
    slides: [
      {
        kind: "title",
        title: "Session 8\nnsPEF, Electroporation, and Dose",
        subtitle:
          "Only after the students have physiology, transport, ODEs, and release biology do we introduce the pulse itself.",
      },
      {
        kind: "bullets",
        title: "What is an nsPEF?",
        bullets: [
          "A nanosecond pulsed electric field is an ultrashort, high-intensity electrical perturbation.",
          "Its duration matters, not only its amplitude.",
          "Its biological effects can extend beyond the plasma membrane.",
        ],
      },
      {
        kind: "figure",
        title: "A cell-level picture of pulse exposure",
        imagePath: FIGURES.humanCellPulse,
        takeaway:
          "Students often need a concrete picture first. This frame helps them visualize the event before equations and dose metrics are introduced.",
      },
      {
        kind: "figure",
        title: "The dosimetry schematic is useful because it separates inputs from interpretations",
        imagePath: FIGURES.dosimetrySchematic,
        takeaway:
          "Pulse width, pulse number, conductivity, geometry, and retained heating do not mean the same thing, and the slide should make that explicit.",
      },
      {
        kind: "figure",
        title: "Heating and charging are different ideas",
        imagePath: FIGURES.temperatureProfiles,
        takeaway:
          "A dose metric is not automatically a temperature metric. That distinction matters for teaching and for model structure.",
      },
      {
        kind: "figure",
        title: "Pulse width changes membrane charging behavior",
        imagePath: FIGURES.membraneVoltage,
        takeaway:
          "This is the cleanest bridge from pulse protocol to the electrodynamic states used later in the codebase.",
      },
      {
        kind: "big-idea",
        title: "Key takeaway",
        statement:
          "The pulse is an external forcing input, but its biological meaning is mediated through membrane and organelle perturbation states.",
        supporting:
          "That idea prepares students for Layers 1 and 2 of the framework.",
      },
    ],
  },
  {
    slug: "09_The_Electroexo_Framework",
    label: "Session 9",
    title: "The Electroexo Framework",
    subtitle: "Now the background is finally assembled into the project model",
    slides: [
      {
        kind: "title",
        title: "Session 9\nThe Electroexo Framework",
        subtitle:
          "At this point students have the background to understand why the framework is layered the way it is.",
      },
      {
        kind: "big-idea",
        title: "Why a layered framework is needed",
        statement:
          "No single state variable can connect an nsPEF protocol to EV quality, yield, subtype balance, and viability.",
        supporting:
          "The model needs multiple linked abstractions because the biology unfolds across different mechanisms and timescales.",
      },
      {
        kind: "bullets",
        title: "Layers 1 and 2",
        bullets: [
          "Pulse descriptors and dose interpretation",
          "Membrane and organelle electrodynamics",
          "The immediate perturbation state",
        ],
      },
      {
        kind: "bullets",
        title: "Layer 3",
        bullets: [
          "Calcium mobilization",
          "Mitochondrial potential and ATP",
          "ROS and osmotic stress",
        ],
      },
      {
        kind: "bullets",
        title: "Layers 4 and 5",
        bullets: [
          "Repair and remodeling",
          "MVB, budding, and release states",
          "Subtype-resolved EV output",
        ],
      },
      {
        kind: "bullets",
        title: "Layers 6, 7, 8, and modifiers",
        bullets: [
          "Cargo and potency",
          "Injury and quality gate",
          "Manufacturing and QC",
          "Cell-state modifiers across all layers",
        ],
      },
      {
        kind: "figure",
        title: "The storyboard is the best single summary of the integrated logic",
        imagePath: FIGURES.storyboard,
        takeaway:
          "This figure helps students see how mild, productive, and injurious regimes can arise from the same framework with different upstream conditions.",
      },
    ],
  },
  {
    slug: "10_Evidence_Calibration_and_Project_Directions",
    label: "Session 10",
    title: "Evidence, Calibration, and Project Directions",
    subtitle: "Teach students how to extend the framework honestly",
    slides: [
      {
        kind: "title",
        title: "Session 10\nEvidence, Calibration, and Project Directions",
        subtitle:
          "The last session teaches students how a computational framework becomes better, not just bigger.",
      },
      {
        kind: "big-idea",
        title: "The framework is useful because it is explicit about its current limits",
        statement:
          "A placeholder model can still be educational and scientifically useful if its assumptions are visible.",
        supporting:
          "That transparency is one of the strengths of the electroexo repository and manuscript pairing.",
      },
      {
        kind: "bullets",
        title: "What calibration means in practice",
        bullets: [
          "Choose observables that map onto actual model terms.",
          "Constrain upstream layers before fitting downstream outputs.",
          "Do not confuse literature support with parameter-identifying evidence.",
        ],
      },
      {
        kind: "table",
        title: "A staged calibration order for students to remember",
        headers: ["Stage", "What to constrain first", "Why"],
        rows: [
          ["1", "Dose and electrodynamics", "Stabilizes the exposure map"],
          ["2", "Calcium, ROS, ATP, stress", "Constrains the main signaling layer"],
          ["3", "Repair and release observables", "Separates recovery from injury"],
          ["4", "Subtype output and cargo", "Adds engineering interpretation"],
        ],
        note: "This ordering follows the model logic instead of chasing whichever dataset is easiest to fit.",
      },
      {
        kind: "bullets",
        title: "Project expectations adapted from BIEN6931",
        bullets: [
          "State the problem clearly.",
          "List assumptions and governing equations clearly.",
          "Show simulations or figures clearly.",
          "Discuss limitations clearly.",
        ],
      },
      {
        kind: "bullets",
        title: "Strong final-project directions for this course",
        bullets: [
          "Improve one submodel.",
          "Build one evidence table.",
          "Propose one validation plan.",
          "Explain one major limitation.",
        ],
      },
      {
        kind: "big-idea",
        title: "Final course takeaway",
        statement:
          "A good computational systems biology course does not only teach students to run a model. It teaches them how to think about model structure, evidence, and biological meaning.",
        supporting:
          "That is the mindset this revised lecture sequence is built to support.",
      },
    ],
  },
];

const FIGURE_LIST = [
  {
    category: "BIEN5700 physiology teaching figures",
    items: [
      ["BIEN5700 page 6", "Lesson 1 homeostasis title slide", "Strong model for simple lesson openers and concept-first pacing."],
      ["BIEN5700 page 11", "Feedback control comparison", "Useful for introducing homeostasis and control logic."],
      ["BIEN5700 page 16", "Physical vs physiological transport", "Useful for transport prerequisites and membrane-transport classification."],
      ["BIEN5700 page 31", "Dynamic steady-state ion distribution table", "Useful for introducing gradients, membrane potential, and bioelectricity."],
    ],
  },
  {
    category: "BIEN6931 modeling-prerequisite figures",
    items: [
      ["Dash_Lect1 page 2", "What is systems biology?", "Useful for defining non-intuitive, nonlinear biological systems."],
      ["Dash_Lect1 page 3", "What is computational systems biology?", "Useful for framing modeling, simulation, and hypothesis generation."],
      ["Compartment modeling handout page 1", "Well-mixed compartment diagram", "Useful for introducing compartment assumptions and mass balance."],
      ["PKPD page 2", "Three phases from dose to effect", "Useful for staged input-response thinking."],
      ["PKPD page 3", "PK versus PD flow", "Useful for separating transport/exposure from effect."],
      ["PBPK page 2", "Whole-body PBPK schematic", "Useful for physiology-informed compartment modeling."],
    ],
  },
  {
    category: "electroexo and manuscript figures",
    items: [
      ["fig01.png", "Integrated electro-exocytosis mechanism", "Best high-level overview of pulse-to-EV pathways."],
      ["figA1_dosimetry_schematic.png", "Dosimetry and pulse interpretation", "Best figure for teaching Layer 1 concepts."],
      ["temperature_profiles.png", "Retained heating profiles", "Useful for separating dose from heating."],
      ["membrane_voltage_by_pulse_width.png", "Voltage vs pulse width", "Useful for Layer 2 and charging intuition."],
      ["layer3_timeseries_comparison.png", "Calcium/ROS/ATP trajectories", "Useful for Layer 3 dynamics and regime comparison."],
      ["remodeling_repair_timeseries.png", "Repair-state trajectories", "Useful for Layer 4 teaching."],
      ["ev_biogenesis_timeseries.png", "EV pool and subtype trajectories", "Useful for Layer 5 teaching."],
      ["multilayer_storyboard.png", "Integrated mild/productive/injurious regimes", "Best single summary figure for the whole framework."],
    ],
  },
  {
    category: "electroexo visualization assets",
    items: [
      ["results/human_cell_3d/storyboard_frames/01_01_nsPEF_exposure.png", "3D pulse exposure frame", "Useful for introducing the perturbation visually."],
      ["results/human_cell_3d/storyboard_frames/02_02_membrane_electrodynamics.png", "3D membrane electrodynamics frame", "Useful for Layer 2 visualization."],
      ["results/human_cell_3d/storyboard_frames/03_03_Ca_influx_ER_release.png", "3D calcium and ER frame", "Useful for calcium-redistribution teaching."],
      ["results/human_cell_3d/storyboard_frames/04_04_ionic_mito_ROS_ATP.png", "3D mitochondrial stress frame", "Useful for bioenergetics and stress coupling."],
      ["results/human_cell_3d/storyboard_frames/05_05_remodeling_repair.png", "3D repair frame", "Useful for membrane repair discussions."],
      ["results/human_cell_3d/storyboard_frames/06_06_EV_release.png", "3D EV release frame", "Useful for final pathway and EV-release visualization."],
      ["results/ev_release_regime_map/ev_release_regime_map.png", "EV regime heatmap", "Useful for parameter-sweep and design-space discussions."],
    ],
  },
];

async function ensureDir(dir) {
  await fs.mkdir(dir, { recursive: true });
}

function addTextBox(slide, { text, left, top, width, height, fontSize = 18, bold = false, color = COLORS.text, alignment = "left" }) {
  const box = slide.shapes.add({
    geometry: "textbox",
    position: { left, top, width, height },
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  box.text = text;
  box.text.style = { fontSize, bold, color, alignment };
  return box;
}

function addHeader(slide, deck, title, index, count) {
  slide.background.fill = COLORS.bg;
  slide.shapes.add({
    geometry: "rect",
    position: { left: 0, top: 0, width: SLIDE_SIZE.width, height: 76 },
    fill: COLORS.teal,
    line: { style: "solid", fill: COLORS.teal, width: 0 },
  });
  addTextBox(slide, {
    text: title,
    left: 74,
    top: 12,
    width: 1020,
    height: 46,
    fontSize: 40,
    bold: true,
    color: "#FFFFFF",
  });
  addTextBox(slide, {
    text: deck.label,
    left: 1080,
    top: 18,
    width: 120,
    height: 24,
    fontSize: 15,
    color: "#E5F5F8",
    alignment: "right",
  });
  addTextBox(slide, {
    text: `Slide ${index + 1} of ${count}`,
    left: 1040,
    top: 682,
    width: 160,
    height: 18,
    fontSize: 12,
    color: COLORS.faint,
    alignment: "right",
  });
}

function addBulletList(slide, bullets, opts) {
  const joined = bullets.map((bullet) => `• ${bullet}`).join("\n\n");
  addTextBox(slide, { ...opts, text: joined });
}

async function addImage(slide, imagePath, position) {
  const bytes = await fs.readFile(imagePath);
  slide.images.add({
    blob: bytes,
    contentType: "image/png",
    alt: path.basename(imagePath),
    fit: "contain",
    position,
  });
}

function addSoftBand(slide, top, text) {
  slide.shapes.add({
    geometry: "rect",
    position: { left: 74, top, width: 1132, height: 62 },
    fill: COLORS.soft,
    line: { style: "solid", fill: COLORS.soft, width: 0 },
  });
  addTextBox(slide, {
    text,
    left: 98,
    top: top + 16,
    width: 1080,
    height: 28,
    fontSize: 22,
    color: COLORS.tealDark,
  });
}

function styleTable(table, rows, columns) {
  table.styleOptions = { headerRow: true, bandedRows: true };
  table.borders.assign({ style: "solid", fill: "#D1D5DB", width: 1 });
  for (let c = 0; c < columns; c += 1) {
    const cell = table.getCell(0, c);
    cell.fill = COLORS.teal;
    cell.text.style = { fontSize: 16, bold: true, color: "#FFFFFF" };
  }
  for (let r = 1; r < rows; r += 1) {
    for (let c = 0; c < columns; c += 1) {
      table.getCell(r, c).text.style = { fontSize: 15, color: COLORS.text };
    }
  }
}

function extractPdfFromZip(zipPath, memberName, outPath) {
  const data = execFileSync("unzip", ["-p", zipPath, memberName], {
    encoding: "buffer",
    maxBuffer: 32 * 1024 * 1024,
  });
  return fs.writeFile(outPath, data);
}

function renderPdfPage(pdfPath, pageNumber, outPrefix) {
  execFileSync(
    "/Users/aghorban/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/pdftoppm",
    ["-f", String(pageNumber), "-l", String(pageNumber), "-singlefile", "-png", pdfPath, outPrefix],
    {
      env: {
        ...process.env,
        XDG_CACHE_HOME: "/private/tmp/fontcache",
        FONTCONFIG_PATH: "/etc/fonts",
      },
    },
  );
  return `${outPrefix}.png`;
}

async function prepareAssets(assetDir) {
  await ensureDir(assetDir);
  await ensureDir("/private/tmp/fontcache");
  const dashIntroPdf = path.join(assetDir, "Dash_Lect1_Intro_Sim_Biol_Syst_Spring_2018.pdf");
  const pkpdPdf = path.join(assetDir, "Compartmental_Kinetic_PKPD_Modeling.pdf");
  const pbpkPdf = path.join(assetDir, "PBPK_spring_2018.pdf");
  const compartmentPdf = path.join(assetDir, "Introduction_to_Compartmental_Modeling_1.pdf");

  await extractPdfFromZip(ZIP_6931, "Dash_Lect1_Intro_Sim_Biol_Syst_Spring_2018.pdf", dashIntroPdf);
  await extractPdfFromZip(ZIP_6931, "Compartmental_Kinetic_PKPD_Modeling.pdf", pkpdPdf);
  await extractPdfFromZip(ZIP_6931, "PBPK_spring_2018.pdf", pbpkPdf);
  await extractPdfFromZip(ZIP_6931, "Introduction to Compartmental Modeling (1).pdf", compartmentPdf);

  return {
    bien5700_intro: renderPdfPage(PDF_5700_HOMEO, 6, path.join(assetDir, "bien5700_intro")),
    bien5700_feedback: renderPdfPage(PDF_5700_HOMEO, 11, path.join(assetDir, "bien5700_feedback")),
    bien5700_transport: renderPdfPage(PDF_5700_HOMEO, 16, path.join(assetDir, "bien5700_transport")),
    bien5700_dynamic_state: renderPdfPage(PDF_5700_HOMEO, 31, path.join(assetDir, "bien5700_dynamic_state")),
    bien6931_systems_biology: renderPdfPage(dashIntroPdf, 2, path.join(assetDir, "bien6931_systems_biology")),
    bien6931_comp_systems_biology: renderPdfPage(dashIntroPdf, 3, path.join(assetDir, "bien6931_comp_systems_biology")),
    bien6931_ode_intro: renderPdfPage(
      path.join(assetDir, "Differential_equations_Spring_2018.pdf"),
      1,
      path.join(assetDir, "bien6931_ode_intro"),
    ),
    bien6931_compartment_intro: renderPdfPage(compartmentPdf, 1, path.join(assetDir, "bien6931_compartment_intro")),
    bien6931_pkpd_phases: renderPdfPage(pkpdPdf, 2, path.join(assetDir, "bien6931_pkpd_phases")),
    bien6931_pk_vs_pd: renderPdfPage(pkpdPdf, 3, path.join(assetDir, "bien6931_pk_vs_pd")),
    bien6931_pbpk: renderPdfPage(pbpkPdf, 2, path.join(assetDir, "bien6931_pbpk")),
  };
}

async function ensureDifferentialEquationsPdf(assetDir) {
  const outPath = path.join(assetDir, "Differential_equations_Spring_2018.pdf");
  await extractPdfFromZip(ZIP_6931, "Differential equations_Spring_2018.pdf", outPath);
  return outPath;
}

async function renderSlide(slide, deck, spec, idx, count, assetMap) {
  if (spec.kind === "title") {
    slide.background.fill = COLORS.bg;
    slide.shapes.add({
      geometry: "rect",
      position: { left: 0, top: 0, width: SLIDE_SIZE.width, height: 720 },
      fill: COLORS.bg,
      line: { style: "solid", fill: COLORS.bg, width: 0 },
    });
    slide.shapes.add({
      geometry: "rect",
      position: { left: 0, top: 0, width: SLIDE_SIZE.width, height: 96 },
      fill: COLORS.teal,
      line: { style: "solid", fill: COLORS.teal, width: 0 },
    });
    addTextBox(slide, {
      text: deck.label,
      left: 78,
      top: 36,
      width: 180,
      height: 22,
      fontSize: 18,
      bold: true,
      color: "#E5F5F8",
    });
    addTextBox(slide, {
      text: spec.title,
      left: 98,
      top: 148,
      width: 1080,
      height: 220,
      fontSize: 52,
      bold: true,
      color: COLORS.text,
    });
    addTextBox(slide, {
      text: spec.subtitle,
      left: 100,
      top: 406,
      width: 980,
      height: 140,
      fontSize: 25,
      color: COLORS.muted,
    });
    addTextBox(slide, {
      text: `${idx + 1} / ${count}`,
      left: 1120,
      top: 678,
      width: 80,
      height: 20,
      fontSize: 12,
      color: COLORS.faint,
      alignment: "right",
    });
    return;
  }

  addHeader(slide, deck, spec.title, idx, count);

  if (spec.kind === "big-idea") {
    addTextBox(slide, {
      text: spec.statement,
      left: FRAME.left,
      top: 176,
      width: FRAME.width,
      height: 170,
      fontSize: 36,
      bold: true,
      color: COLORS.text,
    });
    addSoftBand(slide, 416, spec.supporting);
    return;
  }

  if (spec.kind === "bullets") {
    addBulletList(slide, spec.bullets, {
      left: FRAME.left,
      top: 176,
      width: FRAME.width - 40,
      height: 430,
      fontSize: 28,
      color: COLORS.text,
    });
    return;
  }

  if (spec.kind === "figure") {
    const imagePath = spec.imagePath ?? assetMap[spec.imageKey];
    await addImage(slide, imagePath, { left: 120, top: 154, width: 1040, height: 420 });
    addSoftBand(slide, 610, spec.takeaway);
    return;
  }

  if (spec.kind === "table") {
    const values = [spec.headers, ...spec.rows];
    const table = slide.tables.add({
      rows: values.length,
      columns: spec.headers.length,
      left: 74,
      top: 160,
      width: 1132,
      height: 410,
      values,
    });
    styleTable(table, values.length, spec.headers.length);
    addSoftBand(slide, 602, spec.note);
    return;
  }

  throw new Error(`Unknown slide kind: ${spec.kind}`);
}

async function exportDeck(deck, finalDir, previewDir, assetMap) {
  const presentation = Presentation.create({ slideSize: SLIDE_SIZE });
  for (let i = 0; i < deck.slides.length; i += 1) {
    const slide = presentation.slides.add();
    await renderSlide(slide, deck, deck.slides[i], i, deck.slides.length, assetMap);
  }

  const outPath = path.join(finalDir, `${deck.slug}.pptx`);
  const file = await PresentationFile.exportPptx(presentation);
  await file.save(outPath);

  const montage = await presentation.export({ format: "webp", montage: true, scale: 1 });
  await fs.writeFile(path.join(previewDir, `${deck.slug}-montage.webp`), new Uint8Array(await montage.arrayBuffer()));
  return outPath;
}

async function writeFigureList(outDir) {
  const lines = ["Figure Adaptation List", "======================", ""];
  for (const section of FIGURE_LIST) {
    lines.push(section.category);
    lines.push("-".repeat(section.category.length));
    for (const [source, figure, note] of section.items) {
      lines.push(`* Source: ${source}`);
      lines.push(`  Figure/topic: ${figure}`);
      lines.push(`  Why useful: ${note}`);
      lines.push("");
    }
  }
  await fs.writeFile(path.join(outDir, "figure_adaptation_list.txt"), `${lines.join("\n")}\n`, "utf8");
}

async function writeCourseNotes(outDir) {
  const notes = [
    "General-Audience Course Notes",
    "============================",
    "",
    "Design principles used in this revision:",
    "1. The course begins from physiology and modeling prerequisites, not from the manuscript structure.",
    "2. Each slide is built around one main idea, one comparison, or one figure.",
    "3. BIEN5700 supplied the organism-first, lesson-by-lesson pacing.",
    "4. BIEN6931 supplied the computational prerequisite ladder and project expectations.",
    "5. The electroexo framework is introduced only after students have enough background to read it.",
    "",
    "Core source anchors used:",
    "- BIEN5700 Systems Physiology syllabus and lecture PDFs",
    "- BIEN6931 Modeling and Simulations of Integrated Cellular Systems syllabus and lecture PDFs",
    "- electroexo repository outputs and 3D storyboard frames",
    "- electro-exocytosis manuscript figures",
  ];
  await fs.writeFile(path.join(outDir, "course_design_notes.txt"), `${notes.join("\n")}\n`, "utf8");
}

async function main() {
  const finalDir = process.argv[2];
  const previewDir = process.argv[3];
  const assetDir = process.argv[4];
  if (!finalDir || !previewDir || !assetDir) {
    throw new Error("Usage: node generate_electroexo_course_general_audience.mjs <finalDir> <previewDir> <assetDir>");
  }

  await ensureDir(finalDir);
  await ensureDir(previewDir);
  await ensureDir(assetDir);
  await ensureDifferentialEquationsPdf(assetDir);
  const assetMap = await prepareAssets(assetDir);

  const outPaths = [];
  for (const deck of DECKS) {
    outPaths.push(await exportDeck(deck, finalDir, previewDir, assetMap));
  }
  await writeFigureList(finalDir);
  await writeCourseNotes(finalDir);
  console.log(`Generated ${outPaths.length} decks in ${finalDir}`);
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
