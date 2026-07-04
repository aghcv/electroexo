import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const ODU_BLUE = "#043657";
const ODU_RED = "#C10E20";
const TEXT_DARK = "#2E2C2C";
const TEXT_MUTED = "#636466";

const starterPptx = process.argv[2];
const finalPptx = process.argv[3];
const finalDir = process.argv[4] ?? path.dirname(finalPptx);
const previewDir = process.argv[5] ?? path.join(finalDir, "preview");

if (!starterPptx || !finalPptx) {
  console.error(
    "Usage: node generate_electroexo_syllabus_ag_template.mjs <starter.pptx> <final.pptx> [finalDir] [previewDir]",
  );
  process.exit(1);
}

const footerText = "Computational Systems Biology: Electro-Exocytosis";

const contentSlides = [
  {
    title: "Course Goal",
    kicker: "A cell-biology course that becomes a modeling course.",
    body: [
      "Students learn how nsPEF stimulation perturbs membranes, calcium, organelles, repair, trafficking, and EV release.",
      "Students then translate those mechanisms into explicit systems-biology models.",
      "The course ends with a minimal electro-exocytosis model that is interpretable, testable, and extensible.",
    ],
    note:
      "This opening keeps the course anchored in student preparation, not in the manuscript structure.",
  },
  {
    title: "Audience And Prerequisites",
    kicker: "Mixed backgrounds are expected.",
    body: [
      "Designed for senior undergraduates and beginning graduate students in biology, engineering, or computational science.",
      "Helpful background: introductory cell biology or physiology, calculus, basic ODEs, and basic programming.",
      "The first half builds the shared biological and modeling vocabulary before the electroexo framework appears.",
    ],
    note:
      "The prerequisite bridge follows the spirit of BIEN6931, but lowers the entry barrier for mixed-background students.",
  },
  {
    title: "Learning Objectives",
    kicker: "By the end, students should be able to do six things.",
    body: [
      "Explain membrane structure, transport, organelles, vesicle trafficking, exocytosis, endocytosis, and EV biogenesis.",
      "Describe electroporation, nanoporation, nsPEF dose variables, and major cellular response pathways.",
      "Build compartmental and ODE models for calcium, membrane repair, vesicle pools, and EV release.",
      "Use Python/SciPy and at least one systems-biology platform such as COPASI, VCell, SBML, or BioModels.",
      "Evaluate assumptions, parameters, sensitivity, uncertainty, and model failure modes.",
      "Present a minimal model of nsPEF-induced electro-exocytosis with clear biological interpretation.",
    ],
    note:
      "These objectives deliberately combine biology, biophysics, computation, and communication.",
  },
  {
    title: "Course Architecture",
    kicker: "The syllabus moves from foundations to perturbation to synthesis.",
    body: [
      "Part 1: Cell foundations: membranes, transport, organelles, vesicle traffic, repair, and EV biology.",
      "Part 2: Modeling foundations: compartments, mass balance, ODEs, nonlinear feedback, sensitivity, and simulation.",
      "Part 3: Bioelectric perturbation: electroporation, nsPEF dose, calcium, ROS, mitochondria, and remodeling.",
      "Part 4: Integration: minimal electro-exocytosis model, calibration, uncertainty, and final projects.",
    ],
    note:
      "This is the main organizing logic for the course. It is intentionally not a chapter-by-chapter report outline.",
  },
  {
    title: "Resource Strategy",
    kicker: "Use teaching resources first, then research papers.",
    body: [
      "Core biology: Molecular Biology of the Cell, OpenStax Biology 2e, NCBI Bookshelf, HHMI BioInteractive, EMBL-EBI Training.",
      "Computational biology: MIT OCW systems biology, SBML, BioModels, COPASI, VCell, Python/SciPy, Jupyter notebooks.",
      "Specialized literature: electroporation/nsPEF reviews, calcium signaling papers, membrane-repair papers, and EV standards such as MISEV2023.",
      "Project resources: the electroexo repository, existing model outputs, storyboard visualizations, and the current bibliography.",
    ],
    note:
      "This slide names categories only. The text syllabus contains the first-pass source list; the next phase can assign specific readings by week.",
  },
  {
    title: "Weekly Modules: Foundations",
    kicker: "Weeks 1-7 build a common language.",
    body: [
      "Week 1: Course map: from cells to models to electro-exocytosis.",
      "Week 2: Membrane structure, permeability, transport, and electrical gradients.",
      "Week 3: Organelles, calcium stores, mitochondria, ROS, and ATP.",
      "Week 4: Vesicle trafficking, endocytosis, exocytosis, SNAREs, Rab proteins, and MVBs.",
      "Week 5: Plasma membrane repair, lysosomal exocytosis, annexins, ESCRT, and EV biogenesis.",
      "Week 6: Compartmental modeling, mass balance, ODEs, rate laws, and numerical simulation.",
      "Week 7: Nonlinearity, feedback, regimes, sensitivity analysis, and model credibility.",
    ],
    note:
      "The first seven weeks are the on-ramp. Students should not need to know the project before these topics.",
  },
  {
    title: "Weekly Modules: Electroexo",
    kicker: "Weeks 8-14 connect nsPEF biology to models.",
    body: [
      "Week 8: Bioelectricity, membrane capacitance, transmembrane voltage, and field exposure.",
      "Week 9: Electroporation, nanoporation, nsPEF dosimetry, heating, and membrane charging.",
      "Week 10: nsPEF responses: calcium entry, ER release, mitochondria, ROS, ATP, and cytoskeleton.",
      "Week 11: Coupling repair, remodeling, vesicle trafficking, and EV release.",
      "Week 12: Minimal electro-exocytosis model: inputs, states, parameters, and outputs.",
      "Week 13: Calibration, literature evidence, uncertainty, reproducibility, and model comparison.",
      "Week 14: Final model presentations and peer critique.",
    ],
    note:
      "Week 12 is the first point where the full framework becomes central.",
  },
  {
    title: "Computational Labs",
    kicker: "Labs translate each biological mechanism into a small model.",
    body: [
      "Lab 1: Python warm-up: one-compartment mass balance and exponential recovery.",
      "Lab 2: Membrane transport: ion gradients, permeability, and simple flux models.",
      "Lab 3: Calcium pulse model: influx, buffering, ER release, and clearance.",
      "Lab 4: nsPEF membrane charging and pore-state toy model.",
      "Lab 5: Repair and vesicle-pool model with calcium-dependent activation.",
      "Lab 6: EV release regime map and sensitivity analysis.",
      "Lab 7: Export or compare a small model using SBML, COPASI, VCell, or BioModels.",
    ],
    note:
      "Each lab should be small enough to finish in class or in a short assignment, then become a project building block.",
  },
  {
    title: "Final Project",
    kicker: "Build a minimal systems-biology model of nsPEF-induced electro-exocytosis.",
    body: [
      "Define a biological question, such as productive EV release versus injury after nsPEF exposure.",
      "Choose a minimal state set: membrane perturbation, cytosolic Ca2+, organelle stress, repair state, vesicle pool, and EV output.",
      "Implement the model in Python/SciPy, COPASI, VCell, or SBML-compatible tooling.",
      "Run at least two simulations: a baseline condition and one pulse/dose perturbation.",
      "Report assumptions, parameter sources, sensitivity, limitations, and biological interpretation.",
    ],
    note:
      "The final project mirrors the MIT OCW-style emphasis on biological problem, computational problem, assumptions, and controls.",
  },
  {
    title: "Assessment Plan",
    kicker: "Assessment should reward mechanistic thinking and honest modeling.",
    body: [
      "Short reading responses: identify the mechanism, modelable variables, and missing measurements.",
      "Computational labs: submit notebook, equations, plots, and interpretation.",
      "Model critique: evaluate another model's assumptions, parameter choices, and failure modes.",
      "Final project: specific aims, model diagram, reproducible implementation, results, and presentation.",
      "Suggested weighting: labs 35%, reading/model critiques 20%, final project 35%, participation 10%.",
    ],
    note:
      "This keeps evaluation focused on transparent reasoning rather than whether a model gives an impressive-looking result.",
  },
  {
    title: "Next Resource Pass",
    kicker: "This draft is ready for module-level reading selection.",
    body: [
      "For each week, choose one teaching anchor, one short review, and one optional primary paper.",
      "Prioritize resources with clear figures, accessible language, and mechanisms that can become equations.",
      "Build lecture visuals around one idea per slide: process, assumption, state variable, equation, or result.",
      "Keep a resource matrix: module, concept, source, figure candidate, lab connection, and difficulty level.",
    ],
    note:
      "This is the handoff point for the collaborative resource curation phase requested by the user.",
  },
];

const detailedSyllabus = `Draft Syllabus
==============

Course title
Computational Systems Biology of Electro-Exocytosis

Course goal
This specialized course prepares students to understand how electrical stimulation, especially nanosecond pulsed electric fields (nsPEF), perturbs cellular membranes, calcium signaling, organelles, vesicle trafficking, plasma membrane repair, and extracellular vesicle (EV) release. Students learn the cell biology first, then learn how to represent those processes with compartmental, ODE-based, and systems-biology models.

Audience
Senior undergraduates and beginning graduate students with mixed backgrounds in biology, biomedical engineering, electrical engineering, computational science, or applied mathematics.

Recommended prerequisites
- Introductory biology, cell biology, physiology, or biomedical engineering.
- Calculus and basic ordinary differential equations.
- Basic programming experience, preferably Python.
- No prior nsPEF, extracellular vesicle, or advanced electrophysiology experience is assumed.

Learning objectives
1. Explain the cell biology foundations needed for electro-exocytosis: membrane structure, membrane transport, organelles, calcium signaling, vesicle trafficking, exocytosis, endocytosis, plasma membrane repair, and EV biogenesis.
2. Explain bioelectric stimulation concepts: membrane potential, capacitance, transmembrane voltage, electroporation, nanoporation, pulse duration, pulse number, pulse repetition rate, absorbed energy, heating, and nsPEF-induced cellular responses.
3. Translate biological mechanisms into systems-biology model components: states, compartments, fluxes, rate laws, parameters, feedback loops, and outputs.
4. Build, simulate, and interpret small ODE models in Python/SciPy and compare them with systems-biology tooling such as COPASI, VCell, SBML, and BioModels.
5. Evaluate model assumptions, parameter uncertainty, sensitivity, regime behavior, and limitations.
6. Build and present a minimal model of nsPEF-induced electro-exocytosis.

Recommended resource categories

Core cell biology and physiology
- Molecular Biology of the Cell by Alberts et al.; use current textbook access if available and NCBI Bookshelf chapters as accessible background.
- OpenStax Biology 2e for accessible chapters on cell structure, membranes, energy, and signaling foundations.
- HHMI BioInteractive for short visual teaching assets.
- EMBL-EBI Training for molecular and cellular biology teaching modules.
- FEBS review articles as candidate short reviews for advanced undergraduate reading.
- Local BIEN5700 materials for physiology-first pacing and homeostasis/transport framing.

Computational and systems biology
- MIT OpenCourseWare 7.91J Foundations of Computational and Systems Biology for project expectations, readings, and computational-biology framing.
- SBML for model exchange.
- BioModels for curated model examples.
- COPASI for biochemical network and ODE simulation.
- VCell for spatial and compartmental cellular modeling.
- Python/SciPy/Jupyter for transparent, inspectable modeling workflows.
- CellOrganizer as a possible resource for spatial cell organization; verify access during the next resource pass.
- Local BIEN6931 materials for computational systems-biology expectations, compartment modeling, PK/PD, PBPK, ODEs, and project structure.

Specialized nsPEF, electroporation, repair, calcium, and EV literature
- Kotnik et al. 2019, Membrane electroporation and electropermeabilization: mechanisms and models.
- Krassowska and Filev 2007, Modeling electroporation in a single cell.
- Smith and Weaver 2004, A model of electroporation-induced cellular uptake and membrane conductance changes.
- Beebe et al. 2003, Nanosecond high-intensity pulsed electric fields induce apoptosis in human cells.
- Vernier et al. 2003, Calcium bursts induced by nanosecond electric pulses.
- Orlacchio et al. 2023, nsPEF effects in multicellular spheroid tumor models with pulse duration, repetition rate, absorbed energy, and temperature.
- Reddy et al. 2001, Plasma membrane repair mediated by Ca2+-regulated lysosomal exocytosis.
- McNeil and Steinhardt 2003, Plasma membrane disruption: repair, prevention, adaptation.
- Andrews et al. 2014, Damage control: cellular mechanisms of plasma membrane repair.
- Scheffer et al. 2014 and Shukla et al. 2022 for ESCRT-mediated repair.
- Muratori et al. 2021 for ESCRT-III and Annexin V repair after nsPEF permeabilization.
- Williams et al. 2023 for Annexin A6 and calcium-dependent exosome secretion during plasma membrane repair.
- Messenger et al. 2018 for Ca2+-stimulated exosome release regulated by Munc13-4.
- Welsh et al. 2024, MISEV2023, for EV study standards.

Weekly modules

Week 1: Course map and systems thinking
Theme: Why electro-exocytosis needs both cell biology and modeling.
Concepts: mechanism, perturbation, state variable, output, model assumption.
Lab connection: run and modify a one-state ODE.

Week 2: Membranes and transport
Theme: The plasma membrane as a barrier, capacitor, and transport surface.
Concepts: lipid bilayer, permeability, channels, pumps, gradients, flux.
Lab connection: simple membrane flux and gradient model.

Week 3: Organelles, calcium, and bioenergetics
Theme: Calcium and organelles connect membrane perturbation to cell state.
Concepts: ER, mitochondria, lysosomes, calcium buffering, ROS, ATP.
Lab connection: calcium pulse and clearance model.

Week 4: Vesicle trafficking, exocytosis, and endocytosis
Theme: Cells route material through regulated membrane traffic.
Concepts: vesicle budding, Rab proteins, SNAREs, MVBs, fusion, recycling.
Lab connection: vesicle-pool turnover model.

Week 5: Plasma membrane repair and extracellular vesicles
Theme: Repair and vesicle release are related cellular responses to membrane damage.
Concepts: lysosomal exocytosis, annexins, ESCRT, shedding, exosomes, microvesicles.
Lab connection: calcium-triggered repair-state model.

Week 6: Compartmental modeling and ODEs
Theme: A compartment is a well-mixed bookkeeping device for mass balance.
Concepts: states, fluxes, conservation, rate laws, numerical integration.
Lab connection: two-compartment exchange model in Python/SciPy.

Week 7: Nonlinearity, feedback, and model credibility
Theme: Biological models are useful when their assumptions and failure modes are visible.
Concepts: feedback, thresholds, regimes, sensitivity, parameter uncertainty.
Lab connection: sensitivity scan and regime classification.

Week 8: Bioelectricity and membrane charging
Theme: Electric fields become cellular perturbations through membrane voltage and charging.
Concepts: resting potential, capacitance, transmembrane voltage, cell geometry.
Lab connection: reduced membrane charging model.

Week 9: Electroporation, nanoporation, and nsPEF dose
Theme: Pulse parameters shape membrane permeabilization, intracellular effects, and heating.
Concepts: pulse width, amplitude, pulse number, repetition rate, absorbed energy, thermal confounding.
Lab connection: pore-state toy model with pulse-parameter sweep.

Week 10: nsPEF-induced cellular responses
Theme: A pulse can trigger calcium, organelle stress, ROS, ATP changes, and cytoskeletal remodeling.
Concepts: calcium entry, ER release, mitochondrial perturbation, ROS generation, ATP depletion.
Lab connection: coupled calcium-ROS-ATP model.

Week 11: Coupling repair, remodeling, and EV release
Theme: Productive EV release may sit between mild reversible perturbation and injury.
Concepts: repair capacity, vesicle pool recruitment, blebbing, EV subtype logic, regime maps.
Lab connection: repair-to-EV release model.

Week 12: Minimal electro-exocytosis framework
Theme: Build the smallest model that connects nsPEF input to EV output.
Concepts: model boundary, state selection, input function, observable output, modular coupling.
Lab connection: assemble membrane, calcium, repair, and EV modules.

Week 13: Calibration, validation, and reproducibility
Theme: A useful model must state what evidence can support or falsify it.
Concepts: literature-derived parameters, synthetic data, fitting, uncertainty, reproducible notebooks.
Lab connection: fit or calibrate one submodule and document uncertainty.

Week 14: Final presentations and peer critique
Theme: Students present the model, not just the code.
Deliverables: model diagram, equations, simulation results, sensitivity analysis, limitations, future measurements.

Computational lab arc
1. One-state ODE and exponential recovery.
2. Membrane flux and transport model.
3. Calcium pulse model with buffering and clearance.
4. Vesicle pool turnover model.
5. Membrane charging and pore-state toy model.
6. Repair-state and EV-release model.
7. Regime map and sensitivity analysis.
8. SBML/COPASI/VCell/BioModels comparison or export.
9. Final project model integration.

Possible final project
Students build a minimal systems-biology model of nsPEF-induced electro-exocytosis. The model should include an nsPEF input function, membrane perturbation or pore-state variable, cytosolic Ca2+ state, optional organelle-stress state, repair/remodeling state, vesicle or EV pool state, and EV-release output. Students choose a biological question, justify assumptions, implement simulations, run sensitivity analysis, compare at least two pulse conditions, and explain what experimental data would be needed for calibration.

Suggested final project deliverables
- One-page specific aims.
- Model diagram and state table.
- Equations or reaction/transition rules.
- Parameter table with source or assumption labels.
- Reproducible notebook or model file.
- Simulation figures for baseline and nsPEF perturbation.
- Sensitivity or uncertainty analysis.
- Short written report and oral presentation.

First-pass public resource links
- NCBI Bookshelf, Molecular Biology of the Cell: https://www.ncbi.nlm.nih.gov/books/NBK21054/
- OpenStax Biology 2e: https://openstax.org/details/books/biology-2e
- HHMI BioInteractive classroom resources: https://www.biointeractive.org/classroom-resources
- EMBL-EBI Training: https://www.ebi.ac.uk/training/
- MIT OCW 7.91J Foundations of Computational and Systems Biology: https://ocw.mit.edu/courses/7-91j-foundations-of-computational-and-systems-biology-spring-2014/
- SciPy solve_ivp: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
- SBML: https://sbml.org/
- BioModels: https://www.ebi.ac.uk/biomodels/
- VCell: https://vcell.org/
- COPASI: http://www.copasi.org/
`;

const resourceNotes = `Resource starting points
========================

Teaching foundations
- NCBI Bookshelf / Molecular Biology of the Cell: open background chapters for membranes, organelles, cytoskeleton, signaling, and vesicle traffic.
- OpenStax Biology 2e: accessible review for mixed-background students.
- HHMI BioInteractive: short visual assets for introductory mechanisms.
- EMBL-EBI Training: curated biology and bioinformatics training material.
- FEBS review journals: candidates for short, advanced-undergraduate review readings.
- Local BIEN5700 and BIEN6931 materials: pacing, physiology-first organization, compartmental modeling, ODEs, PK/PD, PBPK, and project expectations.

Computational platforms
- Python/SciPy/Jupyter: transparent baseline workflow for every student.
- COPASI: biochemical reaction networks and deterministic/stochastic simulation.
- VCell: compartmental and spatial cell modeling.
- SBML: exchange format for reusable model definitions.
- BioModels: curated examples and model provenance.
- CellOrganizer: candidate spatial-cell modeling resource; web access should be checked in the next pass.

Specialized literature spine from the current manuscript bibliography
- Kotnik et al. 2019, Annual Review of Biophysics, doi:10.1146/annurev-biophys-052118-115451.
- Krassowska and Filev 2007, Biophysical Journal, doi:10.1529/biophysj.106.094235.
- Smith and Weaver 2004, Biophysical Journal, doi:10.1529/biophysj.104.040733.
- Beebe et al. 2003, FASEB Journal, doi:10.1096/fj.02-0859fje.
- Vernier et al. 2003, BBRC, doi:10.1016/j.bbrc.2003.08.140.
- Orlacchio et al. 2023, IJMS, doi:10.3390/ijms241914999.
- Reddy et al. 2001, Cell, doi:10.1016/S0092-8674(01)00421-4.
- McNeil and Steinhardt 2003, Annual Review of Cell and Developmental Biology, doi:10.1146/annurev.cellbio.19.111301.140101.
- Andrews et al. 2014, Trends in Cell Biology, doi:10.1016/j.tcb.2014.07.008.
- Scheffer et al. 2014, Nature Communications, doi:10.1038/ncomms6646.
- Muratori et al. 2021, Bioelectrochemistry, doi:10.1016/j.bioelechem.2021.107837.
- Williams et al. 2023, eLife, doi:10.7554/eLife.86556.
- Messenger et al. 2018, Journal of Cell Biology, doi:10.1083/jcb.201710132.
- Welsh et al. 2024, MISEV2023, Journal of Extracellular Vesicles, doi:10.1002/jev2.12404.
`;

async function writeBlob(filePath, blob) {
  await fs.writeFile(filePath, new Uint8Array(await blob.arrayBuffer()));
}

async function getRecords(presentation) {
  const inspected = await presentation.inspect({
    kind: "slide,textbox,shape,notes",
    maxChars: 200000,
  });
  return inspected.ndjson
    .trim()
    .split(/\n+/)
    .filter(Boolean)
    .map((line) => JSON.parse(line));
}

function findRecord(records, slide, predicate) {
  return records.find((record) => record.slide === slide && predicate(record));
}

function addTextbox(slide, name, position, text, style = {}) {
  const shape = slide.shapes.add({
    geometry: "textbox",
    name,
    position,
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  shape.text = text;
  shape.text.style = {
    typeface: "Arial",
    fontSize: 22,
    color: TEXT_DARK,
    ...style,
  };
  return shape;
}

function addFooterOverlay(slide, slideNumber) {
  slide.shapes.add({
    geometry: "rect",
    name: `footer-cover-${slideNumber}`,
    position: { left: 258, top: 670.8, width: 935, height: 49 },
    fill: ODU_BLUE,
    line: { style: "solid", fill: ODU_BLUE, width: 0 },
  });
  addTextbox(
    slide,
    `footer-text-${slideNumber}`,
    { left: 274, top: 689, width: 610, height: 18 },
    footerText,
    { fontSize: 10, color: "#FFFFFF" },
  );
  slide.shapes.add({
    geometry: "line",
    name: `footer-separator-${slideNumber}`,
    position: { left: 1193.5, top: 676, width: 0, height: 38.3 },
    fill: "none",
    line: { style: "solid", fill: "#FFFFFF", width: 1 },
  });
  slide.shapes.add({
    geometry: "rect",
    name: `slide-number-cover-${slideNumber}`,
    position: { left: 1195, top: 670.8, width: 85, height: 49 },
    fill: ODU_BLUE,
    line: { style: "solid", fill: ODU_BLUE, width: 0 },
  });
  addTextbox(
    slide,
    `slide-number-${slideNumber}`,
    { left: 1211, top: 688, width: 45, height: 20 },
    String(slideNumber),
    { fontSize: 12, color: "#FFFFFF" },
  );
}

function fillLocalTemplatePlaceholders(slide, slideNumber) {
  for (const shape of slide.shapes.items ?? []) {
    if (shape.name === "Footer Placeholder 1") {
      shape.text = footerText;
      shape.text.style = { typeface: "Arial", fontSize: 10, color: "#FFFFFF" };
    }
    if (shape.name === "Slide Number Placeholder 2") {
      shape.text = String(slideNumber);
      shape.text.style = { typeface: "Arial", fontSize: 12, color: "#FFFFFF" };
    }
  }
}

function addContentSlide(slide, spec, slideNumber) {
  addTextbox(
    slide,
    `title-${slideNumber}`,
    { left: 70, top: 86, width: 1080, height: 58 },
    spec.title.toUpperCase(),
    { typeface: "Arial Black", fontSize: 34, bold: true, color: ODU_BLUE },
  );

  slide.shapes.add({
    geometry: "rect",
    name: `accent-line-${slideNumber}`,
    position: { left: 70, top: 153, width: 188, height: 5 },
    fill: ODU_RED,
    line: { style: "solid", fill: ODU_RED, width: 0 },
  });

  addTextbox(
    slide,
    `kicker-${slideNumber}`,
    { left: 70, top: 184, width: 1090, height: 52 },
    spec.kicker,
    { fontSize: 25, bold: true, color: TEXT_DARK },
  );

  const bodyText = spec.body.map((item) => `- ${item}`).join("\n");
  addTextbox(
    slide,
    `body-${slideNumber}`,
    { left: 88, top: 255, width: 1065, height: 360 },
    bodyText,
    { fontSize: spec.body.length > 5 ? 19 : 21, color: TEXT_DARK },
  );

  slide.speakerNotes.textFrame.setText(spec.note);
  slide.speakerNotes.setVisible(true);
  addFooterOverlay(slide, slideNumber);
}

async function main() {
  await fs.mkdir(finalDir, { recursive: true });
  await fs.mkdir(previewDir, { recursive: true });

  const presentation = await PresentationFile.importPptx(
    await FileBlob.load(starterPptx),
  );
  const records = await getRecords(presentation);

  const titleSlide = presentation.slides.items[0];
  const titleRecord = findRecord(records, 1, (r) => r.placeholder === "title");
  const subtitleRecord = findRecord(
    records,
    1,
    (r) => r.placeholder === "subtitle",
  );
  const affiliationRecord = findRecord(
    records,
    1,
    (r) => r.name === "Content Placeholder 4",
  );
  const pictureRecord = findRecord(
    records,
    1,
    (r) => r.placeholder === "picture",
  );

  if (titleRecord) {
    const titleShape = presentation.resolve(titleRecord.id);
    titleShape.text =
      "COMPUTATIONAL SYSTEMS BIOLOGY\nOF ELECTRO-EXOCYTOSIS";
  }
  if (subtitleRecord) {
    const subtitleShape = presentation.resolve(subtitleRecord.id);
    subtitleShape.text = "Draft syllabus for discussion";
  }
  if (affiliationRecord) {
    const affiliationShape = presentation.resolve(affiliationRecord.id);
    affiliationShape.text =
      "Senior undergraduate / graduate course | Cell biology, bioelectrics, and computational modeling";
  }
  if (pictureRecord) {
    presentation.resolve(pictureRecord.id).delete();
  }
  titleSlide.speakerNotes.textFrame.setText(
    "First-pass syllabus draft using AG_template.pptx. The next pass should select specific readings and figures for each module.",
  );
  titleSlide.speakerNotes.setVisible(true);

  const freshRecords = await getRecords(presentation);
  for (let index = 1; index < presentation.slides.items.length; index += 1) {
    const slideNumber = index + 1;
    const slide = presentation.slides.items[index];
    const footer = findRecord(
      freshRecords,
      slideNumber,
      (r) => r.placeholder === "footer",
    );
    const number = findRecord(
      freshRecords,
      slideNumber,
      (r) => r.placeholder === "slideNumber",
    );
    if (footer) {
      const footerShape = presentation.resolve(footer.id);
      footerShape.text = footerText;
      footerShape.text.style = { typeface: "Arial", fontSize: 10, color: "#FFFFFF" };
    }
    if (number) {
      const numberShape = presentation.resolve(number.id);
      numberShape.text = String(slideNumber);
      numberShape.text.style = { typeface: "Arial", fontSize: 12, color: "#FFFFFF" };
    }
    fillLocalTemplatePlaceholders(slide, slideNumber);
    addContentSlide(slide, contentSlides[index - 1], slideNumber);
  }

  const finalInspect = await presentation.inspect({
    kind: "slide,textbox,shape,image,table,chart,notes",
    maxChars: 300000,
  });
  await fs.writeFile(`${finalPptx}.inspect.ndjson`, finalInspect.ndjson);

  for (const [index, slide] of presentation.slides.items.entries()) {
    const stem = `slide-${String(index + 1).padStart(2, "0")}`;
    await writeBlob(
      path.join(previewDir, `${stem}.png`),
      await presentation.export({ slide, format: "png", scale: 1 }),
    );
    await fs.writeFile(
      path.join(previewDir, `${stem}.layout.json`),
      await (await slide.export({ format: "layout" })).text(),
    );
  }
  await writeBlob(
    path.join(previewDir, "syllabus-montage.webp"),
    await presentation.export({ format: "webp", montage: true, scale: 1 }),
  );

  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(finalPptx);

  await fs.writeFile(path.join(finalDir, "course_syllabus_draft.txt"), detailedSyllabus);
  await fs.writeFile(path.join(finalDir, "resource_starting_points.txt"), resourceNotes);
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
