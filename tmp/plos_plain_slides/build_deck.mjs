import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { Presentation, PresentationFile } from '@oai/artifact-tool';

const ROOT = '/Users/aghorban/code/electroexo';
const PAPER = '/Users/aghorban/code/manuscripts/electroexocytosis-computational-plos-computational-biology';
const SKILL = '/Users/aghorban/.codex/plugins/cache/openai-primary-runtime/presentations/26.1013.11914/skills/presentations';
const BUILD = path.join(ROOT, 'tmp/plos_plain_slides');
const OUT = path.join(ROOT, 'output/presentations');
const PYTHON = '/Users/aghorban/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
const { resolvePresentationFont, applyPresentationChartFont, finalizePresentation } = await import(pathToFileURL(path.join(SKILL, 'container_tools/artifact_tool_utils.mjs')));
const FONT = resolvePresentationFont();
const INK = '#17222E', BLUE = '#315E7B', MUTED = '#53616C', LIGHT = '#D8DEE3';
const pres = Presentation.create({ slideSize: { width: 1280, height: 720 } });
await fs.mkdir(OUT, { recursive: true });
await fs.mkdir(path.join(BUILD, 'renders'), { recursive: true });

function text(slide, str, x, y, w, h, size=30, bold=false, color=INK, align='left') {
  const shape = slide.shapes.add({ geometry:'textbox', name:str.slice(0,60), position:{left:x,top:y,width:w,height:h}, fill:'none', line:{fill:'none',width:0} });
  shape.text = str;
  shape.text.style = { typeface:FONT, fontSize:size, bold, color, alignment:align, verticalAlignment:'top', autoFit:'none', insets:0 };
  return shape;
}
function slide(title, subtitle='') {
  const s=pres.slides.add(); s.background.fill='#FFFFFF';
  text(s,title,64,44,1152,70,43,true);
  if(subtitle) text(s,subtitle,64,124,1152,58,26,false,MUTED);
  text(s,String(pres.slides.items.length),1180,675,36,25,19,false,MUTED,'right');
  return s;
}
async function pic(s,file,alt,box,crop) {
  const bytes=await fs.readFile(file);
  s.images.add({blob:bytes,contentType:'image/png',alt,fit:'contain',position:box,...(crop?{crop}:{})});
}
function notes(s,script,sources,cue='') {
  s.speakerNotes.text=`SPEAKER SCRIPT\n\n${script}\n\nVISUAL CUE\n${cue}\n\nSOURCES\n${sources.join('\n')}`;
}
const scripts=[];
function record(title,script,cue,source) { scripts.push({title,script,cue,source}); }

// 1. A substantive opening rather than a separate cover.
{
  const title='Modeling electro-exocytosis';
  const s=slide(title,'A computational framework for electrical stimulation of extracellular vesicle release');
  text(s,'How do short electric pulses change\nwhat cells release?',64,235,1110,126,48,true,BLUE);
  text(s,'Particle counts can reflect regulated secretion,\nmembrane repair, or cell injury.',64,405,1080,96,34);
  text(s,'Our contribution: connect the cellular mechanisms to the particles an assay measures.',64,555,1100,82,30);
  const script='The central question in this project is how very short electric pulses change the particles that living cells release. A higher particle count can have several explanations. Cells may increase regulated vesicle secretion, shed membrane during repair, or release particles as they become injured. Those outcomes have different biological and engineering meanings. We developed a computational framework that follows the response from electrical exposure through the cell and into the surrounding medium, then connects it to the particle assay. In this talk, I will show what we built, the illustrative behavior it produces, and what the preliminary fit to particle measurements can and cannot tell us.';
  const sources=[`${PAPER}/main.tex, Author summary and Introduction; lines 216–218, 255–275`];
  const cue='Open with the biological question. Define extracellular vesicles as EVs, then keep the focus on the measured particles.';
  notes(s,script,sources,cue); record(title,script,cue,sources);
}

// 2. Eight layers condensed into four readable groups, with the observation bridge explicit.
{
  const title='The framework we developed';
  const s=slide(title,'Eight linked model layers, with cell-state modifiers across the system');
  const groups=[
    ['Electrical exposure','Pulse delivery and dose\nMembrane and organelle response'],
    ['Cell response','Ion transport, calcium, ROS and ATP\nMembrane remodeling and repair'],
    ['Vesicle production','EV pathways and release\nCargo and potency proxies'],
    ['Quality and process','Injury and quality gating\nIsolation and quality-control interfaces'],
  ];
  groups.forEach(([label,body],i)=>{
    const y=205+i*101;
    text(s,label,64,y,285,42,29,true,BLUE);
    text(s,body,365,y,595,78,28);
  });
  text(s,'Measurement bridge',970,205,246,80,29,true,BLUE);
  text(s,'Release from cells\n\nAccumulation and loss\nin the medium\n\nMeasured particle\nsize distribution',970,293,246,295,27);
  text(s,'Python simulator, reproducible parameter snapshots, calibration and sensitivity workflows',64,640,1120,34,23,false,MUTED);
  const script='The software uses a reduced, modular description rather than a whole-cell simulation. The eight layers cover electrical dose and membrane response, ion transport and cellular stress, repair, vesicle release, and downstream quality and process outputs. Cell-state modifiers allow us to vary the shared model for different producer-cell settings. The reduced equations carry information from the short pulse into responses lasting seconds to hours. The crucial additional piece is the measurement bridge on the right. Cellular release first enters an extracellular stock, where particles can disappear, and only then passes through an observation model. This keeps secretion distinct from what the instrument reports. Cargo, potency and process outputs remain exploratory model proxies.';
  const sources=[`${PAPER}/main.tex, Materials and methods; lines 277–340 and module map at 639–647`,`${ROOT}/README.md, Modular solver architecture and Output files`];
  const cue='Read the four model groups briefly, then spend most of the slide on the release, medium and measurement distinction. ROS means reactive oxygen species; ATP denotes cellular energy status.';
  notes(s,script,sources,cue); record(title,script,cue,sources);
}

// 3. Editable numerical summary of the same scenarios plotted in paper Fig 7.
{
  const title='Simulations distinguish secretion from injury';
  const s=slide(title,'Illustrative scenarios from the paper, with deliberately different model settings');
  const values=[
    ['Scenario','Small-EV output','Apoptotic-body output','Viability proxy'],
    ['Mild','1.00','1.00','76%'],
    ['Secretory','2.67','1.37','66%'],
    ['Injury-dominant','0.67','219.49','0%'],
  ];
  const tab=s.tables.add({rows:4,columns:4,left:64,top:220,width:1152,height:272,columnWidths:[330,264,315,243],values});
  tab.styleOptions={headerRow:false,bandedRows:false};
  tab.borders.assign({fill:LIGHT,width:1,style:'solid'});
  for(let r=0;r<4;r++) for(let c=0;c<4;c++) {
    const cell=tab.getCell(r,c); cell.fill='#FFFFFF';
    cell.text.style={typeface:FONT,fontSize:27,bold:r===0,color:r===0?BLUE:INK,alignment:c===0?'left':'center',verticalAlignment:'middle',insets:12};
  }
  text(s,'Cumulative output over 30 min, normalized within each particle pathway. Viability is absolute.',64,515,1140,55,23,false,MUTED);
  text(s,'The injury scenario produces many particles\nwhile the model’s viability proxy collapses.',64,578,1110,80,32,true,BLUE);
  const script='These are constructed model scenarios, not experimental treatment groups. They illustrate how the framework distinguishes different routes to particle release over a thirty-minute simulation. To make the cumulative outputs readable, each particle-pathway column is divided by its own value in the mild scenario. Viability remains an absolute model percentage. In the secretory scenario, small-EV output is about 2.7 times the mild value. In the injury scenario, that output falls, while apoptotic-body output rises by more than two hundred times and the viability proxy reaches zero. The point is that more particles can accompany a very different quality outcome. We changed both exposure settings and selected biological parameters between these scenarios, so this comparison does not identify a pulse optimum or prove a dose-response relationship.';
  const sources=[`${PAPER}/main.tex, integrated preliminary storyboard; lines 443–456`,`${PAPER}/figures/multilayer_storyboard.png, Fig 7`,`${ROOT}/results/multilayer_storyboard/multilayer_storyboard_summary.csv`,`${ROOT}/examples/storyboard_multilayer_overview.py, SCENARIOS`];
  const cue='Compare within each output column. The two normalized columns do not share a common particle-count scale. Viability is a model proxy.';
  notes(s,script,sources,cue); record(title,script,cue,sources);
}

// 4. Exact middle row from the manuscript, including the visible 1-hour underfit.
{
  const title='The fit captures the main particle-size peak';
  const s=slide(title,'26 stimulated histograms across three exposure labels and three harvest times');
  text(s,'Example: 3 pulses, 40-kV source label',64,195,1152,42,29,true,BLUE);
  text(s,'0.5 h',240,252,210,34,26,true,INK,'center');
  text(s,'1 h',618,252,210,34,26,true,INK,'center');
  text(s,'3 h',987,252,210,34,26,true,INK,'center');
  await pic(s,path.join(PAPER,'figures/experimental_size_profile_fit.png'),'Middle row of manuscript Figure 9: three-pulse source-labeled condition at 0.5, 1 and 3 hours, observed means with sample SD and fitted curves.',{left:64,top:290,width:1152,height:253.5},{left:0.07,top:0.347,right:0,bottom:0.39});
  // The paper suppresses x tick labels above its bottom row. Restore the same diameter ticks here.
  for(const centers of [[143,244,344],[529,629,729],[914,1014,1115]]) {
    for(let k=0;k<3;k++) text(s,String((k+1)*100),centers[k]-25,548,50,25,20,false,MUTED,'center');
  }
  text(s,'Particle diameter (nm)',64,584,1152,30,24,false,INK,'center');
  text(s,'Blue: observed mean ± sample SD. Pink: fitted model. Vertical axes: 10⁹ particles/mL per bin.',64,632,1152,28,22,false,MUTED);
  text(s,'80–380 nm in 20-nm bands. Repeated-measurement independence and delivered field remain unconfirmed.',64,675,1100,26,19,false,MUTED);
  const script='The preliminary calibration uses size-resolved particle measurements from CD4-positive cell preparations. The full dataset has three exposure labels, harvests at half an hour, one hour and three hours, and 26 stimulated histograms. We analyzed fifteen common diameter bands from 80 to 380 nanometers. This slide shows the complete three-pulse condition across all three harvests. Blue points and bars are the observed means and sample standard deviations, and pink curves are the selected fit. The model captures the dominant peak around 100 to 160 nanometers, but it visibly misses part of the one-hour peak. There are three repeated measurements at the outer times and two at one hour. Their biological independence is unknown, and the voltage label is not a confirmed delivered field.';
  const sources=[`${PAPER}/main.tex, measurements, size kinetics and calibration; lines 292–294, 474–494`,`${PAPER}/figures/experimental_size_profile_fit.png, Fig 9, middle row; original points and error bars retained`];
  const cue='Point to the common peak, then the underprediction at 1 h. The repeated measurements are not established biological replicates. Y-axis ranges differ across panels, as in the source figure.';
  notes(s,script,sources,cue); record(title,script,cue,sources);
}

// 5. Editable evidence chart for nested calibration, plus two interpretable metrics.
{
  const title='An empirical exposure correction improves the fit';
  const s=slide(title,'In-sample calibration across nine condition–time targets');
  text(s,'Total-concentration error (log₁₀ RMSE)',64,207,730,36,28,true,BLUE);
  const ch=s.charts.add('bar',{
    position:{left:64,top:257,width:740,height:318},
    categories:['Fixed kernels','State-dependent kernels','With exposure correction'],
    series:[{name:'Concentration error',values:[0.1879,0.1879,0.1042],fill:'#6B8AA0',points:[{idx:2,fill:BLUE}]}],
    barOptions:{direction:'bar',grouping:'clustered',gapWidth:80},hasLegend:false,
    xAxis:{textStyle:{fontSize:23,fill:INK},line:{fill:LIGHT,width:1},majorGridlines:null},
    yAxis:{min:0,max:0.22,majorUnit:0.05,numberFormatCode:'0.00',textStyle:{fontSize:22,fill:MUTED},majorGridlines:{fill:LIGHT,width:1}},
    dataLabels:{showValue:true,numberFormatCode:'0.0000',position:'outEnd',textStyle:{fontSize:24,fill:INK}},
    chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF',chartLine:{fill:'none',width:0},
  }); applyPresentationChartFont(ch,{fontFamily:FONT});
  text(s,'1.27-fold',862,257,350,68,48,true,BLUE);
  text(s,'Typical concentration\ndiscrepancy in 80–380 nm',862,332,350,90,27);
  text(s,'4.66 nm',862,448,350,68,48,true,BLUE);
  text(s,'Mean-diameter\nabsolute error',862,521,350,76,27);
  text(s,'The correction signals missing condition dependence in the upstream model.',64,627,1120,42,27,true);
  const script='We compared three observation-model variants. The first uses fixed size kernels. The second lets their medians respond to the upstream cell state. Those variants give essentially the same total-concentration error. Adding a shared empirical exposure correction reduces the log-ten RMSE from 0.1879 to 0.1042. For the selected model, that is a typical multiplicative concentration discrepancy of 1.27-fold within the measured size window. The mean absolute error in histogram mean diameter is 4.66 nanometers, and the mean Hellinger distance for size composition is 0.0857. These are in-sample fit statistics. The improvement tells us that the current upstream model lacks some condition dependence. It does not establish the empirical correction as a biological mechanism or demonstrate prediction on new conditions.';
  const sources=[`${PAPER}/main.tex, Preliminary calibration; lines 488–494`,`${ROOT}/results/ffrci_size_resolved_fit_v1_1/model_comparison.csv`];
  const cue='Explain that smaller RMSE is better. The fold metric is 10 raised to the log10 RMSE. Concentration refers to measured-window particles, not total authenticated EV yield.';
  notes(s,script,sources,cue); record(title,script,cue,sources);
}

// 6. Original local-information diagnostic with a plain-language interpretation.
{
  const title='Some parameters remain poorly constrained';
  const s=slide(title,'Exploratory screening and local diagnostics answer different questions');
  text(s,'Local information near the fit (excerpt)',64,210,650,39,26,true,BLUE);
  await pic(s,path.join(BUILD,'identifiability.png'),'Lower rows of panel A of manuscript Figure 16, retaining the four AB parameters, local residual-curvature axis and near-zero-information threshold.',{left:64,top:260,width:650,height:277.3},{left:0,top:0.267,right:0.13,bottom:0.52});
  text(s,'All four apoptotic-body parameters\nfall below the local information\nthreshold.',750,222,465,130,30,true,BLUE);
  text(s,'Release scaling, extracellular loss\nand assay recovery remain\ndifficult to separate.',750,395,465,130,29);
  text(s,'198 runtime parameters and 16 observation parameters screened separately',64,635,1150,31,23,false,MUTED);
  text(s,'Dashed line: near-zero-information threshold. Curvature is a local diagnostic, not a confidence interval.',64,675,1100,26,19,false,MUTED);
  const script='We also asked which parameters affect the outputs and which parameters the fitted observations can constrain. Those are different questions. The exploratory sensitivity analysis screened 198 runtime parameters and 16 observation parameters in separate analyses. Source-conversion factors and effective particle persistence can influence the measured concentration. However, near the selected fit, all four apoptotic-body parameters fall below the local information threshold shown here. This is consistent with the limited measurement window and a fitted apoptotic contribution near its lower bound. More broadly, release scaling, extracellular loss and assay recovery remain difficult to disentangle. Recovery was provisionally fixed at one, with uncertainty absorbed into fitted source scales, rather than independently estimated. The analysis guides measurement choices. Its conditional ranges and small sampling design do not justify a definitive biological ranking, and the local curvature is not a confidence interval.';
  const sources=[`${PAPER}/main.tex, Global sensitivity and local constraint; lines 342–376, 543–562`,`${PAPER}/figures/global_sensitivity_identifiability.pdf, Fig 16, panel A, lower rows`];
  const cue='Point out the four AB labels below the dashed threshold. AB means apoptotic body; GSD means geometric standard deviation. A large curvature diagonal alone also does not prove individual identifiability.';
  notes(s,script,sources,cue); record(title,script,cue,sources);
}

// 7. The contribution and the experiments it motivates.
{
  const title='The framework makes the next measurements clear';
  const s=slide(title,'A model for experimental design and further calibration');
  const rows=[
    ['Exposure and assay controls','Confirm pulse delivery, culture volume and recovery.'],
    ['Release versus disappearance','Add matched shams and more harvests between 1 and 3 h.'],
    ['Biological identity and cell state','Pair particle profiles with EV characterization, calcium and injury/viability.'],
    ['Independent validation','Test the calibrated model on additional exposure conditions.'],
  ];
  rows.forEach(([label,body],i)=>{
    const y=212+i*88;
    text(s,label,64,y,436,42,28,true,BLUE);
    text(s,body,525,y,690,72,28);
  });
  text(s,'A working mechanistic framework, a preliminary fit,\nand explicit targets for predictive calibration',64,590,1150,78,34,true);
  const script='The practical outcome is a clear set of measurements for the next study. First, confirm the delivered pulse and the assay scaling, including medium volume and recovery controls. Second, add matched sham harvests and more observations between one and three hours, so we can better distinguish continued release from disappearance in the medium. Third, combine size profiles with orthogonal EV characterization and measurements of calcium and cell injury or viability. Those measurements connect particle behavior back to the biological pathways. Finally, reserve additional exposure conditions for independent validation. What we have now is a working mechanistic framework and a preliminary observation-level fit. It reproduces the dominant size behavior and exposes the measurements needed to develop a predictive model.';
  const sources=[`${PAPER}/main.tex, measurement priorities and discussion; lines 528, 555–562, 689–692`];
  const cue='End on the scientific contribution: the framework connects mechanisms to measurement and makes unresolved explanations testable. These measurement proposals are not a powered experimental plan.';
  notes(s,script,sources,cue); record(title,script,cue,sources);
}

// Script is also provided as plain text for rehearsal, including the precise sources and cues.
let scriptText='ELECTRO-EXOCYTOSIS: SEVEN PLAIN SLIDES\n\nAudience: computational biology research audience\nSuggested duration: approximately 8 minutes\nBasis: the supplied PLOS Computational Biology manuscript draft\n\n';
scripts.forEach((s,i)=>{
  scriptText+=`SLIDE ${i+1}: ${s.title}\n${'='.repeat(72)}\n\n${s.script}\n\nVisual cue: ${s.cue}\n\nSources:\n${s.source.map(x=>'  '+x).join('\n')}\n\n\n`;
});
scriptText+='FIVE- OR SIX-SLIDE OPTIONS\n\nFor six slides, merge the opening question into the framework slide.\nFor five slides, also combine the particle-size example and calibration comparison.\nKeep the distinction between illustrative simulation, in-sample calibration, and independent validation visible in either version.\n';
await fs.writeFile(path.join(OUT,'electroexocytosis_speaker_script.txt'),scriptText);
await fs.writeFile(path.join(BUILD,'content.json'),JSON.stringify(scripts,null,2));

for(let i=0;i<pres.slides.items.length;i++) {
  const s=pres.slides.items[i];
  const png=await pres.export({slide:s,format:'png',scale:1.5});
  await fs.writeFile(path.join(BUILD,'renders',`slide-${i+1}.png`),new Uint8Array(await png.arrayBuffer()));
}
const candidate=path.join(BUILD,'candidate_v3.pptx');
await (await PresentationFile.exportPptx(pres)).save(candidate);
const result=await finalizePresentation({
  workspaceDir:ROOT,candidatePath:candidate,finalPath:path.join(OUT,'electroexocytosis_seven_plain_slides_v3.pptx'),
  explicitTotalSlideCount:7,requiredNativeTableOwnerSlides:[3],requiredNativeChartOwnerSlides:[5],materializeLiteralChartWorkbooks:true,
  pythonExecutable:PYTHON,integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),
  layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),
  layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit','--require-native-table-slide','3'],
  fontPolicy:{basis:'design',families:[FONT]},verifyArtifactToolImport:true,
  receiptPath:path.join(BUILD,'validation_v3.json'),
});
console.log(JSON.stringify({font:FONT,result}));
