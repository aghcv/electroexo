import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {createHash} from 'node:crypto';
import {FileBlob,PresentationFile} from '@oai/artifact-tool';

const ROOT='/Users/aghorban/code/electroexo';
const BUILD=ROOT+'/tmp/eight_module_slides';
const SOURCE=BUILD+'/source.pptx';
const OUT=ROOT+'/output/presentations';
const RESULTS=ROOT+'/results/eight_module_storyboards';
const SKILL='/Users/aghorban/.codex/plugins/cache/openai-primary-runtime/presentations/26.1013.11914/skills/presentations';
const PYTHON='/Users/aghorban/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
const {finalizePresentation,applyPresentationChartFont}=await import(pathToFileURL(SKILL+'/container_tools/artifact_tool_utils.mjs'));
const manifest=JSON.parse(await fs.readFile(RESULTS+'/manifest.json','utf8'));
const pres=await PresentationFile.importPptx(await FileBlob.load(SOURCE));
const originalIds=pres.slides.items.map(s=>s.id);
const FONT='Helvetica Neue',INK='#17222E',BLUE='#315E7B',MUTED='#53616C',MODEL='#4477AA';
const labels=[
 '1 Pulse delivery\n& dosimetry',
 '2 Membrane & organelle\nelectrodynamics',
 '3 Ion transport,\nCa²⁺, ROS & energy',
 '4 Remodeling\n& repair',
 '5 EV biogenesis\n& subtype release',
 '6 Cargo, composition\n& potency',
 '7 Injury, debris\n& quality',
 '8 Manufacturing,\nisolation & QC',
];
const metricLabels=[
 'Dose (J/mL)','Voltage proxy (V)','Calcium (µM)','Repair activation',
 'Small EVs (rel. units)','Protein cargo score','Viability score (%)','Recovered sEVs (rel.)',
];
const inputs=['Pulse settings','From 1','From 2','From 3','From 2, 3, 4, 7','From 3, 5','From 3','From 5 + recovery'];
const slideTitles=[
 'Mild, secretory and injury responses',
 'What changes when we add more pulses?',
 'Calcium clearance changes the cell response',
 'Vesicle machinery can change release and cargo',
 'Isolation efficiency changes recovered yield',
];
// Portable Excel chart caches retain 12 significant figures. Full simulation
// precision is preserved in the accompanying CSVs and manifest.
const portable=v=>Number(Number(v).toPrecision(12));
function text(s,str,x,y,w,h,size=20,bold=false,color=INK,align='left') {
 const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 sh.text=str;
 sh.text.style={typeface:FONT,fontSize:size,bold,color,alignment:align,autoFit:'none',insets:0};
 return sh;
}
function numFormat(max) {
 if(max>=10)return '0';
 if(max>=1)return '0.0';
 if(max>=.1)return '0.00';
 if(max>=.01)return '0.000';
 return '0.0E+00';
}
function axisFormat(max) {
 const half=max/2;
 if(Number.isInteger(half))return '0';
 if(max>=10)return '0.0';
 if(max>=1)return Number.isInteger(half*10)?'0.0':'0.00';
 return max>=.02?'0.000':'0.0000';
}
function valueFormat(v) {
 if(v===0)return '0';
 if(Math.abs(v)<.1)return '0.0000';
 if(Math.abs(v)<10)return '0.00';
 return '0.0';
}
const newSlides=[];
const footers=[];
for(let si=0;si<manifest.studies.length;si++) {
 const study=manifest.studies[si];
 const slide=pres.slides.add();slide.background.fill='#FFFFFF';
 newSlides.push(slide);
 text(slide,slideTitles[si],32,24,1210,56,37,true);
 text(slide,study.takeaway,32,88,1210,52,23,false,BLUE);
 const left=100,cellW=144;
 for(let c=0;c<8;c++) {
  const x=left+c*cellW;
  text(slide,labels[c],x+1,149,cellW-3,59,17,true,INK,'center');
  text(slide,inputs[c],x+1,207,cellW-3,20,12,false,MUTED,'center');
  text(slide,metricLabels[c],x+1,231,cellW-3,22,13,true,BLUE,'center');
 }
 for(let r=0;r<study.rows.length;r++) {
  const row=study.rows[r];
  const y=258+r*137;
  text(slide,row.label,5,y+27,94,86,17,true,INK);
  for(let c=0;c<8;c++) {
   const module=manifest.modules[c],panel=row.panels[c];
   const isLine=module.kind==='line';
   const x=left+c*cellW;
   const max=module.y_max;
   const axisBase={textStyle:{fontSize:11.5,fill:MUTED},line:{fill:'#AAB3BA',width:.7},majorGridlines:null};
   const xAxis={...axisBase,visible:isLine,tickLabelPosition:isLine?'nextTo':'none'};
   if(isLine) {
    // Leave room beyond the final tick so its label stays inside the chart frame.
    xAxis.min=module.x_min??0;xAxis.max=portable(module.x_max*1.1);
    xAxis.majorUnit=(module.x_max-(module.x_min??0))/2;
    xAxis.numberFormatCode='0';
    if(r===2)xAxis.title={text:module.x_label,textStyle:{fontSize:11,fill:MUTED}};
   }
   const series={name:'Model output',values:(isLine?panel.y:[panel.value??panel.y?.[0]]).map(portable),fill:MODEL};
   if(isLine) {
    series.xValues=panel.x.map(portable);
    series.line={fill:MODEL,width:1.5,style:'solid'};
    series.marker={symbol:'none'};
   }
   const chart=slide.charts.add(isLine?'scatter':'bar',{
    position:{left:x,top:y,width:cellW-7,height:126},
    categories:isLine?undefined:['Model'],series:[series],hasLegend:false,
    ...(isLine?{scatterOptions:{style:'line'}}:{barOptions:{direction:'column',grouping:'clustered',gapWidth:170}}),
    xAxis,
    yAxis:{...axisBase,min:module.y_min??0,max,majorUnit:(max-(module.y_min??0))/2,numberFormatCode:axisFormat(max),majorGridlines:{fill:'#E2E6EA',width:.55}},
    ...(isLine?{}:{dataLabels:{showValue:true,position:'outEnd',numberFormatCode:valueFormat(panel.value??panel.y?.[0]),textStyle:{fontSize:11.5,fill:INK}}}),
    chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF',chartLine:{fill:'none',width:0},
   });
   applyPresentationChartFont(chart,{fontFamily:FONT});
  }
 }
 text(slide,'Illustrative simulations. Cargo, viability and recovered output are model scores.',32,685,1175,23,16,false,MUTED);
 footers.push(text(slide,String(originalIds.length+si+1),1224,685,30,23,16,false,MUTED,'right'));
 slide.speakerNotes.text=`SPEAKER SCRIPT\n\n${study.script}\n\nWHAT CHANGES\n${study.rows.map(r=>r.label+': '+r.detail).join('\n')}\n\nFIGURE CAPTION\n${study.caption}\n\nMODULES AND LINKS\n${manifest.modules.map(m=>`${m.number}. ${m.title}. Readout: ${m.metric_label} (${m.unit}). Inputs: ${Array.isArray(m.inputs)?m.inputs.join(', '):m.inputs}.`).join('\n')}\n\nINTERPRETATION\nBars in modules 1 and 2 summarize exposure and charging. Modules 3 and 4 show the first 60 s. Module 5 shows cumulative small-EV release over 30 min in relative pathway units. Cargo is an exploratory endpoint score, viability uses peak simulated damage, and manufacturing applies a recovery multiplier to cumulative release rather than harvesting the extracellular stock. Couplings branch: column order is not a strictly serial pathway. Heating is diagnostic and does not drive cellular injury in this implementation. Repair activation does not dynamically reseal pores or change the damage-repair coefficient.\n\nSOURCES\n${ROOT}/examples/storyboard_eight_modules.py\n${RESULTS}/manifest.json\n${RESULTS}/storyboard_summary.csv\n`;
 await fs.writeFile(BUILD+`/study-${si+1}.layout.json`,await(await slide.export({format:'layout'})).text());
 const png=await pres.export({slide,format:'png',scale:1.5});
 await fs.writeFile(BUILD+`/study-${si+1}.png`,new Uint8Array(await png.arrayBuffer()));
 console.log(`Prepared study ${si+1}: ${study.title}`);
}

const hash=createHash('sha256').update(await fs.readFile(SOURCE)).digest('hex');
const finalOptions={
 workspaceDir:ROOT,pythonExecutable:PYTHON,
 integrityValidatorPath:SKILL+'/container_tools/inspect_presentation_package_integrity.py',
 layoutValidatorPath:SKILL+'/container_tools/inspect_presentation_layout_geometry.py',
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit'],
 fontPolicy:{basis:'reference',families:[FONT],referencePath:SOURCE,referenceSha256:hash},
 materializeLiteralChartWorkbooks:true,verifyArtifactToolImport:true,
};
const combinedPath=OUT+'/electroexocytosis_with_eight_module_studies.pptx';
let hasCombined=false;
try { await fs.access(combinedPath);hasCombined=true; } catch {}
if(!hasCombined) {
 const merged=BUILD+'/candidate_combined_v3.pptx';
 await(await PresentationFile.exportPptx(pres)).save(merged);
 const full=await finalizePresentation({...finalOptions,candidatePath:merged,finalPath:combinedPath,explicitTotalSlideCount:originalIds.length+5,requiredNativeChartOwnerSlides:[7,8,9,10,11],receiptPath:BUILD+'/validation_combined_v3.json'});
 console.log('Combined deck: '+full.finalPath);
} else {
 const receipt=JSON.parse(await fs.readFile(BUILD+'/validation_combined_v2.json','utf8'));
 const actual=createHash('sha256').update(await fs.readFile(combinedPath)).digest('hex');
 if(actual!==receipt.finalSha256)throw new Error('Previously finalized combined deck changed');
 console.log('Reusing validated combined deck.');
}
pres.slides.keep([6,7,8,9,10]);
footers.forEach((s,i)=>s.text=String(i+1));
const standalone=BUILD+'/candidate_five_studies_v2.pptx';
await(await PresentationFile.exportPptx(pres)).save(standalone);
const five=await finalizePresentation({...finalOptions,candidatePath:standalone,finalPath:OUT+'/electroexocytosis_five_eight_module_studies.pptx',explicitTotalSlideCount:5,requiredNativeChartOwnerSlides:[1,2,3,4,5],receiptPath:BUILD+'/validation_five_studies_v2.json'});
console.log('Study deck: '+five.finalPath);
const script=manifest.studies.map((s,i)=>`SLIDE ${i+1}: ${s.title}\n\n${s.script}\n\nScenario details:\n${s.rows.map(r=>r.label+': '+r.detail).join('\n')}\n\n${s.caption}`).join('\n\n\n');
await fs.writeFile(OUT+'/electroexocytosis_eight_module_study_scripts.txt',script);
