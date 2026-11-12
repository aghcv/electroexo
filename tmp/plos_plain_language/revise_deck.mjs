import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {FileBlob, PresentationFile} from '@oai/artifact-tool';

const ROOT='/Users/aghorban/code/electroexo';
const BUILD=ROOT+'/tmp/plos_plain_language';
const OUT=ROOT+'/output/presentations';
const SOURCE=OUT+'/electroexocytosis_seven_plain_slides_v3.pptx';
const PAPER='/Users/aghorban/code/manuscripts/electroexocytosis-computational-plos-computational-biology';
const SKILL='/Users/aghorban/.codex/plugins/cache/openai-primary-runtime/presentations/26.1013.11914/skills/presentations';
const PYTHON='/Users/aghorban/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
const {finalizePresentation}=await import(pathToFileURL(SKILL+'/container_tools/artifact_tool_utils.mjs'));
const pres=await PresentationFile.importPptx(await FileBlob.load(SOURCE));
const FONT='Helvetica Neue', INK='#17222E', BLUE='#315E7B', MUTED='#53616C';

function change(id,value,style={},position) {
 const shape=pres.resolve(id);
 shape.text=value;
 shape.text.style={typeface:FONT,autoFit:'none',insets:0,...style};
 if(position)shape.position=position;
 return shape;
}
function addText(s,value,x,y,w,h,size=28,bold=false,color=INK) {
 const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 sh.text=value;
 sh.text.style={typeface:FONT,fontSize:size,bold,color,insets:0,autoFit:'none',verticalAlignment:'top'};
 return sh;
}
const titles=[
 'Electro-exocytosis: a model to guide experiments',
 'The biological steps included in the model',
 'Linking cell responses to laboratory measurements',
 'The model captures the main particle-size pattern',
 'The model helps us choose what to measure next',
 'A focused plan for the next experiment',
];

// Preserve the user's six-slide order and overview figure.
change('sh/sryl4zqx',titles[0],{fontSize:43,bold:true,color:INK});
change('sh/7qp4be9c','How do short electric pulses change EV release?\nHow can models guide the next experiment?',{fontSize:44,bold:true,color:BLUE},{left:64,top:235,width:1152,height:165});

change('sh/xk7qlczu',titles[1],{fontSize:43,bold:true,color:INK});

change('sh/a5c3ql8z',titles[2],{fontSize:43,bold:true,color:INK});
change('sh/cza94vmx','The model follows the steps that can affect the final particle count',{fontSize:26,color:MUTED});
change('sh/298ryl4v','Pulse treatment');
change('sh/3ah8rqlg','Pulse settings and exposure\nEffects on cell membranes');
change('sh/pcjqtg36','Calcium, stress and energy\nMembrane damage and repair');
change('sh/m5cra54z','EV release');
change('sh/n6ls3alk','Vesicle formation and release\nEstimates of cargo and activity');
change('sh/n2l4fq98','Sample quality');
change('sh/m1c3mlsn','Cell injury and unwanted particles\nEffects of collection and isolation');
change('sh/943mhgre','What we measure');
change('sh/83ulovat','Particles leave cells\n\nSome disappear\nfrom the medium\n\nThe assay counts\nparticles by size');
change('sh/v6l4jq94','The model explores how release, particle loss and measurement shape the final count.',{fontSize:23,color:MUTED});

change('sh/oz29krqh',titles[3],{fontSize:43,bold:true,color:INK});
change('sh/1cj2d8b6','The comparison also shows where our explanation needs improvement',{fontSize:26,color:MUTED});
change('sh/rm1k7yt4','Example: 3 pulses, labeled 40 kV in the source data');
change('sh/ql8jytsj','30 min');
change('sh/doj29oba','1 hour');
change('sh/sna103ap','3 hours');
change('sh/h4fa9wfy','Blue: measured particles and variation between measurements. Pink: model estimates.',{fontSize:22,color:MUTED});
change('sh/p0batw72','Concentration: billions of particles/mL in each size band. Diameter: 80–380 nm.',{fontSize:19,color:MUTED});

// Replace the specialist diagnostic with the practical experimental meaning.
const s5=pres.resolve('sl/i107q5of');
change('sh/w32dkbuh',titles[4],{fontSize:43,bold:true,color:INK});
change('sh/dgbulwnm','We ask what could change EV counts and which measurements would explain the change',{fontSize:26,color:MUTED});
for(const id of ['sh/z2tcnm5s','im/vm1cv6p0','sh/l4bupwny','sh/032tgr6d','sh/n6dcr65o','sh/m5kbi1oj'])pres.delete(id);
const table=s5.tables.add({
 rows:4,columns:2,left:64,top:218,width:1152,height:294,columnWidths:[560,592],
 values:[
  ['What the model helps us ask','How this helps an experiment'],
  ['What could change the particle count?','Focus on release, particle loss and recovery.'],
  ['Could different processes explain the same data?','Choose measurements that tell them apart.'],
  ['Where do we need more biological data?','Add cell-injury measurements to assess particles from dying cells.'],
 ],
});
table.styleOptions={headerRow:false,bandedRows:false};
table.borders.assign({fill:'#D8DEE3',width:1,style:'solid'});
for(let r=0;r<4;r++)for(let c=0;c<2;c++) {
 const cell=table.getCell(r,c);cell.fill='#FFFFFF';
 cell.text.style={typeface:FONT,fontSize:27,bold:r===0,color:r===0?BLUE:INK,insets:12,verticalAlignment:'middle'};
}
addText(s5,'The current size data cannot tell us how much\ncomes from dying cells.',64,555,1148,86,32,true,BLUE);
addText(s5,'Early guidance from the current model and its assumptions',64,675,1100,26,19,false,MUTED);

change('sh/u1kbu1ov',titles[5],{fontSize:43,bold:true,color:INK});
change('sh/yhg7epsj','Each measurement helps separate possible explanations for the EV response',{fontSize:26,color:MUTED});
change('sh/87ipkzal','Check particle recovery');
change('sh/98rqt4r6','Measure how collection and processing affect the count.');
change('sh/ml07i9sv','Collect more time points');
change('sh/7m98ru9g','Add matched controls and more samples between 1 and 3 hours.');
change('sh/oryp8fah','Measure cell and EV biology');
change('sh/pc76hkr2','Pair size profiles with EV markers, calcium and cell injury.');
change('sh/h4bupgn6','Test new pulse conditions');
change('sh/w32twb6l','Check the model against data it has not seen.');
change('sh/v2tcn650','Use the model to compare explanations\nand plan the next measurements.',{fontSize:34,bold:true,color:INK});

// Refresh numbering after the user's slide additions and removals.
['sh/k3yl0zql','sh/ozy1ofad','sh/d0jax03i','sh/0ba143al','sh/cf2tcr61','sh/zi98nu94'].forEach((id,i)=>change(id,String(i+1)));

const scripts=[
 'We want to understand how short electric pulses change EV release, and how a model can help us plan the next experiment. A particle count on its own can have several explanations. Cells may release more vesicles, repair their membranes, or shed particles as they become injured. Some particles may also disappear before we collect the sample. We built a model that brings these processes together and compares their combined effect with laboratory measurements. The goal is to explore possible explanations for the response and identify which additional measurements would be most useful. I will show the biology we included, how the model compares with our current data, and how it can guide the next study.',
 'This overview shows the biology we brought into the model. The electric pulse affects the cell membrane and internal compartments. That can change calcium signals, cell stress and repair. Those changes can then affect how vesicles form and leave the cell, as well as how much unwanted material accompanies them. Cell type and cell condition can influence the whole response. Bringing these steps together lets us ask how a change in one part could affect what we eventually collect. The diagram also includes possible changes in vesicle contents and activity. Those parts remain early estimates rather than experimentally confirmed outcomes in this study.',
 'We organize the model around the same questions an experimentalist would ask. What pulse did the cells receive? How did the cells respond? What did they release? What reached the final sample? The distinction on the right is especially useful. Particles can leave a cell and then disappear from the medium before collection. Collection and measurement can change the count again. That means a change in the assay result does not automatically tell us how secretion changed. The model gives us a way to explore these steps separately, compare their combined effect with the data, and decide which part needs a direct measurement.',
 'Here, blue shows the measured size profiles and pink shows the model after we adjusted it using the available data. The example follows the three-pulse condition at thirty minutes, one hour and three hours. The model reproduces the main group of particles around one hundred to one hundred sixty nanometers. It also leaves a clear gap at the one-hour peak. That gap is useful because it shows where our current explanation is incomplete and where extra measurements could help. The full dataset contains twenty-six repeated size measurements across three pulse settings and three collection times. Their biological independence still needs confirmation. The forty-kilovolt label comes from the source file, and the delivered field still needs confirmation. This comparison describes the data used to adjust the model. We still need new conditions to test its predictions.',
 'For an experimentalist, the useful outcome is help choosing the next measurement. We vary the model inputs to see which changes could alter the particle count, then ask whether the available data can tell those explanations apart. In this dataset, release, disappearance from the medium, and particle recovery all need attention. The size measurements also cannot reliably tell us how much of the signal comes from dying cells. That gives us a reason to collect cell-injury measurements alongside particle counts. More generally, this approach helps focus experiments on processes that may change the time course and on gaps that make the results difficult to interpret. These are early priorities to test. Their order depends on the model assumptions and does not yet establish which biological process is truly most important.',
 'The next study can make the model much more useful by adding measurements that separate the possible explanations. First, measure how collection and processing affect particle recovery, and confirm the pulse settings and culture volume. Second, add matched controls and more collection times between one and three hours. That will help us examine whether cells keep releasing particles or whether particles disappear from the medium. Third, pair size profiles with EV markers, calcium responses and cell-injury measurements, so the particle signal has a clearer biological meaning. Finally, use additional pulse conditions to test predictions on new data. The value of the model is that it gives us a structured way to compare explanations and choose the next measurements.',
];
const sourceNotes=[
 `${PAPER}/main.tex, Author summary and Introduction`,
 `${PAPER}/main.tex, Materials and methods and module map. The overview illustration is retained from the user-edited deck.`,
 `${PAPER}/main.tex, Materials and methods, size-resolved observation pipeline and module map`,
 `${PAPER}/main.tex, In vitro particle-size measurements and Preliminary calibration; manuscript Figure 9`,
 `${PAPER}/main.tex, Global sensitivity analysis and Local parameter constraint and priorities for further measurements`,
 `${PAPER}/main.tex, measurement-priority discussion and Preliminary calibration`,
];
const cues=[
 'Introduce extracellular vesicles as EVs. Emphasize what the model can provide to an experimental study.',
 'Walk through the broad biological stages. The original overview contains detailed labels that can stay as reference.',
 'Emphasize the difference between what cells release and what reaches the assay.',
 'Point first to the shared size peak, then to the one-hour gap. Error bars show sample standard deviation. Size alone does not establish EV identity or origin.',
 'Lead with the experimental benefit. Particle recovery was provisionally fixed at one in fitting, with uncertainty absorbed into source scales. Recovery itself was not independently estimated. The original Gii diagnostic is intentionally replaced by its experimental implications.',
 'Connect each measurement with the explanation it would help test. This is an early measurement plan, not a statistically powered study design.',
];
for(let i=0;i<6;i++)pres.slides.items[i].speakerNotes.text=`SPEAKER SCRIPT\n\n${scripts[i]}\n\nPRESENTER CUE\n${cues[i]}\n\nSOURCE\n${sourceNotes[i]}`;
await fs.writeFile(OUT+'/electroexocytosis_experimentalist_speaker_script.txt',
 'ELECTRO-EXOCYTOSIS: SCRIPT FOR EXPERIMENTALISTS\n\n'+scripts.map((s,i)=>`SLIDE ${i+1}: ${titles[i]}\n\n${s}\n\nPresenter cue: ${cues[i]}\n\nSource: ${sourceNotes[i]}\n`).join('\n\n'));
await fs.writeFile(BUILD+'/revised_content.json',JSON.stringify(titles.map((t,i)=>({title:t,script:scripts[i]})),null,2));
for(let i=0;i<6;i++) {
 const b=await pres.export({slide:pres.slides.items[i],format:'png',scale:1.25});
 await fs.writeFile(BUILD+`/revised-${i+1}.png`,new Uint8Array(await b.arrayBuffer()));
}
const inspect=await pres.inspect({select:'$.slides[*]',include:['id','index','title','elements.id','elements.text','elements.position','notes.text']});
await fs.writeFile(BUILD+'/revised_inspect.ndjson',inspect.ndjson);
const candidate=BUILD+'/candidate_v2.pptx';
await(await PresentationFile.exportPptx(pres)).save(candidate);
const hash=createHash('sha256').update(await fs.readFile(SOURCE)).digest('hex');
const result=await finalizePresentation({
 workspaceDir:ROOT,candidatePath:candidate,finalPath:OUT+'/electroexocytosis_for_experimentalists_v2.pptx',
 explicitTotalSlideCount:6,requiredNativeTableOwnerSlides:[5],requiredNativeChartOwnerSlides:[],
 pythonExecutable:PYTHON,integrityValidatorPath:SKILL+'/container_tools/inspect_presentation_package_integrity.py',
 layoutValidatorPath:SKILL+'/container_tools/inspect_presentation_layout_geometry.py',
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit','--require-native-table-slide','5'],
 fontPolicy:{basis:'reference',families:[FONT],referencePath:SOURCE,referenceSha256:hash},
 verifyArtifactToolImport:true,receiptPath:BUILD+'/validation_v2.json',
});
console.log(JSON.stringify({finalPath:result.finalPath,slides:6,sourceSha256:hash}));
