import fs from 'node:fs/promises';
import {pathToFileURL} from 'node:url';
import {FileBlob,PresentationFile} from '@oai/artifact-tool';
const ROOT='/Users/aghorban/code/electroexo',BUILD=ROOT+'/tmp/eight_module_slides',OUT=ROOT+'/output/presentations';
const SKILL='/Users/aghorban/.codex/plugins/cache/openai-primary-runtime/presentations/26.1013.11914/skills/presentations';
const PYTHON='/Users/aghorban/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
const {finalizePresentation}=await import(pathToFileURL(SKILL+'/container_tools/artifact_tool_utils.mjs'));
const pres=await PresentationFile.importPptx(await FileBlob.load(OUT+'/electroexocytosis_with_eight_module_studies.pptx'));
pres.slides.keep([6,7,8,9,10]);
const snap=await pres.inspect({select:'$.slides[*]',include:['index','title','elements.id','elements.kind','elements.text','elements.position']});
for(const rec of snap.records) {
 const footer=rec.elements.find(e=>e.kind==='textbox'&&e.position?.left===1224&&e.position?.top===685);
 if(!footer)throw new Error('Missing study slide footer');
 pres.resolve(footer.id).text=String(rec.index+1);
}
const candidate=BUILD+'/candidate_five_studies_imported.pptx';
await(await PresentationFile.exportPptx(pres)).save(candidate);
const result=await finalizePresentation({
 workspaceDir:ROOT,candidatePath:candidate,finalPath:OUT+'/electroexocytosis_five_eight_module_studies.pptx',
 explicitTotalSlideCount:5,requiredNativeChartOwnerSlides:[1,2,3,4,5],
 pythonExecutable:PYTHON,integrityValidatorPath:SKILL+'/container_tools/inspect_presentation_package_integrity.py',
 layoutValidatorPath:SKILL+'/container_tools/inspect_presentation_layout_geometry.py',
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit'],
 fontPolicy:{basis:'design',families:['Helvetica Neue']},materializeLiteralChartWorkbooks:true,verifyArtifactToolImport:true,
 receiptPath:BUILD+'/validation_five_studies_imported.json',
});
console.log('Study deck: '+result.finalPath);
const manifest=JSON.parse(await fs.readFile(ROOT+'/results/eight_module_storyboards/manifest.json','utf8'));
const script=manifest.studies.map((s,i)=>`SLIDE ${i+1}: ${snap.records[i].title}\n\n${s.script}\n\nScenario details:\n${s.rows.map(r=>r.label+': '+r.detail).join('\n')}\n\n${s.caption}`).join('\n\n\n');
await fs.writeFile(OUT+'/electroexocytosis_eight_module_study_scripts.txt',script);
for(let i=0;i<5;i++) {
 const b=await pres.export({slide:pres.slides.items[i],format:'png',scale:1.5});
 await fs.writeFile(BUILD+`/final-study-${i+1}.png`,new Uint8Array(await b.arrayBuffer()));
}
