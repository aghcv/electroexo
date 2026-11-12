import fs from 'node:fs/promises';
import { FileBlob, PresentationFile } from '@oai/artifact-tool';
const ROOT='/Users/aghorban/code/electroexo';
const BUILD=ROOT+'/tmp/plos_plain_language';
const p=await PresentationFile.importPptx(await FileBlob.load(ROOT+'/output/presentations/electroexocytosis_seven_plain_slides_v3.pptx'));
const r=await p.inspect({select:'$.slides[*]',include:['id','index','layout','title','elements.id','elements.kind','elements.text','elements.alt','elements.position','elements.style','notes.text']});
await fs.writeFile(BUILD+'/source_inspect.ndjson',r.ndjson);
console.log(r.ndjson);
console.log('MASTERS',p.masters.items.length);
for(let i=0;i<p.slides.items.length;i++) {
 const b=await p.export({slide:p.slides.items[i],format:'png',scale:1.25});
 await fs.writeFile(BUILD+`/source-${i+1}.png`,new Uint8Array(await b.arrayBuffer()));
}
