import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { Presentation, PresentationFile } from '@oai/artifact-tool';
import sharp from 'sharp';

const workspaceDir = '/Users/mac_seo/Research_local/7DTonCOSMOS';
const taskDir = path.join(workspaceDir, 'presentation/2026-10-02');
const TMP_DIR = path.join(taskDir, '.build');
const outputDir = path.join(taskDir, 'deliverables');
const FINAL_PPTX = path.join(outputDir, '7DT_COSMOS1_matching_presentation_clean.pptx');
const SKILL_DIR = '/Users/mac_seo/.codex/plugins/cache/openai-primary-runtime/presentations/26.909.12148/skills/presentations';
const RUNTIME_PYTHON = '/Users/mac_seo/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
const { finalizePresentation } = await import(pathToFileURL(path.join(SKILL_DIR, 'container_tools/artifact_tool_utils.mjs')).href);

await fs.mkdir(TMP_DIR, { recursive: true });
await fs.mkdir(outputDir, { recursive: true });

const W = 1280, H = 720;
const C = {
  ink: '#152334', muted: '#526272', blue: '#1765A3', red: '#C33C35',
  pale: '#E8F1F8', line: '#D7DEE5', white: '#FFFFFF', paper: '#FBFCFE'
};
const FONT = 'Arial';
const p = Presentation.create({ slideSize: { width: W, height: H } });
const figs = path.join(workspaceDir, 'output/COSMOS1/singlebandphot/matching/figures');
const img = {};
for (const name of ['singleband_input_summary.png','7dt_matching_summary.png','finalcat_radius_diagnostics.png','master_population_summary.png']) {
  img[name] = new Uint8Array(await fs.readFile(path.join(figs,name)));
}

function textbox(slide, text, x, y, w, h, opts={}) {
  const shape = slide.shapes.add({
    geometry: 'textbox', name: opts.name || text.slice(0,35),
    position: { left:x, top:y, width:w, height:h },
    fill: 'none', line: { fill:'none', width:0 }
  });
  shape.text = text;
  shape.text.style = {
    typeface: FONT, fontSize: opts.size || 22, bold: !!opts.bold,
    color: opts.color || C.ink, alignment: opts.align || 'left',
    verticalAlignment: opts.valign || 'middle', autoFit: 'none',
    wrap: true
  };
  return shape;
}
function line(slide,x1,y1,x2,y2,color=C.line,width=1) {
  return slide.shapes.add({geometry:'line',position:{left:x1,top:y1,width:x2-x1,height:y2-y1},fill:'none',line:{style:'solid',fill:color,width}});
}
function base(title, n, subtitle='') {
  const slide = p.slides.add();
  slide.background.fill = C.paper;
  textbox(slide,title,56,28,1168,59,{size:34,bold:true,name:`slide-${n}-title`});
  if (subtitle) textbox(slide,subtitle,58,91,1160,34,{size:19,color:C.muted,name:`slide-${n}-subtitle`});
  line(slide,58,119,1222,119,C.line,1);
  textbox(slide,`7DT · COSMOS1`,58,680,300,24,{size:14,color:C.muted});
  textbox(slide,String(n),1180,680,40,24,{size:14,color:C.muted,align:'right'});
  return slide;
}
async function cropPanel(slide, key, pos, crop, alt) {
  const meta=await sharp(img[key]).metadata();
  const left=Math.round(meta.width*crop.left), top=Math.round(meta.height*crop.top);
  const right=Math.round(meta.width*crop.right), bottom=Math.round(meta.height*crop.bottom);
  const panel=new Uint8Array(await sharp(img[key]).extract({left,top,width:meta.width-left-right,height:meta.height-top-bottom}).png().toBuffer());
  slide.images.add({blob:panel,contentType:'image/png',alt,fit:'contain',position:pos});
}
function nativeTable(slide, rows, pos, widths, opts={}) {
  const t = slide.tables.add({rows:rows.length,columns:rows[0].length,left:pos.x,top:pos.y,width:pos.w,height:pos.h,columnWidths:widths,values:rows});
  t.borders.assign({style:'solid',fill:C.line,width:0.65});
  for (let r=0;r<rows.length;r++) {
    t.rows[r].height = pos.h/rows.length;
    for (let c=0;c<rows[r].length;c++) {
      const cell=t.getCell(r,c);
      cell.fill = r===0 ? C.pale : (opts.highlightRow===r ? '#EEF6FC' : C.white);
      cell.text.style = {typeface:FONT,fontSize:opts.fontSize||17,bold:r===0 || opts.highlightRow===r,color:opts.highlightRow===r?C.blue:C.ink,alignment:c===0?'left':'right',verticalAlignment:'middle',autoFit:'none'};
    }
  }
  return t;
}

// Slide 1: integrated title and workflow.
{
  const s=base('7DT Single-band Matching in COSMOS1',1,'From 20 band catalogs to a 75,259-object union catalog');
  const steps=[
    ['20 band catalogs','78,882 selected detections'],
    ['7DT cross-band match','19,489 7DT masters'],
    ['FINALCAT sky match','2,977 shared objects'],
    ['Union catalog','75,259 total objects']
  ];
  steps.forEach((a,i)=>{
    const y=163+i*116;
    textbox(s,a[0],66,y,505,36,{size:25,bold:true,color:i===3?C.blue:C.ink});
    textbox(s,a[1],66,y+37,505,30,{size:20,color:C.muted});
    if(i<3) { line(s,85,y+78,85,y+107,C.blue,2); textbox(s,'↓',75,y+78,22,30,{size:22,color:C.blue,align:'center'}); }
  });
  await cropPanel(s,'singleband_input_summary.png',{left:625,top:150,width:565,height:225},{left:0,top:0.03,right:0.51,bottom:0.52},'Matching input count by 7DT band');
  await cropPanel(s,'singleband_input_summary.png',{left:625,top:391,width:565,height:225},{left:0,top:0.49,right:0.51,bottom:0.01},'Adopted position uncertainty by 7DT band');
  textbox(s,'Matching uses a position error for each band.',625,630,557,28,{size:18,color:C.muted});
  s.speakerNotes.textFrame.setText('Hello. Today I will show how I matched sources in the 7DT COSMOS1 data. This catalog is one step toward testing redshift estimates. I started with 20 separate band catalogs. I found and measured sources in each band. After quality cuts, I used 78,882 detections from the shared area. I first linked the 7DT detections across bands. Then I compared the new 7DT catalog with FINALCAT.\n\nSource: presentation/2026-10-02/matching_presentation_plan.md; output/COSMOS1/singlebandphot/matching/figures/singleband_input_summary.png.');
}

// Slide 2: method and requested editable schematic.
{
  const s=base('7DT Cross-band Matching',2,'Use position error to judge whether detections belong to one object');
  const labels=['New band source','3σ candidate search','One-to-one match','Weighted master update'];
  const xs=[60,365,670,975];
  labels.forEach((label,i)=>{
    textbox(s,String(i+1).padStart(2,'0'),xs[i],170,235,42,{size:23,bold:true,color:C.blue});
    textbox(s,label,xs[i],219,250,65,{size:23,bold:true});
    if(i<3){line(s,xs[i]+255,244,xs[i]+290,244,C.blue,2);textbox(s,'›',xs[i]+265,220,34,48,{size:33,color:C.blue,align:'center'});}
  });
  line(s,58,310,1222,310,C.line,1);
  textbox(s,'Band position error',62,344,330,34,{size:19,bold:true,color:C.muted});
  textbox(s,'σb = FWHMb / (2.355 × S/Nmin,b)',62,384,550,44,{size:27});
  textbox(s,'Pair score',672,344,280,34,{size:19,bold:true,color:C.muted});
  textbox(s,'D = d / √(σmaster² + σb²) < 3',672,384,550,44,{size:27,color:C.blue});
  textbox(s,'Updated master position',62,481,440,34,{size:19,bold:true,color:C.muted});
  textbox(s,'xmaster = Σ(xb / σb²) / Σ(1 / σb²)',62,521,750,48,{size:27});
  textbox(s,'High-precision positions get more weight.',62,590,1010,34,{size:21,color:C.muted});
  s.speakerNotes.textFrame.setText('The images share the same pixel grid. I processed the bands in wavelength order, starting with m400. A source position is less certain when the image is less sharp or the signal is weaker. So I gave each band its own position error. For each new source, I compared its distance from the current master with their combined error. I kept pairs within three times that error. When a source joined a master, I updated the master position. Better positions received more weight. This made one position for each linked object.\n\nSource: presentation/2026-10-02/matching_presentation_plan.md.');
}

// Slide 3: internal consistency.
{
  const s=base('7DT Matching Validation',3,'Accepted members agree with their final master position');
  textbox(s,'99.85%',60,150,450,85,{size:70,bold:true,color:C.blue});
  textbox(s,'within 3σ of the final master',64,233,680,37,{size:24,bold:true});
  textbox(s,'117 suspect members of 78,882 detections',770,180,460,42,{size:21,color:C.muted,align:'right'});
  await cropPanel(s,'7dt_matching_summary.png',{left:56,top:293,width:565,height:330},{left:0.666,top:0.012,right:0,bottom:0.508},'Normalized separation for accepted assignments');
  await cropPanel(s,'7dt_matching_summary.png',{left:655,top:293,width:565,height:330},{left:0,top:0.54,right:0.666,bottom:0.005},'Residual distance from final weighted master position');
  textbox(s,'19,489 masters · 1,157 detected in all 20 bands · R = dfinal / √(σfinal² + σb²)',59,631,1165,34,{size:17,color:C.muted});
  s.speakerNotes.textFrame.setText('Next, I checked the final 7DT masters. For each member, I measured its distance from the final master position. I divided this distance by the expected position error. Most members are close to their final master. In total, 99.85 percent are within three times the error. Only 117 of 78,882 members are outside this range. The final catalog has 19,489 7DT masters. This check shows good position agreement within the catalog. It does not prove that every match is correct.\n\nAdditional results: median normalized residual 0.364; 95th percentile 1.647; 99th percentile 2.333. The near-zero spike in assignment separation contains sources used to initialize masters. Source: presentation/2026-10-02/matching_presentation_plan.md; output/COSMOS1/singlebandphot/matching/figures/7dt_matching_summary.png.');
}

// Slide 4: radius evidence and editable table.
{
  const s=base('FINALCAT Match Radius',4,'The random-shift test supports a 0.7″ sky-match radius');
  const xs=[46,442,838];
  for(let i=0;i<3;i++) await cropPanel(s,'finalcat_radius_diagnostics.png',{left:xs[i],top:157,width:390,height:304},{left:i/3,top:0,right:(2-i)/3,bottom:0},['Matches by radius','Random association share','Separation of selected pairs'][i]);
  nativeTable(s,[
    ['Radius','Actual matches','Shifted random','Random / actual'],
    ['0.5″','2,682','55.8 ± 5.2','2.08%'],
    ['0.7″','2,977','107.9 ± 7.2','3.62%'],
    ['1.0″','3,212','211.4 ± 16.9','6.58%']
  ],{x:104,y:494,w:1072,h:154},[235,245,290,302],{highlightRow:2,fontSize:18});
  textbox(s,'Eight 10″ or 20″ shifts estimate chance associations.',104,653,1055,28,{size:17,color:C.muted});
  s.speakerNotes.textFrame.setText('I then matched the 7DT masters to FINALCAT using sky positions. I tested three maximum distances. To estimate chance matches, I shifted the FINALCAT positions and ran the match again. I used eight shifted samples. At 0.5 arcseconds, there were 2,682 real matches. At 0.7 arcseconds, there were 2,977. The share of chance matches was about 3.6 percent. At 1.0 arcseconds, the number of real matches rose to 3,212. But the share of chance matches rose faster, to about 6.6 percent. I therefore chose 0.7 arcseconds. The median distance of the selected pairs is 0.205 arcseconds.\n\nAdditional result: 95th percentile separation 0.578 arcseconds. Source: presentation/2026-10-02/matching_presentation_plan.md; output/COSMOS1/singlebandphot/matching/figures/finalcat_radius_diagnostics.png.');
}

// Slide 5: final populations and editable table.
{
  const s=base('Final Union Catalog',5,'2,977 shared objects in a 75,259-object catalog');
  await cropPanel(s,'master_population_summary.png',{left:54,top:148,width:555,height:330},{left:0,top:0,right:0.676,bottom:0},'Sky positions of matched and unmatched sources');
  await cropPanel(s,'master_population_summary.png',{left:666,top:148,width:555,height:330},{left:0.342,top:0,right:0.334,bottom:0},'Counts of matched, 7DT-only, and FINALCAT-only sources');
  nativeTable(s,[
    ['Master class','Objects'],
    ['7DT + FINALCAT','2,977'],
    ['7DT only','16,512'],
    ['FINALCAT only','55,770'],
    ['Total union catalog','75,259']
  ],{x:238,y:493,w:804,h:147},[550,254],{highlightRow:4,fontSize:17});
  s.speakerNotes.textFrame.setText('The final result has 2,977 objects seen in both 7DT and FINALCAT. It also has 16,512 objects found only in 7DT, and 55,770 found only in FINALCAT. Together, they form a union catalog of 75,259 objects. The main points are simple. The 7DT matches use the position error of each band. The internal position check is strong. The FINALCAT radius has a chance-match test behind it. This catalog is ready for the next photo-z analysis. Thank you.\n\nSource: presentation/2026-10-02/matching_presentation_plan.md; output/COSMOS1/singlebandphot/matching/figures/master_population_summary.png.');
}

// Save a private draft and run required checks before delivery.
const stagingDir=path.join(taskDir,'.codex-finalizer');
await fs.mkdir(stagingDir,{recursive:true});
const candidatePath=path.join(stagingDir,'candidate.pptx');
await (await PresentationFile.exportPptx(p)).save(candidatePath);
const result=await finalizePresentation({
  explicitTotalSlideCount:5,
  requiredNativeTableOwnerSlides:[4,5],
  requiredNativeChartOwnerSlides:[],
  workspaceDir,candidatePath,finalPath:FINAL_PPTX,
  pythonExecutable:RUNTIME_PYTHON,
  integrityValidatorPath:path.join(SKILL_DIR,'container_tools/inspect_presentation_package_integrity.py'),
  layoutValidatorPath:path.join(SKILL_DIR,'container_tools/inspect_presentation_layout_geometry.py'),
  layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit','--require-native-table-slide','4','--require-native-table-slide','5'],
  fontPolicy:{basis:'design',families:[FONT]},
  verifyArtifactToolImport:true,
  receiptPath:path.join(stagingDir,'matching_deck_clean.validation.json')
});
console.log(JSON.stringify({final:FINAL_PPTX,result},null,2));
