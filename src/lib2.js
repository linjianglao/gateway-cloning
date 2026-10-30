/* ============================ att detection & recombination ============================ */
const M1 = "TTTGTACAAAAAAGC";   // att1 core, top strand
const M2 = "GCTTTCTTGTACAAA";   // att2 core, top strand (rev-comp orientation)
const B1L="ACAAG", B1R="AGGCT", B2L="ACCCA", B2R="GTGGT";

function siteType(seq, i, which){
  const L=seq.slice(Math.max(0,i-5), i), R=seq.slice(i+15, i+20);
  const bl = which===1 ? L===B1L : L===B2L;
  const br = which===1 ? R===B1R : R===B2R;
  if(bl&&br) return "B";
  if(!bl&&br) return which===1?"L":"R";
  if(bl&&!br) return which===1?"R":"L";
  return "P";
}
function rotate(rec, off){
  const n=rec.seq.length;
  const seq=rec.seq.slice(off)+rec.seq.slice(0,off);
  const feats=(rec.features||[]).map(f=>{
    let s=f.start-1-off, e=f.end-1-off;
    if(s<0)s+=n; if(e<0)e+=n;
    if(e<s) return null;                       // feature straddles the new origin: drop it
    return {type:f.type,label:f.label,start:s+1,end:e+1,strand:f.strand};
  }).filter(Boolean);
  return Object.assign({},rec,{seq,features:feats});
}
/* Locate the two att cores; rotate a circular vector so site 1 precedes site 2. */
function findSites(rec){
  let r=rec;
  let i1=r.seq.indexOf(M1), i2=r.seq.indexOf(M2);
  if(i1<0||i2<0){
    const rcRec=Object.assign({},rec,{seq:rc(rec.seq),features:(rec.features||[]).map(f=>({
      type:f.type,label:f.label,start:rec.seq.length-f.end+1,end:rec.seq.length-f.start+1,strand:-f.strand}))});
    const j1=rcRec.seq.indexOf(M1), j2=rcRec.seq.indexOf(M2);
    if(j1>=0&&j2>=0){ r=rcRec; i1=j1; i2=j2; r.flipped=true; }
    else return {error:(i1<0&&i2<0)?"No Gateway att core sequence found in this sequence."
      :(i1<0?"Found an att2 core but no att1 core.":"Found an att1 core but no att2 core.")};
  }
  if(i1>i2){
    if(r.topology!=="circular") return {error:"att sites are in the wrong order for a linear sequence — att1 must precede att2."};
    r=rotate(r, Math.max(0,i1-60));
    i1=r.seq.indexOf(M1); i2=r.seq.indexOf(M2);
    if(i1<0||i2<0||i1>i2) return {error:"Could not orient the att cassette on this circular sequence."};
  }
  const t1=siteType(r.seq,i1,1), t2=siteType(r.seq,i2,2);
  const role = (t1==="P"&&t2==="P")?"donor" : (t1==="L"&&t2==="L")?"entry"
             : (t1==="R"&&t2==="R")?"dest"  : (t1==="B"&&t2==="B")?"expression" : "mixed";
  return {rec:r, i1, i2, t1, t2, role, cassetteBp:(i2+15)-(i1)};
}
const SITE_ARM = {B:5, L:80, R:110, P:110};   // bp of each site lying outside the 15-bp core, upstream side
/* Swap the region between the two cores for a new insert. Works for BP (donor) and LR (destination). */
function recombine(rec, sites, mid, insertLabel){
  const s=sites.rec.seq, i1=sites.i1, i2=sites.i2;
  const left=s.slice(0,i1), right=s.slice(i2+15);
  const core=M1+B1R+mid+B2L+M2;
  const seq=left+core+right;
  const oldLen=(i2+15)-i1, newLen=core.length, delta=newLen-oldLen;
  const feats=[];
  (sites.rec.features||[]).forEach(f=>{
    if(f.end<=i1) feats.push(Object.assign({},f));
    else if(f.start-1>=i2+15) feats.push(Object.assign({},f,{start:f.start+delta,end:f.end+delta}));
    // features overlapping the cassette are the ccdB/Cm markers being replaced — dropped
  });
  // annotate the two new att sites and the insert (1-based inclusive)
  const p1 = sites.t1==="B" ? "attL1" : (sites.t1==="P" ? "attL1" : (sites.t1==="R" ? "attB1" : "attB1"));
  const p2 = sites.t1==="P" ? "attL2" : "attB2";
  const name1 = sites.role==="donor" ? "attL1" : sites.role==="dest" ? "attB1" : (sites.t1==="P"?"attL1":"attB1");
  const name2 = sites.role==="donor" ? "attL2" : sites.role==="dest" ? "attB2" : (sites.t2==="P"?"attL2":"attB2");
  const arm1 = sites.role==="donor" ? ATT.L1.indexOf(M1) : 5;
  const arm2 = sites.role==="donor" ? (ATT.L2.length-ATT.L2.indexOf(M2)-15) : 5;
  const s1s=Math.max(1,i1-arm1+1), s1e=i1+20;
  const insS=i1+20+1, insE=i1+20+mid.length;
  const s2s=insE+1, s2e=Math.min(seq.length, insE+5+15+arm2);
  feats.push({type:"misc_recomb",label:name1,start:s1s,end:s1e,strand:1});
  feats.push({type:"CDS",label:insertLabel||"insert",start:insS,end:insE,strand:1});
  feats.push({type:"misc_recomb",label:name2,start:s2s,end:s2e,strand:1});
  feats.sort((a,b)=>a.start-b.start);
  // reciprocal by-product (linear fragment carrying ccdB)
  const byL = sites.role==="donor" ? B1L : ATT.L1.slice(0, ATT.L1.indexOf(M1));
  const byR = sites.role==="donor" ? B2R : ATT.L2.slice(ATT.L2.indexOf(M2)+15);
  const byproduct = byL+M1+s.slice(i1+15,i2)+M2+byR;
  return {seq, features:feats, delta, byproduct,
          site1:{start:s1s,end:s1e,name:name1}, insert:{start:insS,end:insE},
          site2:{start:s2s,end:s2e,name:name2}};
}

/* ============================ primer design (Won Do Heo protocol) ============================ */
const P_TAIL_B1 = "GGGGACAAGTTTGTACAAAAAAGCAGGCTTC";      // GGGG + attB1 + TC
const P_TAIL_B2 = "GGGGACCACTTTGTACAAGAAAGCTGGGT";        // GGGG + attB2 (rev)
const ADAPT_B1  = "GGGGACAAGTTTGTACAAAAAAGCAGGCT";
const ADAPT_B2  = "GGGGACCACTTTGTACAAGAAAGCTGGGT";
const SHORT_B1  = "AAAAAGCAGGCTTC";                        // 12-base adapter-PCR tail
const SHORT_B2  = "AGAAAGCTGGGTG";
const SD_KOZAK  = "GAAGGAGATAGAACCATG";                    // Shine-Dalgarno + Kozak + ATG
const KOZAK     = "ACCATG";
const TEV_FWD   = "GAAAACCTGTATTTTCAGGGCATG";              // ENLYFQG + ATG
const FORMS = {
  A:{side:"B1",add:SD_KOZAK,short:SHORT_B1,
     title:"A — Shine–Dalgarno + Kozak + ATG",
     blurb:"Carries both ribosome-recognition sequences and the start codon, so one clone expresses in E. coli and in eukaryotes. The ATG stays in frame with attB1, keeping an N-terminal fusion possible later."},
  B:{side:"B1",add:KOZAK,short:SHORT_B1,
     title:"B — Kozak + ATG",
     blurb:"Eukaryotic ribosome recognition only. The ATG stays in frame with attB1."},
  C:{side:"B1",add:"",short:SHORT_B1,
     title:"C — N-terminal fusion",
     blurb:"No rrs and no start codon — both come from the destination vector. The gene-specific sequence must continue the attB1 reading frame; including your own ATG is optional."},
  K:{side:"B2",add:"TTA",short:SHORT_B2,
     title:"K — stop codon added",
     blurb:"Adds a stop codon (TAA on the coding strand) so translation terminates before attB2. CTA and TCA work here too."},
  L:{side:"B2",add:"",extraBase:"G",short:SHORT_B2,
     title:"L — C-terminal fusion",
     blurb:"No stop codon, plus one extra base to complete the codon that runs into attB2, so the reading frame carries on into the vector's C-terminal tag."}
};
const MATRIX = {
  ecoli:      {native:["A","K"], nterm:["C","K"], cterm:["A","L"], both:["A","K"]},
  yeast:      {native:["B","K"], nterm:["C","K"], cterm:["B","L"], both:["B","K"]},
  baculo:     {native:["B","K"], nterm:["C","K"], cterm:["B","L"], both:["B","K"]},
  mammalian:  {native:["B","K"], nterm:["C","K"], cterm:["B","L"], both:["B","K"]},
  ybm:        {native:["B","K"], nterm:["C","K"], cterm:["B","L"], both:["B","K"]},
  all:        {native:["A","K"], nterm:["C","K"], cterm:["A","L"], both:["A","K"]}
};
const ORGANISMS={ecoli:"E. coli only",yeast:"Yeast only",baculo:"Baculovirus only",
  mammalian:"Mammalian only",ybm:"Yeast, baculovirus and mammalian",all:"E. coli, yeast, baculovirus and mammalian"};
const FORMATS={native:"Native protein",nterm:"N-terminal fusion",cterm:"C-terminal fusion",both:"Native and N-terminal fusion"};
const STOPS=["TAA","TAG","TGA"];

function pickGsp(seq, fromStart, minLen, maxLen){
  // choose 18-24 nt with the Tm closest to 60 C, preferring a G/C 3' end
  let best=null;
  for(let n=minLen;n<=maxLen;n++){
    const s=fromStart?seq.slice(0,n):seq.slice(seq.length-n);
    if(s.length<n) break;
    const t=tm(s), last=fromStart?s[s.length-1]:s[s.length-1];
    const score=Math.abs(t-60) - (("GC".indexOf(last)>=0)?0.6:0);
    if(!best||score<best.score) best={s,t,score,n};
  }
  return best||{s:seq.slice(0,minLen),t:0,n:minLen};
}
function designProtocolPrimers(cfg){
  const orf=clean(cfg.orf);
  const notes=[];
  const hasStart=orf.slice(0,3)==="ATG";
  const endStop=STOPS.indexOf(orf.slice(-3))>=0;
  const body=(hasStart?orf.slice(3):orf).slice(0, (hasStart?orf.slice(3):orf).length - (endStop?3:0));
  if(!hasStart) notes.push({k:"warn",t:"<b>No ATG at the 5′ end.</b> The primer supplies one; check that this really is your start codon."});
  if(body.length%3) notes.push({k:"bad",t:"<b>The coding sequence is not a multiple of three</b> ("+body.length+" nt after the ATG and stop codon are removed). Fix the ORF before ordering primers — every downstream reading frame depends on it."});
  const pair=MATRIX[cfg.organism][cfg.format];
  const f1=FORMS[pair[0]], f2=FORMS[pair[1]];

  // forward
  const fwdGsp=pickGsp(body,true,18,24);
  const includeAtg = pair[0]==="C" && cfg.cAtg;
  const fAdd = pair[0]==="C" ? (includeAtg?"ATG":"") : f1.add;
  const fTev = cfg.tev && pair[0]==="C" ? TEV_FWD : (cfg.tev ? TEV_FWD : "");
  const fCore = cfg.tev ? TEV_FWD : fAdd;
  const fwd = (cfg.protocol==="adapter" ? SHORT_B1 : P_TAIL_B1) + fCore + fwdGsp.s;
  // stop-codon check at the attB1 junction
  const junctionCodon = (fCore+fwdGsp.s).slice(0,3);
  if(STOPS.indexOf(junctionCodon)>=0)
    notes.push({k:"bad",t:"<b>The codon straight after attB1 is "+junctionCodon+" — a stop codon.</b> Shift the gene-specific sequence by one codon, or include the ATG, so TTA, TAG and TGA cannot form at this junction."});

  // reverse
  const revGspSrc = pair[1]==="K" ? body : body;
  const revPick=pickGsp(revGspSrc,false,18,24);
  const revGsp=rc(revPick.s);
  const rTail=(cfg.protocol==="adapter" ? SHORT_B2 : P_TAIL_B2 + (f2.extraBase||""));
  const rev = rTail + (f2.add||"") + revGsp;

  // top strand between attB1 and attB2
  const mid = "TC" + fCore + body + (pair[1]==="K" ? "TAA" : "C");
  const inFrame = ((ATT.B1.length + mid.length) % 3)===0;

  if(pair[1]==="L"){
    notes.push({k:body.length%3?"bad":"ok",
      t:"<b>C-terminal fusion.</b> The added "+(f2.extraBase||"G")+" completes the codon that runs into attB2 (read as CAC on the coding strand), so the frame carries into the vector tag. The gene-specific sequence must end exactly on a codon boundary."});
  } else {
    notes.push({k:"ok",t:"<b>Stop codon supplied by the reverse primer</b> (TTA on the primer = TAA on the coding strand). Nothing downstream of attB2 is translated."});
  }
  if(cfg.tev) notes.push({k:"ok",t:"<b>TEV site included.</b> Cleavage leaves a single glycine on the N-terminus; replace the ATG with GGC if you would rather have a Met→Gly substitution than an extra residue."});
  if(cfg.protocol==="adapter") notes.push({k:"ok",t:"<b>Adapter PCR.</b> Run a first PCR with these short-tailed gene primers, then a second with the universal adapter primers below. Cheaper and far less prone to primer dimers than full-length attB primers."});

  return {orf, body, hasStart, endStop, pair, f1, f2, fwd, rev, mid, inFrame, notes,
    fwdGsp:fwdGsp.s, revGsp, fwdTm:fwdGsp.t, revTm:revPick.t, fCore,
    fwdParts:[{s:(cfg.protocol==="adapter"?SHORT_B1:P_TAIL_B1),c:"sB"},{s:fCore,c:"sK"},{s:fwdGsp.s,c:"sO"}],
    revParts:[{s:rTail,c:"sB"},{s:(f2.add||""),c:"sK"},{s:revGsp,c:"sO"}]};
}
