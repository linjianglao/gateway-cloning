/* ============================ state ============================ */
const state = {
  tab:"bp",
  bp:{ mode:"design", orf:SAMPLE_ORF, name:"EGFP", donor:"pDONR221",
       format:"native", organism:"all", protocol:"nonstop", tev:false, cAtg:false, product:"" },
  lr:{ mode:"carry", cassette:"", name:"EGFP", dest:"pcDNA-DEST40" },
  lib:{ kind:"all", app:"all", q:"", sel:"pDONR221", adding:false,
        draft:{name:"",addgene:"",kind:"dest",app:"mammalian",prom:"",tag:"none",tagpos:"none",
               bact:"Ampicillin",sel:"—",note:"",rec:null} },
  custom:[], history:[], caps:{db:"pending",dl:"pending"}, msg:null
};
let DB=null, DL=null;

/* ============================ vectors ============================ */
function allVectors(){ return VECTORS.concat(state.custom); }
function getVec(id){ return allVectors().find(v=>v.id===id) || null; }
function vecRec(id){
  const v=getVec(id);
  if(!v||!v.seq) return null;
  return {name:v.id, seq:v.seq, features:v.features||[], topology:v.topology||"circular"};
}

/* ============================ capabilities ============================ */
function boot(){
  const use = (window.claude && typeof window.claude.use==="function") ? window.claude.use.bind(window.claude) : null;
  if(!use){ state.caps.db="off"; state.caps.dl="off"; rerenderAll(); return; }
  use("db").then(d=>{
    DB=d; state.caps.db = d?"on":"off";
    if(d){ loadCustom(); loadHistory(); } else rerenderAll();
  }).catch(()=>{ state.caps.db="off"; rerenderAll(); });
  use("downloads").then(d=>{ DL=d; state.caps.dl=d?"on":"off"; rerenderAll(); })
    .catch(()=>{ state.caps.dl="off"; rerenderAll(); });
}
function loadCustom(){
  DB.collection("vectors").limit(200).get().then(snap=>{
    state.custom = snap.docs.map(d=>Object.assign({custom:true,id:d.id}, d.data()));
    rerenderAll();
  }).catch(()=>{});
}
function loadHistory(){
  DB.collection("runs").orderBy("createdAt","desc").limit(60).get().then(snap=>{
    state.history = snap.docs.map(d=>Object.assign({docId:d.id}, d.data()));
    render("hist");
  }).catch(()=>{});
}
function flash(t,kind){ state.msg={t,kind:kind||"ok"}; rerenderAll();
  setTimeout(()=>{ state.msg=null; rerenderAll(); }, 4200); }

/* ============================ reactions ============================ */
function midFromAttB(text){
  const p=clean(text);
  const i1=p.indexOf(ATT.B1), i2=p.indexOf(ATT.B2);
  if(i1<0) return {error:"No attB1 site found. The product must contain "+ATT.B1+"."};
  if(i2<0) return {error:"No attB2 site found. The product must contain "+ATT.B2+"."};
  if(i2<=i1+ATT.B1.length) return {error:"attB2 appears before attB1 — check the orientation of your product."};
  return {mid:p.slice(i1+ATT.B1.length,i2)};
}
function midFromAttL(text){
  const s=clean(text);
  const i1=s.indexOf(M1), i2=s.indexOf(M2);
  if(i1<0||i2<0) return {error:"Could not find both att cores. Paste an attL1…attL2 entry cassette, or upload the entry clone file."};
  if(i2<=i1+15) return {error:"att sites are out of order — att1 must come before att2."};
  return {mid:s.slice(i1+20, i2-5)};
}
function reaction(kind){
  const isBP = kind==="BP";
  const c = isBP?state.bp:state.lr;
  let mid=null, design=null, err=null;
  if(isBP){
    if(c.mode==="design"){ design=designProtocolPrimers(c); mid=design.mid; }
    else { const r=midFromAttB(c.product); if(r.error) err=r.error; else mid=r.mid; }
  } else {
    if(c.mode==="carry"){ const b=reaction("BP"); if(b.error) err="Finish the BP reaction first — "+b.error; else mid=b.mid; design=b.design; }
    else { const r=midFromAttL(c.cassette); if(r.error) err=r.error; else mid=r.mid; }
  }
  if(err) return {error:err};
  const vecId = isBP?c.donor:c.dest;
  const vec = getVec(vecId);
  const rec = vecRec(vecId);
  const name = isBP ? ("pENTR-"+(state.bp.name||"insert"))
                    : (String(vecId).replace(/\s*\(.*\)$/,"")+"-"+((state.lr.mode==="carry"?state.bp.name:c.name)||"insert"));
  const cassette = isBP ? (ATT.L1+mid+ATT.L2) : (ATT.B1+mid+ATT.B2);
  let full=null, sites=null, siteWarn=null;
  if(rec){
    sites=findSites(rec);
    if(sites.error) siteWarn=sites.error;
    else if(sites.role!==(isBP?"donor":"dest"))
      siteWarn="This vector reads as a"+(sites.role==="mixed"?" mixed-site":"n "+sites.role)+" vector ("+sites.t1+"/"+sites.t2+" sites), not a "+(isBP?"donor (attP × attP)":"destination (attR × attR)")+" vector.";
    else full=recombine(rec,sites,mid,(isBP?state.bp.name:(state.lr.mode==="carry"?state.bp.name:c.name))||"insert");
  }
  const plasmid = full ? full.seq.length
    : (vec ? vec.size - (vec.cass||0) + cassette.length : cassette.length);
  return {kind, mid, design, vec, rec, sites, siteWarn, full, cassette, plasmid, name, insertName:(isBP?state.bp.name:(state.lr.mode==="carry"?state.bp.name:c.name))};
}

/* ============================ files out ============================ */
function bundleFor(r){
  const base=slug(r.name);
  const seq = r.full ? r.full.seq : r.cassette;
  const feats = r.full ? r.full.features : (r.kind==="BP"
      ? [{type:"misc_recomb",label:"attL1",start:1,end:100,strand:1},
         {type:"CDS",label:r.insertName||"insert",start:101,end:100+r.mid.length,strand:1},
         {type:"misc_recomb",label:"attL2",start:101+r.mid.length,end:200+r.mid.length,strand:1}]
      : [{type:"misc_recomb",label:"attB1",start:1,end:25,strand:1},
         {type:"CDS",label:r.insertName||"insert",start:26,end:25+r.mid.length,strand:1},
         {type:"misc_recomb",label:"attB2",start:26+r.mid.length,end:50+r.mid.length,strand:1}]);
  const gb=toGenBank({name:base, desc:(r.kind==="BP"?"Gateway entry clone":"Gateway expression clone")+
    " — "+(r.insertName||"insert")+" in "+(r.vec?r.vec.id:"vector"), seq, features:feats,
    topology:r.full?"circular":"linear"});
  const files=[{name:base+".gb", text:gb},{name:base+".fasta", text:toFastaText(base,seq)}];
  if(r.design){
    files.push({name:base+"-primers.csv", text:
      "name,sequence,length,tm_gsp_c,form,note\n"+
      '"'+base+'_attB1_F","'+r.design.fwd+'",'+r.design.fwd.length+","+r.design.fwdTm.toFixed(1)+","+r.design.pair[0]+',"'+FORMS[r.design.pair[0]].title.replace(/"/g,"")+'"\n'+
      '"'+base+'_attB2_R","'+r.design.rev+'",'+r.design.rev.length+","+r.design.revTm.toFixed(1)+","+r.design.pair[1]+',"'+FORMS[r.design.pair[1]].title.replace(/"/g,"")+'"\n'+
      (state.bp.protocol==="adapter"
        ? '"adapter_attB1","'+ADAPT_B1+'",'+ADAPT_B1.length+',,"universal adapter primer"\n'+
          '"adapter_attB2","'+ADAPT_B2+'",'+ADAPT_B2.length+',,"universal adapter primer"\n' : "")+
      SEQ_PRIMERS.map(p=>'"'+p.name+'","'+p.seq+'",'+p.seq.length+',,"'+p.note+'"').join("\n")+"\n"});
  }
  files.push({name:"README.txt", text:
    "Gateway Bench export\n====================\n\n"+
    "Construct : "+r.name+"\n"+
    "Reaction  : "+r.kind+" ("+(r.kind==="BP"?"attB x attP -> attL + attR":"attL x attR -> attB + attP")+")\n"+
    "Vector    : "+(r.vec?r.vec.id:"(none)")+"\n"+
    "Insert    : "+(r.insertName||"insert")+", "+r.mid.length+" bp between the att cores\n"+
    "Length    : "+seq.length+" bp "+(r.full?"(complete recombined plasmid)":"(cassette only — no vector sequence was supplied)")+"\n"+
    "Generated : "+new Date().toISOString()+"\n\n"+
    (r.full?"":"To get the complete plasmid, add the vector's GenBank file on the Vector library tab and run the reaction again.\n")});
  return {base, files, seq, gb};
}
function saveFile(filename, data){
  if(!DL){ flash("Downloads are not available in this view.","bad"); return; }
  DL.save({filename, data}).then(()=>flash("Saved "+filename,"ok"))
    .catch(e=>{
      const c=e&&e.code;
      if(c==="declined") return;
      flash(c==="rate_limited"?"A download prompt is already open — try again in a moment."
        :c==="too_large"?"That file is too large for this download."
        :"Could not save the file ("+(c||"unknown")+").","bad");
    });
}
function downloadBundle(r){ const b=bundleFor(r); saveFile(b.base+".zip", zipStore(b.files)); }
function downloadGb(r){ const b=bundleFor(r); saveFile(b.base+".gb.txt", b.gb); }

/* ============================ history ============================ */
function saveRun(r){
  if(!DB){ flash("History needs the artifact's data store, which is not available in this view.","bad"); return; }
  const id=String(Date.now())+"-"+Math.random().toString(36).slice(2,7);
  const doc={kind:r.kind, name:r.name, insert:r.insertName||"insert", vector:r.vec?r.vec.id:"(none)",
    insertBp:r.mid.length, plasmidBp:r.plasmid, complete:!!r.full,
    mid:r.mid.slice(0,60000), createdAt:new Date().toISOString(),
    format:r.kind==="BP"?state.bp.format:"", organism:r.kind==="BP"?state.bp.organism:"",
    fwd:r.design?r.design.fwd:"", rev:r.design?r.design.rev:""};
  DB.collection("runs").doc(id).set(doc).then(()=>{
    state.history.unshift(Object.assign({docId:id},doc)); flash("Saved to history","ok");
  }).catch(e=>flash(e&&e.code==="invalid_argument"
    ? "You have read-only access to this page, so history cannot be written."
    : "Could not save to history.","bad"));
}
function deleteRun(docId){
  if(!DB) return;
  DB.collection("runs").doc(docId).delete().then(()=>{
    state.history=state.history.filter(h=>h.docId!==docId); render("hist");
  }).catch(()=>flash("Could not delete that entry.","bad"));
}
function restoreRun(h){
  if(h.kind==="BP"){ state.bp.mode="paste"; state.bp.product=ATT.B1+h.mid+ATT.B2; state.bp.name=h.insert; switchTab("bp"); }
  else { state.lr.mode="paste"; state.lr.cassette=ATT.L1+h.mid+ATT.L2; state.lr.name=h.insert; switchTab("lr"); }
}

/* ============================ vector add ============================ */
function saveCustomVector(){
  const d=state.lib.draft;
  if(!d.name){ flash("Give the vector a name first.","bad"); return; }
  const id=slug(d.name);
  const doc={name:d.name, addgene:d.addgene||"", kind:d.kind, app:d.app, size:d.rec?d.rec.seq.length:(d.size||0),
    cass:d.cass||0, bact:d.bact||"—", sel:d.sel||"—", prom:d.prom||"—", tag:d.tag||"none", tagpos:d.tagpos||"none",
    att:d.att||"—", src:d.addgene?("Addgene #"+d.addgene):"user upload", note:d.note||"",
    seq:d.rec?d.rec.seq:"", features:d.rec?d.rec.features.slice(0,300):[], topology:d.rec?d.rec.topology:"circular"};
  if(JSON.stringify(doc).length>240000){ flash("That sequence is too large to store (256 KB limit per vector).","bad"); return; }
  if(!DB){
    state.custom.push(Object.assign({custom:true,id,session:true},doc));
    state.lib.adding=false; state.lib.sel=id;
    flash("Added for this session only — the data store is unavailable, so it will not persist.","warn");
    rerenderAll(); return;
  }
  DB.collection("vectors").doc(id).set(doc).then(()=>{
    state.custom=state.custom.filter(v=>v.id!==id);
    state.custom.push(Object.assign({custom:true,id},doc));
    state.lib.adding=false; state.lib.sel=id;
    flash("Added "+d.name+" to the library","ok"); rerenderAll();
  }).catch(e=>flash(e&&e.code==="invalid_argument"
    ? "You have read-only access to this page, so the library cannot be changed."
    : "Could not save the vector.","bad"));
}
function deleteCustomVector(id){
  if(!DB){ state.custom=state.custom.filter(v=>v.id!==id); rerenderAll(); return; }
  DB.collection("vectors").doc(id).delete().then(()=>{
    state.custom=state.custom.filter(v=>v.id!==id);
    if(state.lib.sel===id) state.lib.sel="pDONR221";
    rerenderAll();
  }).catch(()=>flash("Could not remove that vector.","bad"));
}
function readFile(file, cb){
  const fr=new FileReader();
  fr.onload=()=>cb(null,String(fr.result));
  fr.onerror=()=>cb(new Error("read failed"));
  fr.readAsText(file);
}
function handleSeqFile(file, target){
  if(/\.dna$/i.test(file.name)){
    flash("SnapGene .dna files are binary — open it in SnapGene and export as GenBank (.gb) first.","bad"); return;
  }
  readFile(file,(err,text)=>{
    if(err){ flash("Could not read that file.","bad"); return; }
    const recs=parseSeqFile(file.name,text);
    if(!recs.length){ flash("No sequence found in "+file.name+" — expected FASTA, GenBank or plain sequence.","bad"); return; }
    const rec=recs[0];
    if(target==="orf"){
      state.bp.orf=rec.seq; state.bp.name=rec.name||state.bp.name; state.bp.mode="design";
      flash("Loaded "+rec.name+" — "+bp(rec.seq.length),"ok");
    } else if(target==="entry"){
      state.lr.mode="paste"; state.lr.cassette=rec.seq; state.lr.name=rec.name||state.lr.name;
      flash("Loaded "+rec.name+" — "+bp(rec.seq.length),"ok");
    } else {
      const d=state.lib.draft;
      d.rec=rec; d.name=d.name||rec.name;
      const s=findSites(rec);
      if(s.error){ d.att="—"; flash("Loaded "+bp(rec.seq.length)+", but no Gateway att sites were found: "+s.error,"warn"); }
      else {
        d.rec=s.rec; d.kind = s.role==="donor"?"donor":s.role==="entry"?"entry":s.role==="dest"?"dest":d.kind;
        d.att="att"+s.t1+"1 / att"+s.t2+"2"; d.cass=s.cassetteBp;
        const ccdb=/ATGCAGTTTAAGGTTTACACCTATAAAAGA/.test(rec.seq);
        const labels=(rec.features||[]).map(f=>String(f.label||"").toLowerCase());
        const find=(re)=>{const f=(rec.features||[]).find(x=>re.test(String(x.label||"")));return f?f.label:null;};
        d.bact = find(/^(AmpR|bla|Amp)/i)?"Ampicillin" : find(/KanR|nptII|neo/i)?"Kanamycin"
               : find(/SmR|SpecR|aadA/i)?"Spectinomycin" : find(/GmR|aacC1/i)?"Gentamicin" : d.bact;
        d.prom = find(/CMV|EF-1|EF1|PGK|T7|SV40|35S|TRE|UBC|ADH1|polyhedrin/i) || d.prom;
        d.note = (d.note||"") + (ccdb?" ccdB cassette detected in the uploaded sequence.":"");
        flash("Read "+bp(rec.seq.length)+" · "+(rec.features||[]).length+" features · "+
              (s.role==="mixed"?"mixed att sites":s.role+" vector")+" ("+d.att+")"+(ccdb?" · ccdB present":""),"ok");
      }
      state.lib.adding=true;
    }
    rerenderAll();
  });
}

/* ============================ shared view bits ============================ */
function copyBtn(id,label){return '<button class="btn ghost sm" data-copy="'+id+'">'+(label||"Copy")+"</button>";}
const APPNAMES={universal:"Universal",multisite:"MultiSite",bacterial:"E. coli",mammalian:"Mammalian",
  lentiviral:"Lentiviral",retroviral:"Retroviral",adenoviral:"Adenoviral",insect:"Insect",yeast:"Yeast",
  plant:"Plant",custom:"Custom"};
function vecOptions(kind, sel){
  const groups={};
  allVectors().filter(v=>v.kind===kind).forEach(v=>{(groups[v.app]=groups[v.app]||[]).push(v);});
  return Object.keys(groups).map(g=>'<optgroup label="'+esc(APPNAMES[g]||g)+'">'+
    groups[g].map(v=>'<option value="'+esc(v.id)+'"'+(v.id===sel?" selected":"")+">"+esc(v.id)+(v.seq?" ✓":"")+"</option>").join("")+
    "</optgroup>").join("");
}
function msgBar(){
  if(!state.msg) return "";
  return '<div class="panel" style="margin-bottom:16px"><div class="panel-bd"><div class="note '+
    (state.msg.kind==="bad"?"bad":state.msg.kind==="warn"?"warn":"ok")+'">'+state.msg.t+"</div></div></div>";
}
function productPanel(r){
  const isBP=r.kind==="BP";
  const cls=isBP?"L":"B";
  const armBp=isBP?100:25;
  const feats=[feat("att"+cls+"1",0,armBp,cls),
               feat(r.insertName||"insert",armBp,armBp+r.mid.length,"O",true),
               feat("att"+cls+"2",armBp+r.mid.length,2*armBp+r.mid.length,cls)];
  let circ;
  if(r.full){
    const chosen=r.full.features.filter(f=>(f.end-f.start)>150 || /^att/.test(f.label))
      .sort((a,b)=>(b.end-b.start)-(a.end-a.start)).slice(0,8).sort((a,b)=>a.start-b.start);
    circ=circularMap(r.name, r.full.seq.length, chosen.map(f=>
      feat(f.label, f.start, f.end, /^att/.test(f.label)?featColorFor(f.label):(f.label===r.insertName?"O":"N"), (f.end-f.start)>200)));
  } else {
    const v=r.vec, b=armBp, m=r.mid.length;
    circ=circularMap(r.name, r.plasmid, [
      feat("att"+cls+"1",0,b,cls), feat(r.insertName||"insert",b,b+m,"O",true),
      feat("att"+cls+"2",b+m,2*b+m,cls),
      feat((v&&v.bact?v.bact.split(" ")[0]:"marker")+"R", 2*b+m+400, 2*b+m+1200,"N",true),
      feat("ori", 2*b+m+1500, 2*b+m+2100,"N",true)]);
  }
  return '<div class="panel"><div class="panel-hd"><h2>'+(isBP?"Entry clone":"Expression clone")+' · '+esc(r.name)+
    '</h2><span class="spacer"></span><span class="pill'+(r.full?" on":"")+'">'+bp(r.plasmid)+
    (r.full?" · complete sequence":" · calculated")+'</span></div><div class="panel-bd stack">'+
    linearMap(feats, 2*armBp+r.mid.length)+
    '<div class="grid2">'+circ+
    '<div class="stack"><dl class="kv">'+
      "<dt>Insert</dt><dd class=\"mono\">"+bp(r.mid.length)+"</dd>"+
      "<dt>Plasmid</dt><dd class=\"mono\">"+bp(r.plasmid)+" ("+kb(r.plasmid)+")</dd>"+
      (r.vec?"<dt>Vector</dt><dd>"+esc(r.vec.id)+"</dd>":"")+
      (r.vec?"<dt>E. coli</dt><dd>"+esc(r.vec.bact)+"</dd>":"")+
      (r.vec&&!isBP?"<dt>Promoter</dt><dd>"+esc(r.vec.prom)+"</dd><dt>Tag</dt><dd>"+esc(r.vec.tag)+
        (r.vec.tagpos!=="none"?" · "+r.vec.tagpos+"-term":"")+"</dd>":"")+
    "</dl>"+
    (r.siteWarn?'<div class="note warn"><b>Vector sequence check:</b> '+esc(r.siteWarn)+"</div>":"")+
    (r.full
      ? '<div class="note ok"><b>Recombined from the real vector sequence.</b> Backbone features were carried across with shifted coordinates; the ccdB cassette was removed.</div>'
      : '<div class="note"><b>Cassette only.</b> The plasmid size is the published vector length minus its att–ccdB cassette. Add this vector’s GenBank file on the Vector library tab to build the complete plasmid.</div>')+
    '<div class="row">'+
      '<button class="btn sm" data-dl="'+r.kind+'" '+(state.caps.dl==="on"?"":"disabled")+'>Download .zip</button>'+
      '<button class="btn ghost sm" data-dlgb="'+r.kind+'" '+(state.caps.dl==="on"?"":"disabled")+'>GenBank only</button>'+
      '<button class="btn ghost sm" data-save="'+r.kind+'" '+(state.caps.db==="on"?"":"disabled")+'>Save to history</button>'+
    "</div>"+
    (state.caps.dl==="off"?'<p class="hint">Downloads are unavailable in this view.</p>':
      '<p class="hint">The zip holds a GenBank file, a FASTA file'+(r.design?", the primer table":"")+' and a README.</p>')+
    "</div></div>"+
    '<div class="legend"><span><i style="background:var(--att'+cls+'-bg);border:1px solid var(--att'+cls+')"></i>att'+cls+' site</span>'+
    '<span><i style="background:var(--orf-bg);border:1px solid var(--orf)"></i>insert</span>'+
    '<span><i style="background:var(--surface2);border:1px solid var(--ink3)"></i>backbone</span></div>'+
    "</div></div>";
}

/* ============================ BP view ============================ */
function viewBP(){
  const c=state.bp, r=reaction("BP");
  const donor=getVec(c.donor);
  let left='<div class="panel"><div class="panel-hd"><h2>1 · Insert</h2></div><div class="panel-bd stack">'+
    '<div class="seg" role="group" aria-label="Insert source">'+
      '<button data-bpmode="design" aria-pressed="'+(c.mode==="design")+'">Design attB primers</button>'+
      '<button data-bpmode="paste" aria-pressed="'+(c.mode==="paste")+'">Paste attB product</button>'+
    "</div>";
  if(c.mode==="design"){
    left+='<label class="fld">Insert name<input type="text" id="bpName" value="'+esc(c.name)+'"></label>'+
      '<label class="fld">ORF / insert sequence (5′→3′)<textarea id="bpOrf" spellcheck="false">'+esc(c.orf)+"</textarea></label>"+
      '<div class="row"><span class="pill">'+bp(clean(c.orf).length)+'</span>'+
        '<span class="pill">'+gcPct(clean(c.orf)).toFixed(1)+'% GC</span>'+
        '<button class="btn ghost sm" id="bpSample">Load EGFP sample</button></div>'+
      '<label class="file"><input type="file" id="bpFile" accept=".fa,.fasta,.fas,.gb,.gbk,.genbank,.txt,.seq">'+
        "<span>⤒ Upload FASTA or GenBank</span></label>";
  } else {
    left+='<label class="fld">attB-flanked PCR product<textarea id="bpProd" spellcheck="false" placeholder="ACAAGTTTGTACAAAAAAGCAGGCT… insert …ACCCAGCTTTCTTGTACAAAGTGGT">'+esc(c.product)+"</textarea></label>"+
      '<label class="file"><input type="file" id="bpFile2" accept=".fa,.fasta,.gb,.gbk,.txt,.seq"><span>⤒ Upload FASTA or GenBank</span></label>'+
      '<p class="hint">Only the sequence between attB1 and attB2 recombines. Everything outside the att sites is discarded.</p>';
  }
  left+="</div></div>";

  if(c.mode==="design"){
    left+='<div class="panel"><div class="panel-hd"><h2>2 · Expression plan</h2></div><div class="panel-bd stack">'+
      '<label class="fld">Protein form'+
        '<select id="bpFormat">'+Object.keys(FORMATS).map(k=>'<option value="'+k+'"'+(c.format===k?" selected":"")+">"+FORMATS[k]+"</option>").join("")+"</select></label>"+
      '<label class="fld">Expression host'+
        '<select id="bpOrg">'+Object.keys(ORGANISMS).map(k=>'<option value="'+k+'"'+(c.organism===k?" selected":"")+">"+ORGANISMS[k]+"</option>").join("")+"</select></label>"+
      '<label class="fld">PCR strategy<div class="seg" role="group">'+
        '<button data-proto="nonstop" aria-pressed="'+(c.protocol==="nonstop")+'">Non-stop PCR</button>'+
        '<button data-proto="adapter" aria-pressed="'+(c.protocol==="adapter")+'">Adapter PCR</button></div></label>'+
      '<p class="hint">Non-stop PCR uses full-length attB primers — simplest for a few reactions. Adapter PCR uses short-tailed gene primers plus universal adapter primers — cheaper in bulk and less prone to primer dimers.</p>'+
      '<label class="check"><input type="checkbox" id="bpTev" '+(c.tev?"checked":"")+'><span>Insert a TEV protease site after the N-terminal tag</span></label>'+
      (MATRIX[c.organism][c.format][0]==="C"
        ? '<label class="check"><input type="checkbox" id="bpCatg" '+(c.cAtg?"checked":"")+'><span>Include my own ATG (optional for form C)</span></label>' : "")+
      "</div></div>";
  }

  left+='<div class="panel"><div class="panel-hd"><h2>'+(c.mode==="design"?"3":"2")+' · Donor vector</h2></div><div class="panel-bd stack">'+
    '<label class="fld">Donor (attP × attP)<select id="bpDonor">'+vecOptions("donor",c.donor)+"</select></label>"+
    (donor?'<dl class="kv"><dt>Size</dt><dd class="mono">'+bp(donor.size)+"</dd>"+
      "<dt>Selection</dt><dd>"+esc(donor.bact)+'</dd><dt>Sites</dt><dd class="mono">'+esc(donor.att)+"</dd>"+
      (donor.seq?'<dt>Sequence</dt><dd class="mono">loaded ✓</dd>':"")+"</dl>"+
      '<p class="hint">'+esc(donor.note||"")+"</p>":"")+
    '<p class="hint">A ✓ in the list means that vector has a full sequence loaded, so the reaction returns a complete plasmid.</p>'+
    "</div></div>";

  let right="";
  if(r.error){
    right='<div class="panel"><div class="panel-bd"><div class="note bad"><b>Cannot run the reaction.</b><br>'+esc(r.error)+"</div></div></div>";
    return '<div class="split"><div>'+left+"</div><div>"+msgBar()+right+"</div></div>";
  }
  right='<div class="panel"><div class="panel-hd"><h2>BP reaction</h2><span class="spacer"></span>'+
    '<span class="pill">BP Clonase II · Int + IHF</span></div><div class="panel-bd stack">'+
    '<div class="eq">'+
      '<div class="box"><div class="t">substrate</div><div class="v"><span class="chip c-B">attB1</span> insert <span class="chip c-B">attB2</span></div></div>'+
      '<span class="op">×</span>'+
      '<div class="box"><div class="t">donor</div><div class="v"><span class="chip c-P">attP1</span> '+esc(c.donor)+' <span class="chip c-P">attP2</span></div></div>'+
      '<span class="op">→</span>'+
      '<div class="box"><div class="t">entry clone</div><div class="v"><span class="chip c-L">attL1</span> insert <span class="chip c-L">attL2</span></div></div>'+
      '<span class="op">+</span>'+
      '<div class="box"><div class="t">by-product</div><div class="v"><span class="chip c-R">attR1</span> ccdB <span class="chip c-R">attR2</span></div></div>'+
    "</div>"+
    '<div class="note ok"><b>attB × attP → attL + attR.</b> Integrase and IHF cut and rejoin inside the 7-bp overlap of the 15-bp core, so the core is preserved and only the flanking arms are exchanged.</div>'+
    "</div></div>";

  if(r.design){
    const d=r.design;
    right+='<div class="panel"><div class="panel-hd"><h2>attB primers</h2><span class="spacer"></span>'+
      '<span class="pill on">form '+d.pair[0]+" + "+d.pair[1]+'</span><span class="pill">'+
      (c.protocol==="adapter"?"adapter PCR":"non-stop PCR")+'</span></div><div class="panel-bd stack">'+
      '<div class="note"><b>'+esc(FORMATS[c.format])+" · "+esc(ORGANISMS[c.organism])+"</b> selects "+
      "attB1 primer <b>"+d.pair[0]+"</b> and attB2 primer <b>"+d.pair[1]+"</b>.</div>"+
      '<div class="grid2">'+
        '<div class="note">'+esc(d.f1.title)+" — "+esc(d.f1.blurb)+"</div>"+
        '<div class="note">'+esc(d.f2.title)+" — "+esc(d.f2.blurb)+"</div>"+
      "</div>"+
      '<div><div class="row" style="justify-content:space-between"><b>Forward — attB1, form '+d.pair[0]+'</b>'+
        '<span class="row"><span class="pill">'+d.fwd.length+' nt</span><span class="pill">gsp Tm '+d.fwdTm.toFixed(1)+'°C</span>'+copyBtn("p-fwd")+"</span></div>"+
        renderSeq(d.fwdParts.map(p=>({seq:p.s,cls:p.c})),70)+
        '<textarea id="p-fwd" hidden>'+d.fwd+"</textarea></div>"+
      '<div><div class="row" style="justify-content:space-between"><b>Reverse — attB2, form '+d.pair[1]+'</b>'+
        '<span class="row"><span class="pill">'+d.rev.length+' nt</span><span class="pill">gsp Tm '+d.revTm.toFixed(1)+'°C</span>'+copyBtn("p-rev")+"</span></div>"+
        renderSeq(d.revParts.map(p=>({seq:p.s,cls:p.c})),70)+
        '<textarea id="p-rev" hidden>'+d.rev+"</textarea></div>"+
      '<div class="legend"><span><i style="background:var(--attB-bg);border:1px solid var(--attB)"></i>attB tail</span>'+
      '<span><i style="background:var(--rrs-bg);border:1px solid var(--rrs)"></i>rrs / start / stop</span>'+
      '<span><i style="background:var(--orf-bg);border:1px solid var(--orf)"></i>gene-specific (18–24 nt)</span></div>'+
      (c.protocol==="adapter"
        ? '<div class="stack" style="gap:7px"><b>Universal adapter primers — second PCR</b>'+
          '<div class="row"><span class="pill">attB1</span><span class="mono" style="font-size:12px">'+ADAPT_B1+"</span>"+copyBtn("p-ad1")+"</div>"+
          '<div class="row"><span class="pill">attB2</span><span class="mono" style="font-size:12px">'+ADAPT_B2+"</span>"+copyBtn("p-ad2")+"</div>"+
          '<textarea id="p-ad1" hidden>'+ADAPT_B1+'</textarea><textarea id="p-ad2" hidden>'+ADAPT_B2+"</textarea>"+
          '<p class="hint">Amplify with the gene primers above (15–20 cycles), then amplify that product with these two adapter primers to build the full attB sites.</p></div>' : "")+
      d.notes.map(n=>'<div class="note '+n.k+'">'+n.t+"</div>").join("")+
      '<ul class="tight"><li>Keep the gene-specific portion 18–24 nt; the attB tail contributes almost nothing to annealing.</li>'+
      "<li>Use a proofreading polymerase and gel- or column-purify the product — a single indel in the core kills recombination.</li>"+
      "<li>TTA, TAG and TGA must not form at the attB1 junction; the tool checks this for you.</li></ul>"+
      "</div></div>";
  }

  right+=productPanel(r);
  right+='<div class="panel"><div class="panel-hd"><h2>'+(r.full?"Entry clone sequence":"Entry cassette sequence")+
    '</h2><span class="spacer"></span><span class="pill">'+bp(r.full?r.full.seq.length:r.cassette.length)+"</span>"+
    copyBtn("seq-entry")+'</div><div class="panel-bd stack">'+
    (r.full
      ? renderSeq([{seq:r.full.seq.slice(0,Math.max(0,r.full.site1.start-1)),cls:"sX"},
                   {seq:r.full.seq.slice(r.full.site1.start-1,r.full.site1.end),cls:"sL"},
                   {seq:r.full.seq.slice(r.full.insert.start-1,r.full.insert.end),cls:"sO"},
                   {seq:r.full.seq.slice(r.full.site2.start-1,r.full.site2.end),cls:"sL"},
                   {seq:r.full.seq.slice(r.full.site2.end),cls:"sX"}],60,4200)
      : renderSeq([{seq:ATT.L1,cls:"sL"},{seq:r.mid,cls:"sO"},{seq:ATT.L2,cls:"sL"}]))+
    '<textarea id="seq-entry" hidden>'+(r.full?r.full.seq:r.cassette)+"</textarea>"+
    '<p class="hint">Teal = attL sites · blue = your insert'+(r.full?" · grey = donor backbone":"")+".</p>"+
    "</div></div>";

  right+='<div class="panel"><div class="panel-hd"><h2>Sequencing primers</h2></div><div class="panel-bd stack">'+
    SEQ_PRIMERS.map((p,i)=>'<div class="row" style="justify-content:space-between"><span><b>'+esc(p.name)+"</b> "+
      '<span class="hint">'+esc(p.note)+'</span></span><span class="row"><span class="mono" style="font-size:12px">'+p.seq+"</span>"+
      copyBtn("sp"+i)+'</span><textarea id="sp'+i+'" hidden>'+p.seq+"</textarea></div>").join("")+
    "</div></div>";

  right+='<div class="panel"><div class="panel-hd"><h2>Bench protocol</h2></div><div class="panel-bd stack">'+
    '<ol class="proto">'+
    "<li>Combine <b>15–150 ng</b> attB PCR product (1–7 µL), <b>150 ng</b> "+esc(c.donor)+" and TE pH 8.0 to <b>8 µL</b>.</li>"+
    "<li>Thaw BP Clonase II on ice, vortex twice for 2 s, add <b>2 µL</b> and mix.</li>"+
    "<li>Incubate <b>1 h at 25 °C</b> — overnight for inserts over 5 kb or low yields.</li>"+
    "<li>Add <b>1 µL proteinase K</b>, 10 min at 37 °C.</li>"+
    "<li>Transform into a <b>ccdB-sensitive</b> strain (DH5α, TOP10, Mach1) and plate on "+esc(donor?donor.bact:"the donor marker")+".</li>"+
    "<li>Sequence across both junctions with SeqL-A and SeqL-B before going on to the LR reaction.</li>"+
    "</ol>"+
    '<div class="note warn"><b>ccdB is the counter-selection.</b> Propagate the donor vector itself only in DB3.1 or ccdB Survival 2 T1<sup>R</sup>. Background colonies almost always mean a ccdB-resistant host.</div>'+
    "</div></div>";
  return '<div class="split"><div>'+left+"</div><div>"+msgBar()+right+"</div></div>";
}

/* ============================ LR view ============================ */
function viewLR(){
  const c=state.lr, r=reaction("LR");
  const dest=getVec(c.dest);
  let left='<div class="panel"><div class="panel-hd"><h2>1 · Entry clone</h2></div><div class="panel-bd stack">'+
    '<div class="seg" role="group" aria-label="Entry clone source">'+
      '<button data-lrmode="carry" aria-pressed="'+(c.mode==="carry")+'">Use BP product</button>'+
      '<button data-lrmode="paste" aria-pressed="'+(c.mode==="paste")+'">Paste or upload entry clone</button>'+
    "</div>";
  if(c.mode==="carry"){
    left+='<p class="hint">Reading the entry clone built on the BP tab. Change the insert there and this follows.</p>'+
      '<label class="fld">Insert name<input type="text" id="lrName" value="'+esc(state.bp.name)+'"></label>';
  } else {
    left+='<label class="fld">Entry clone or attL cassette<textarea id="lrCass" spellcheck="false" placeholder="CAAATAATGATTTTATTTTGACTGATAGTG… attL1 … insert … attL2">'+esc(c.cassette)+"</textarea></label>"+
      '<label class="file"><input type="file" id="lrFile" accept=".fa,.fasta,.gb,.gbk,.txt,.seq"><span>⤒ Upload FASTA or GenBank</span></label>'+
      '<p class="hint">A whole entry plasmid works — the tool finds the att cores itself.</p>'+
      '<label class="fld">Insert name<input type="text" id="lrName" value="'+esc(c.name)+'"></label>';
  }
  left+="</div></div>";
  left+='<div class="panel"><div class="panel-hd"><h2>2 · Destination vector</h2></div><div class="panel-bd stack">'+
    '<label class="fld">Destination (attR × attR)<select id="lrDest">'+vecOptions("dest",c.dest)+"</select></label>"+
    (dest?'<dl class="kv"><dt>Host</dt><dd>'+esc(APPNAMES[dest.app]||dest.app)+"</dd>"+
      "<dt>Promoter</dt><dd>"+esc(dest.prom)+"</dd>"+
      "<dt>Tag</dt><dd>"+esc(dest.tag)+(dest.tagpos!=="none"?" ("+dest.tagpos+"-terminal)":"")+"</dd>"+
      "<dt>E. coli</dt><dd>"+esc(dest.bact)+"</dd><dt>Host sel.</dt><dd>"+esc(dest.sel)+"</dd>"+
      '<dt>Size</dt><dd class="mono">'+bp(dest.size)+"</dd>"+
      (dest.seq?'<dt>Sequence</dt><dd class="mono">loaded ✓</dd>':"")+"</dl>"+
      '<p class="hint">'+esc(dest.note||"")+"</p>":"")+
    "</div></div>";

  let right="";
  if(r.error){
    right='<div class="panel"><div class="panel-bd"><div class="note bad"><b>Cannot run the reaction.</b><br>'+esc(r.error)+"</div></div></div>";
    return '<div class="split"><div>'+left+"</div><div>"+msgBar()+right+"</div></div>";
  }
  right='<div class="panel"><div class="panel-hd"><h2>LR reaction</h2><span class="spacer"></span>'+
    '<span class="pill">LR Clonase II · Int + IHF + Xis</span></div><div class="panel-bd stack">'+
    '<div class="eq">'+
      '<div class="box"><div class="t">entry clone</div><div class="v"><span class="chip c-L">attL1</span> insert <span class="chip c-L">attL2</span></div></div>'+
      '<span class="op">×</span>'+
      '<div class="box"><div class="t">destination</div><div class="v"><span class="chip c-R">attR1</span> '+esc(c.dest)+' <span class="chip c-R">attR2</span></div></div>'+
      '<span class="op">→</span>'+
      '<div class="box"><div class="t">expression clone</div><div class="v"><span class="chip c-B">attB1</span> insert <span class="chip c-B">attB2</span></div></div>'+
      '<span class="op">+</span>'+
      '<div class="box"><div class="t">by-product</div><div class="v"><span class="chip c-P">attP1</span> ccdB <span class="chip c-P">attP2</span></div></div>'+
    "</div>"+
    '<div class="note ok"><b>attL × attR → attB + attP.</b> Excisionase drives the reverse of BP: the insert moves into the destination backbone and the ccdB cassette leaves on the attP by-product.</div>'+
    "</div></div>";
  right+=productPanel(r);

  // reading frame
  const expr=ATT.B1+r.mid+ATT.B2;
  const prot=translate(expr);
  const stopAt=prot.indexOf("*");
  const b1pep=translate(ATT.B1);
  const junctionCodons=Math.floor((ATT.B1.length+r.mid.length)/3);
  const offset=(ATT.B1.length+r.mid.length)%3;
  const junctionCodon=offset===1?(r.mid.slice(-1)+ATT.B2.slice(0,2)):null;
  const inFrameB2=offset===1;
  const notes=[];
  notes.push({k:b1pep.indexOf("TSLYKKAG")===0?"ok":"warn",
    t:"<b>attB1 peptide:</b> <span class=\"mono\">"+esc(b1pep)+"</span> — the scar between an N-terminal tag and your first residue. The attB1 primer sets this frame, so no offset guessing is needed."});
  if(dest&&dest.tagpos==="N"){
    notes.push({k:(stopAt<0||stopAt*3>=ATT.B1.length)?"ok":"bad",
      t:"<b>N-terminal "+esc(dest.tag)+" fusion:</b> the insert must have no ATG of its own and must continue the attB1 frame. "+
        (stopAt<0?"No in-frame stop before attB2.":"First in-frame stop at codon "+(stopAt+1)+".")});
  }
  if(dest&&dest.tagpos==="C"){
    notes.push({k:inFrameB2?"ok":"bad",
      t:inFrameB2
        ? "<b>C-terminal "+esc(dest.tag)+" fusion: frame carries through.</b> The last base of the insert plus the first two of attB2 make the codon <span class=\"mono\">"+esc(junctionCodon)+"</span> ("+translate(junctionCodon)+"), exactly the geometry reverse primer form L is built for."
        : "<b>C-terminal "+esc(dest.tag)+" fusion: frame is wrong.</b> The insert ends "+(offset===0?"flush with a codon boundary":"two bases into a codon")+", so the tag would be translated "+(offset===0?"+1":"−1")+" base out of frame. Rebuild the insert with reverse primer form L, which adds the one base that completes the codon running into attB2."});
    notes.push({k:(stopAt<0||stopAt>=junctionCodons)?"ok":"bad",
      t:(stopAt<0||stopAt>=junctionCodons)
        ? "<b>No premature stop</b> — translation reaches attB2 and carries on into the tag."
        : "<b>Stop codon at codon "+(stopAt+1)+"</b> ends translation before attB2, so the C-terminal tag will not be made. Rebuild the insert with reverse primer form L (C-terminal fusion)."});
  }
  if(!dest||dest.tagpos!=="C"){
    notes.push({k:"ok",t:"<b>Translation stops inside the insert</b>"+(stopAt<0?" — no in-frame stop was found, so the ribosome will read into attB2 and beyond.":" at codon "+(stopAt+1)+", "+(junctionCodons-stopAt)+" codons before attB2. Nothing downstream of the insert is translated.")});
  }
  if(dest&&dest.tagpos==="none"){
    notes.push({k:r.mid.indexOf("ATG")>=0?"ok":"warn",
      t:"<b>Untagged vector</b> — the insert must supply its own start codon"+(dest.app==="bacterial"?" and a ribosome-binding site (primer form A)":"")+"."});
  }
  const shown=prot.slice(0,stopAt<0?prot.length:stopAt);
  right+='<div class="panel"><div class="panel-hd"><h2>Reading-frame check</h2><span class="spacer"></span>'+
    '<span class="pill'+(offset===1?" on":"")+'">'+(offset===1?"read-through frame ✓":offset===0?"terminates at attB2":"frame shifted")+'</span></div><div class="panel-bd stack">'+
    notes.map(n=>'<div class="note '+n.k+'">'+n.t+"</div>").join("")+
    '<div><div class="row" style="justify-content:space-between"><b>Predicted translation from attB1</b>'+
    '<span class="row"><span class="pill">'+shown.length+' aa</span>'+copyBtn("aa-prot")+"</span></div>"+
    '<div class="seqwrap"><div class="aa">'+esc(shown.replace(/(.{10})/g,"$1 "))+"</div></div>"+
    '<textarea id="aa-prot" hidden>'+shown+"</textarea></div>"+
    "</div></div>";

  right+='<div class="panel"><div class="panel-hd"><h2>'+(r.full?"Expression clone sequence":"Expression cassette sequence")+
    '</h2><span class="spacer"></span><span class="pill">'+bp(r.full?r.full.seq.length:expr.length)+"</span>"+
    copyBtn("seq-expr")+'</div><div class="panel-bd stack">'+
    (r.full
      ? renderSeq([{seq:r.full.seq.slice(0,Math.max(0,r.full.site1.start-1)),cls:"sX"},
                   {seq:r.full.seq.slice(r.full.site1.start-1,r.full.site1.end),cls:"sB"},
                   {seq:r.full.seq.slice(r.full.insert.start-1,r.full.insert.end),cls:"sO"},
                   {seq:r.full.seq.slice(r.full.site2.start-1,r.full.site2.end),cls:"sB"},
                   {seq:r.full.seq.slice(r.full.site2.end),cls:"sX"}],60,4200)
      : renderSeq([{seq:ATT.B1,cls:"sB"},{seq:r.mid,cls:"sO"},{seq:ATT.B2,cls:"sB"}]))+
    '<textarea id="seq-expr" hidden>'+(r.full?r.full.seq:expr)+"</textarea>"+
    '<p class="hint">Amber = attB sites (25 bp each) · blue = your insert'+(r.full?" · grey = destination backbone":"")+".</p>"+
    "</div></div>";

  right+='<div class="panel"><div class="panel-hd"><h2>Bench protocol</h2></div><div class="panel-bd stack">'+
    '<ol class="proto">'+
    "<li>Combine <b>50–150 ng</b> entry clone, <b>150 ng</b> "+esc(c.dest)+" and TE pH 8.0 to <b>8 µL</b>.</li>"+
    "<li>Add <b>2 µL LR Clonase II</b> (vortexed twice, kept on ice) and mix.</li>"+
    "<li>Incubate <b>1 h at 25 °C</b>"+(dest&&dest.size>10000?" — use <b>overnight at 25 °C</b> for this vector; large destinations recombine slowly.":" — overnight for larger constructs.")+"</li>"+
    "<li>Add <b>1 µL proteinase K</b>, 10 min at 37 °C.</li>"+
    "<li>Transform into a ccdB-sensitive strain"+(dest&&(dest.app==="lentiviral"||dest.app==="retroviral")?" (use a <b>recombination-deficient</b> host such as Stbl3 — LTRs are unstable in DH5α)":"")+
      " and select on "+esc(dest?dest.bact:"the vector marker")+".</li>"+
    "<li>Verify both attB junctions by sequencing before moving into cells.</li></ol>"+
    (dest&&(dest.app==="lentiviral"||dest.app==="retroviral"||dest.app==="adenoviral")
      ? '<div class="note warn"><b>Viral backbone.</b> Grow at 30 °C in a recombination-deficient strain and check the LTRs or ITRs by restriction digest.</div>':"")+
    (dest?'<div class="note"><b>Entry and destination must carry different markers.</b> '+esc(dest.id)+" is "+esc(dest.bact)+".</div>":"")+
    "</div></div>";
  return '<div class="split"><div>'+left+"</div><div>"+msgBar()+right+"</div></div>";
}

/* ============================ library view ============================ */
function viewLib(){
  const f=state.lib, d=f.draft;
  const rows=allVectors().filter(v=>
    (f.kind==="all"||v.kind===f.kind) && (f.app==="all"||v.app===f.app) &&
    (!f.q || (v.id+" "+v.tag+" "+v.prom+" "+v.src+" "+(v.note||"")+" "+v.bact).toLowerCase().includes(f.q.toLowerCase())));
  const sel=allVectors().find(v=>v.id===f.sel)||rows[0];
  const apps=[...new Set(allVectors().map(v=>v.app))];
  let s=msgBar();
  s+='<div class="panel"><div class="panel-hd"><h2>Vector library</h2><span class="spacer"></span>'+
    '<span class="pill">'+rows.length+" of "+allVectors().length+"</span>"+
    '<button class="btn sm" id="addToggle">'+(f.adding?"Close":"+ Add vector")+"</button></div>"+
    '<div class="panel-bd stack">';
  if(f.adding){
    s+='<div class="stack" style="border:1px dashed var(--line);border-radius:10px;padding:13px">'+
      "<b>Add a plasmid from Addgene, or any Gateway vector of your own</b>"+
      '<p class="hint">On the plasmid’s Addgene page open <b>Sequences</b> → <b>Download full sequence</b> (GenBank), then drop the file here. The att sites, ccdB cassette, resistance marker and promoter are read straight out of the file.</p>'+
      '<label class="file"><input type="file" id="vecFile" accept=".fa,.fasta,.gb,.gbk,.genbank,.txt,.seq"><span>⤒ Upload the plasmid’s GenBank or FASTA file</span></label>'+
      (d.rec?'<div class="note ok"><b>'+esc(d.rec.name)+"</b> — "+bp(d.rec.seq.length)+" · "+(d.rec.features||[]).length+
        " features"+(d.att&&d.att!=="—"?' · <span class="mono">'+esc(d.att)+"</span>":"")+"</div>":"")+
      '<div class="grid2">'+
        '<label class="fld">Name<input type="text" id="dName" value="'+esc(d.name)+'" placeholder="pLX304"></label>'+
        '<label class="fld">Addgene number<input type="text" id="dAdd" value="'+esc(d.addgene)+'" placeholder="25890"></label>'+
        '<label class="fld">Role<select id="dKind">'+
          ['donor','entry','dest'].map(k=>'<option value="'+k+'"'+(d.kind===k?" selected":"")+">"+
            {donor:"Donor (attP)",entry:"Entry (attL)",dest:"Destination (attR)"}[k]+"</option>").join("")+"</select></label>"+
        '<label class="fld">Host<select id="dApp">'+Object.keys(APPNAMES).map(k=>'<option value="'+k+'"'+(d.app===k?" selected":"")+">"+APPNAMES[k]+"</option>").join("")+"</select></label>"+
        '<label class="fld">Promoter<input type="text" id="dProm" value="'+esc(d.prom)+'" placeholder="CMV"></label>'+
        '<label class="fld">Tag<input type="text" id="dTag" value="'+esc(d.tag)+'" placeholder="V5"></label>'+
        '<label class="fld">Tag position<select id="dTagpos">'+["none","N","C"].map(k=>'<option value="'+k+'"'+(d.tagpos===k?" selected":"")+">"+
          {none:"none",N:"N-terminal",C:"C-terminal"}[k]+"</option>").join("")+"</select></label>"+
        '<label class="fld">E. coli selection<input type="text" id="dBact" value="'+esc(d.bact)+'"></label>'+
        '<label class="fld">Host selection<input type="text" id="dSel" value="'+esc(d.sel)+'"></label>'+
      "</div>"+
      '<label class="fld">Notes<input type="text" id="dNote" value="'+esc(d.note)+'" placeholder="What this vector is for"></label>'+
      '<div class="row"><button class="btn sm" id="dSave">Add to library</button>'+
      '<button class="btn ghost sm" id="dCancel">Cancel</button>'+
      (state.caps.db==="on"?'<span class="hint">Saved with this page, so it is here next time.</span>'
        :'<span class="hint">The data store is unavailable — this vector will last only for this session.</span>')+"</div></div>";
  }
  s+='<div class="filters">'+
    '<select id="fKind"><option value="all">All roles</option>'+
      ['donor','entry','dest'].map(k=>'<option value="'+k+'"'+(f.kind===k?" selected":"")+">"+
        {donor:"Donor (attP)",entry:"Entry (attL)",dest:"Destination (attR)"}[k]+"</option>").join("")+"</select>"+
    '<select id="fApp"><option value="all">All hosts</option>'+
      apps.map(a=>'<option value="'+a+'"'+(f.app===a?" selected":"")+">"+esc(APPNAMES[a]||a)+"</option>").join("")+"</select>"+
    '<input type="text" id="fQ" placeholder="Search name, tag, promoter…" value="'+esc(f.q)+'">'+
  "</div>"+
  '<div class="tablewrap"><table><thead><tr><th>Vector</th><th>Role</th><th>Host</th><th>Promoter</th><th>Tag</th>'+
  "<th>E. coli</th><th>Host sel.</th><th>Size</th></tr></thead><tbody>"+
  (rows.length? rows.map(v=>'<tr data-vec="'+esc(v.id)+'" aria-selected="'+(sel&&v.id===sel.id)+'">'+
    '<td class="vname">'+esc(v.id)+(v.seq?' <span class="pill on">seq</span>':"")+(v.custom?' <span class="pill">added</span>':"")+"</td>"+
    "<td>"+{donor:'<span class="chip c-P">attP</span>',entry:'<span class="chip c-L">attL</span>',dest:'<span class="chip c-R">attR</span>'}[v.kind]+"</td>"+
    "<td>"+esc(APPNAMES[v.app]||v.app)+"</td><td>"+esc(v.prom)+"</td>"+
    "<td>"+esc(v.tag)+(v.tagpos!=="none"?' <span class="pill">'+v.tagpos+"</span>":"")+"</td>"+
    "<td>"+esc(v.bact)+"</td><td>"+esc(v.sel)+'</td><td class="num">'+bp(v.size)+"</td></tr>").join("")
    : '<tr><td colspan="8" class="empty">No vectors match those filters.</td></tr>')+
  "</tbody></table></div></div></div>";

  if(sel){
    s+='<div class="panel"><div class="panel-hd"><h2>'+esc(sel.id)+'</h2><span class="spacer"></span>'+
      '<span class="chip c-'+{donor:"P",entry:"L",dest:"R"}[sel.kind]+'">'+esc(sel.att)+"</span></div>"+
      '<div class="panel-bd stack"><div class="grid2">'+
      '<dl class="kv"><dt>Role</dt><dd>'+{donor:"Donor — accepts attB inserts in a BP reaction",
        entry:"Entry — donates the insert in an LR reaction",dest:"Destination — receives the insert in an LR reaction"}[sel.kind]+"</dd>"+
      "<dt>Host</dt><dd>"+esc(APPNAMES[sel.app]||sel.app)+'</dd><dt>Size</dt><dd class="mono">'+bp(sel.size)+"</dd>"+
      "<dt>Promoter</dt><dd>"+esc(sel.prom)+"</dd><dt>Tag</dt><dd>"+esc(sel.tag)+(sel.tagpos!=="none"?" ("+sel.tagpos+"-terminal)":"")+"</dd>"+
      "<dt>E. coli</dt><dd>"+esc(sel.bact)+"</dd><dt>Host sel.</dt><dd>"+esc(sel.sel)+"</dd>"+
      "<dt>Source</dt><dd>"+esc(sel.src)+(sel.addgene?' · <a href="https://www.addgene.org/'+encodeURIComponent(sel.addgene)+'/" target="_blank" rel="noopener">Addgene #'+esc(sel.addgene)+"</a>":"")+"</dd>"+
      (sel.seq?'<dt>Sequence</dt><dd class="mono">'+bp(sel.seq.length)+" · "+(sel.features||[]).length+" features</dd>":"")+
      "</dl>"+
      '<div class="stack"><p style="margin:0">'+esc(sel.note||"")+"</p>"+
      (sel.kind==="donor"?'<button class="btn sm" data-use="donor" data-id="'+esc(sel.id)+'">Use in BP reaction</button>':"")+
      (sel.kind==="dest"?'<button class="btn sm" data-use="dest" data-id="'+esc(sel.id)+'">Use in LR reaction</button>':"")+
      (sel.custom?'<button class="btn danger sm" data-delvec="'+esc(sel.id)+'">Remove from library</button>':"")+
      (!sel.seq?'<div class="note"><b>No sequence loaded.</b> Add this vector’s GenBank file to get complete recombined plasmids and downloads.</div>'
        :'<div class="note ok"><b>Sequence loaded.</b> Reactions using this vector return the complete plasmid.</div>')+
      ((sel.kind==="donor"||sel.kind==="dest")?'<div class="note warn"><b>ccdB vector.</b> Maintain in DB3.1 or ccdB Survival 2 T1<sup>R</sup> only.</div>':"")+
      "</div></div></div></div>";
  }
  return s;
}

/* ============================ history view ============================ */
function viewHist(){
  let s=msgBar();
  s+='<div class="panel"><div class="panel-hd"><h2>Saved runs</h2><span class="spacer"></span>'+
    '<span class="pill">'+state.history.length+'</span></div><div class="panel-bd">';
  if(state.caps.db==="pending") s+='<div class="empty">Connecting to this page’s data store…</div>';
  else if(state.caps.db==="off") s+='<div class="note warn"><b>History is unavailable in this view.</b> Saved runs are kept with the artifact itself, which needs the page to be opened from its own link while signed in.</div>';
  else if(!state.history.length) s+='<div class="empty">Nothing saved yet. Run a BP or LR reaction and choose <b>Save to history</b>.</div>';
  else {
    s+='<div class="hist">'+state.history.map(h=>{
      const when=new Date(h.createdAt);
      return '<div class="hrow"><div class="meta"><div><b>'+esc(h.name)+"</b> "+
        '<span class="chip c-'+(h.kind==="BP"?"L":"B")+'">'+h.kind+"</span></div>"+
        '<div class="when">'+when.toLocaleString()+" · "+esc(h.vector)+" · insert "+bp(h.insertBp)+
        " · "+bp(h.plasmidBp)+(h.complete?" complete":" calculated")+"</div></div>"+
        '<div class="row"><button class="btn ghost sm" data-restore="'+esc(h.docId)+'">Restore</button>'+
        (h.fwd?'<button class="btn ghost sm" data-hprimers="'+esc(h.docId)+'">Primers</button>':"")+
        '<button class="btn danger sm" data-delrun="'+esc(h.docId)+'">Delete</button></div></div>';
    }).join("")+"</div>";
  }
  s+="</div></div>";
  s+='<div class="panel"><div class="panel-bd stack">'+
    '<div class="note"><b>What a saved run keeps:</b> the reaction type, the insert sequence between the att cores, the vector, both calculated sizes and the primer pair. Restoring loads that insert back into the matching tab, so a construct from last month picks up where it left off.</div>'+
    (state.caps.db==="on"?'<p class="hint">History is stored with this artifact and is visible to anyone it is shared with.</p>':"")+
    "</div></div>";
  return s;
}

/* ============================ att view ============================ */
function attRow(label,cls,seq,core,desc){
  const i=seq.indexOf(core);
  const segs=i<0?[{seq,cls}]:[{seq:seq.slice(0,i),cls},{seq:core,cls:cls+" sC"},{seq:seq.slice(i+core.length),cls}];
  const id="att-"+label.replace(/\W/g,"");
  return '<div class="stack" style="gap:7px"><div class="row" style="justify-content:space-between">'+
    "<b>"+label+'</b><span class="row"><span class="pill">'+seq.length+" bp</span>"+copyBtn(id)+"</span></div>"+
    '<p class="hint">'+desc+"</p>"+renderSeq(segs,60)+
    '<textarea id="'+id+'" hidden>'+seq+"</textarea></div>";
}
function viewAtt(){
  let s='<div class="panel"><div class="panel-hd"><h2>How the chemistry works</h2></div><div class="panel-bd stack">'+
    '<div class="grid2">'+
    '<div class="stack"><b>BP reaction</b>'+
      '<div class="eq"><span class="chip c-B">attB</span><span class="op">×</span><span class="chip c-P">attP</span>'+
      '<span class="op">→</span><span class="chip c-L">attL</span><span class="op">+</span><span class="chip c-R">attR</span></div>'+
      '<p class="hint">BP Clonase II = λ <b>Integrase</b> + <i>E. coli</i> <b>IHF</b>. Moves an attB-flanked PCR product into a donor vector to make an entry clone.</p></div>'+
    '<div class="stack"><b>LR reaction</b>'+
      '<div class="eq"><span class="chip c-L">attL</span><span class="op">×</span><span class="chip c-R">attR</span>'+
      '<span class="op">→</span><span class="chip c-B">attB</span><span class="op">+</span><span class="chip c-P">attP</span></div>'+
      '<p class="hint">LR Clonase II adds λ <b>Excisionase</b>. Moves the insert from an entry clone into any destination vector.</p></div>'+
    "</div>"+
    '<div class="note"><b>Specificity lives in the core.</b> Strand exchange happens inside a 7-bp overlap at the centre of a 15-bp core. The att1 and att2 sites differ inside that overlap and in the bases flanking it, so attB1 only ever reacts with attP1 — that is what keeps the insert from going in backwards.</div>'+
    (function(){
      const a=ATT.B1, b=rc(ATT.B2);
      let m=""; for(let i=0;i<a.length;i++) m+=(a[i]===b[i]?" ":"^");
      const i1=a.indexOf(CORE[1]);
      let core=""; for(let i=0;i<a.length;i++) core+=(i>=i1&&i<i1+15)?"─":" ";
      return '<div class="seqwrap"><div class="seq">'+
        'attB1        <span class="sB">'+a+"</span>\n"+
        'attB2 (rev)  <span class="sB">'+b+"</span>\n"+
        '             <span class="n">'+m+"</span>  mismatches\n"+
        '             <span class="n">'+core+"</span>  15-bp core\n"+
      "</div></div>";
    })()+"</div></div>";

  s+='<div class="panel"><div class="panel-hd"><h2>Site sequences</h2><span class="spacer"></span>'+
    '<span class="hint">dotted underline = 15-bp core</span></div><div class="panel-bd stack" style="gap:20px">'+
    attRow("attB1","sB",ATT.B1,CORE[1],"On the forward PCR primer, after four 5′ G residues. Translated in frame it reads TSLYKKAGS.")+
    attRow("attB2","sB",ATT.B2,CORE[2],"Top strand as it appears in the expression clone; the reverse primer carries its reverse complement.")+
    attRow("attL1","sL",ATT.L1,CORE[1],"Formed in the entry clone. attB1 is retained at its 3′ end and the λ arm is added upstream.")+
    attRow("attL2","sL",ATT.L2,CORE[2],"The attL2 end of the entry cassette — the mirror image of attL1.")+
    attRow("attR1","sR",ATT.R1,CORE[1],"Flanks the ccdB–Cmᴿ cassette in every destination vector.")+
    attRow("attR2","sR",ATT.R2,CORE[2],"The downstream side of the destination cassette.")+
    '<div class="note"><b>attP1 / attP2</b> are the ~200 bp λ sites in donor vectors. They share the same cores but carry both λ arms, so they are not reproduced here in full — upload your donor’s GenBank file and the tool reads them directly.</div>'+
    "</div></div>";

  s+='<div class="panel"><div class="panel-hd"><h2>Primer forms</h2><span class="spacer"></span>'+
    '<span class="hint">after Won Do Heo’s attB primer-design protocol</span></div><div class="panel-bd stack">'+
    '<div class="tablewrap"><table style="min-width:560px"><thead><tr><th>Form</th><th>Side</th><th>Adds</th><th>Use</th></tr></thead><tbody>'+
    Object.keys(FORMS).map(k=>{const f=FORMS[k];
      return "<tr><td class=\"vname\">"+k+'</td><td>'+(f.side==="B1"?"attB1 (forward)":"attB2 (reverse)")+
      '</td><td class="mono">'+(f.add||f.extraBase||"—")+"</td><td>"+esc(f.title.replace(/^.*— /,""))+"</td></tr>";}).join("")+
    "</tbody></table></div>"+
    '<div class="tablewrap"><table style="min-width:620px"><thead><tr><th>Host</th>'+
      Object.keys(FORMATS).map(k=>"<th>"+FORMATS[k]+"</th>").join("")+"</tr></thead><tbody>"+
      Object.keys(ORGANISMS).map(o=>"<tr><td>"+ORGANISMS[o]+"</td>"+
        Object.keys(FORMATS).map(fm=>'<td class="mono">'+MATRIX[o][fm].join(" + ")+"</td>").join("")+"</tr>").join("")+
    "</tbody></table></div>"+
    '<ul class="tight">'+
    "<li>The start codon and ribosome-recognition sequence come from the destination vector for an N-terminal fusion, and from the PCR product otherwise.</li>"+
    "<li>The stop codon comes from the vector for a C-terminal fusion, and from the reverse primer otherwise.</li>"+
    "<li>TTA, TAG and TGA must never form at the attB1 junction.</li>"+
    "<li>Form L adds one base so the codon running into attB2 is complete; the gene-specific sequence must end on a codon boundary.</li>"+
    "</ul></div></div>";

  s+='<div class="panel"><div class="panel-hd"><h2>What ends up in your protein</h2></div><div class="panel-bd stack">'+
    '<div class="tablewrap"><table style="min-width:520px"><thead><tr><th>Site</th><th>Nucleotide</th><th>In-frame peptide</th><th>Consequence</th></tr></thead><tbody>'+
    '<tr><td class="vname">attB1</td><td class="mono">'+ATT.B1+'</td><td class="mono">'+translate(ATT.B1)+"</td><td>Linker between an N-terminal tag and your first residue.</td></tr>"+
    '<tr><td class="vname">attB2</td><td class="mono">'+ATT.B2+'</td><td class="mono">'+translate(ATT.B2)+"</td><td>Linker to a C-terminal tag — translated only when the insert has no stop codon.</td></tr>"+
    "</tbody></table></div>"+
    '<ul class="tight">'+
    "<li>The attB scar is unavoidable in Gateway. If an untouched terminus matters, use a scarless cloning method.</li>"+
    "<li>attB1 is 25 bp — not a multiple of three — so the forward primer adds TC to bring the ORF into frame.</li>"+
    "<li>Both att sites are AT-rich, so the gene-specific portion carries the annealing.</li>"+
    "</ul></div></div>";
  return s;
}

/* ============================ wiring ============================ */
const VIEWS={bp:viewBP, lr:viewLR, lib:viewLib, hist:viewHist, att:viewAtt};
const TABS=["bp","lr","lib","hist","att"];
function render(tab){ const el=document.getElementById("view-"+tab); if(el) el.innerHTML=VIEWS[tab](); }
function rerenderAll(){ TABS.forEach(render); }
function switchTab(t){
  state.tab=t;
  document.querySelectorAll("#tabs .tab").forEach(b=>b.setAttribute("aria-selected",String(b.dataset.tab===t)));
  TABS.forEach(k=>{const el=document.getElementById("view-"+k); if(el) el.hidden=(k!==t);});
  render(t);
}
document.getElementById("tabs").addEventListener("click",e=>{
  const b=e.target.closest(".tab"); if(b) switchTab(b.dataset.tab);
});
document.getElementById("themeBtn").addEventListener("click",()=>{
  const cur=document.documentElement.getAttribute("data-theme");
  const dark=cur?cur==="dark":matchMedia("(prefers-color-scheme: dark)").matches;
  document.documentElement.setAttribute("data-theme",dark?"light":"dark");
});
document.addEventListener("click",e=>{
  const cp=e.target.closest("[data-copy]");
  if(cp){
    const el=document.getElementById(cp.dataset.copy); if(!el) return;
    const txt=el.value!==undefined?el.value:el.textContent;
    const done=()=>{const o=cp.textContent;cp.textContent="Copied";cp.classList.add("copied");
      setTimeout(()=>{cp.textContent=o;cp.classList.remove("copied");},1400);};
    try{ navigator.clipboard.writeText(txt).then(done,()=>{el.hidden=false;el.select();}); }
    catch(_){ el.hidden=false; el.select(); }
    return;
  }
  let t;
  if((t=e.target.closest("[data-bpmode]"))){ state.bp.mode=t.dataset.bpmode; rerenderAll(); return; }
  if((t=e.target.closest("[data-lrmode]"))){ state.lr.mode=t.dataset.lrmode; render("lr"); return; }
  if((t=e.target.closest("[data-proto]"))){ state.bp.protocol=t.dataset.proto; rerenderAll(); return; }
  if((t=e.target.closest("tr[data-vec]"))){ state.lib.sel=t.dataset.vec; render("lib"); return; }
  if((t=e.target.closest("[data-use]"))){
    if(t.dataset.use==="donor"){ state.bp.donor=t.dataset.id; switchTab("bp"); }
    else { state.lr.dest=t.dataset.id; switchTab("lr"); } return;
  }
  if((t=e.target.closest("[data-dl]"))){ downloadBundle(reaction(t.dataset.dl)); return; }
  if((t=e.target.closest("[data-dlgb]"))){ downloadGb(reaction(t.dataset.dlgb)); return; }
  if((t=e.target.closest("[data-save]"))){ saveRun(reaction(t.dataset.save)); render("hist"); return; }
  if((t=e.target.closest("[data-delvec]"))){ deleteCustomVector(t.dataset.delvec); return; }
  if((t=e.target.closest("[data-delrun]"))){ deleteRun(t.dataset.delrun); return; }
  if((t=e.target.closest("[data-restore]"))){
    const h=state.history.find(x=>x.docId===t.dataset.restore); if(h) restoreRun(h); return;
  }
  if((t=e.target.closest("[data-hprimers]"))){
    const h=state.history.find(x=>x.docId===t.dataset.hprimers);
    if(h) flash("<b>"+esc(h.name)+"</b><br><span class='mono' style='font-size:11.5px'>F "+esc(h.fwd)+"<br>R "+esc(h.rev)+"</span>","ok");
    return;
  }
  if(e.target.id==="bpSample"){ state.bp.orf=SAMPLE_ORF; state.bp.name="EGFP"; rerenderAll(); return; }
  if(e.target.id==="addToggle"){ state.lib.adding=!state.lib.adding; render("lib"); return; }
  if(e.target.id==="dCancel"){ state.lib.adding=false; state.lib.draft.rec=null; render("lib"); return; }
  if(e.target.id==="dSave"){ saveCustomVector(); return; }
});
document.addEventListener("change",e=>{
  const id=e.target.id;
  if(e.target.type==="file"){
    const f=e.target.files&&e.target.files[0]; if(!f) return;
    handleSeqFile(f, id==="vecFile"?"vector":(id==="lrFile"?"entry":"orf"));
    e.target.value=""; return;
  }
  const d=state.lib.draft;
  const M={bpDonor:()=>state.bp.donor=e.target.value, bpFormat:()=>state.bp.format=e.target.value,
    bpOrg:()=>state.bp.organism=e.target.value, bpTev:()=>state.bp.tev=e.target.checked,
    bpCatg:()=>state.bp.cAtg=e.target.checked, lrDest:()=>state.lr.dest=e.target.value,
    fKind:()=>state.lib.kind=e.target.value, fApp:()=>state.lib.app=e.target.value,
    dKind:()=>d.kind=e.target.value, dApp:()=>d.app=e.target.value, dTagpos:()=>d.tagpos=e.target.value};
  if(M[id]){ M[id]();
    if(id==="lrDest") render("lr");
    else if(id==="fKind"||id==="fApp"||id==="dKind"||id==="dApp"||id==="dTagpos") render("lib");
    else rerenderAll();
  }
});
let reTimer=null, pending=null;
function scheduleRe(tab){
  pending=tab; clearTimeout(reTimer);
  reTimer=setTimeout(()=>{
    const a=document.activeElement, id=a&&a.id, ss=a&&a.selectionStart, se=a&&a.selectionEnd;
    if(pending==="all") rerenderAll(); else { render(pending); if(pending==="bp") render("lr"); }
    if(id){ const n=document.getElementById(id); if(n){ n.focus(); try{ n.setSelectionRange(ss,se); }catch(_){}} }
  },320);
}
document.addEventListener("input",e=>{
  const id=e.target.id, v=e.target.value, d=state.lib.draft;
  const B={bpOrf:()=>state.bp.orf=v, bpProd:()=>state.bp.product=v, bpName:()=>state.bp.name=v};
  const L={lrCass:()=>state.lr.cassette=v,
    lrName:()=>{ if(state.lr.mode==="carry") state.bp.name=v; else state.lr.name=v; }};
  const V={fQ:()=>state.lib.q=v, dName:()=>d.name=v, dAdd:()=>d.addgene=v, dProm:()=>d.prom=v,
    dTag:()=>d.tag=v, dBact:()=>d.bact=v, dSel:()=>d.sel=v, dNote:()=>d.note=v};
  if(B[id]){ B[id](); scheduleRe("bp"); }
  else if(L[id]){ L[id](); scheduleRe("lr"); }
  else if(V[id]){ V[id](); scheduleRe("lib"); }
});

switchTab("bp");
TABS.forEach(render);
boot();
