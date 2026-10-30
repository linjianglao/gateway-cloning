/* ============================ att site sequences ============================ */
const ATT = {
  B1:"ACAAGTTTGTACAAAAAAGCAGGCT",
  B2:"ACCCAGCTTTCTTGTACAAAGTGGT",
  L1:"CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCT",
  L2:"ACCCAGCTTTCTTGTACAAAGTTGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTG",
  R1:"ACAAGTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATG",
  R2:"CATAGTGACTGGATATGTTGTGTTTTACAGTATTATGTAGTCTGTTTTTTATGCAAAATCTAATTTAATATATTGATATTTATATCATTTTACGTTTCTCGTTCAGCTTTCTTGTACAAAGTGGT"
};
const CORE = {1:"TTTGTACAAAAAAGC", 2:"TTTGTACAAGAAAGC"};
const SEQ_PRIMERS = [
  {name:"SeqL-A", seq:"TCGCGTTAACGCTAGCATGGATCTC", note:"proximal to attL1 · 106 nt from the cloned ORF"},
  {name:"SeqL-B", seq:"GTAACATCAGAGATTTTGAGACAC", note:"proximal to attL2 · 123 nt from the cloned ORF"},
  {name:"M13 forward", seq:"GTAAAACGACGGCCAG", note:"pDONR221 and most pENTR backbones"},
  {name:"M13 reverse", seq:"CAGGAAACAGCTATGAC", note:"pDONR221 and most pENTR backbones"}
];

/* ============================ sample ORF ============================ */
const SAMPLE_ORF =
"ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGC"+
"GAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGC"+
"GTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAG"+
"GACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAAC"+
"ATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCAC"+
"AACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACC"+
"CAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTAC"+
"AAGTAA";

/* ============================ sequence utils ============================ */
const CODON_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const B4 = "TCAG";
function clean(s){return (s||"").toUpperCase().replace(/[^ACGTUN]/g,"").replace(/U/g,"T");}
function rc(s){const m={A:"T",T:"A",G:"C",C:"G",N:"N"};return s.split("").reverse().map(c=>m[c]||"N").join("");}
function translate(s){let out="";for(let i=0;i+2<s.length;i+=3){const a=B4.indexOf(s[i]),b=B4.indexOf(s[i+1]),c=B4.indexOf(s[i+2]);
  out += (a<0||b<0||c<0) ? "X" : CODON_AA[16*a+4*b+c];}return out;}
function gcPct(s){if(!s.length)return 0;return 100*(s.match(/[GC]/g)||[]).length/s.length;}
function tm(s){const n=s.length;if(!n)return 0;const gc=(s.match(/[GC]/g)||[]).length;
  return n<14 ? 2*(n-gc)+4*gc : 64.9+41*(gc-16.4)/n;}
function bp(n){return n.toLocaleString("en-US")+" bp";}
function kb(n){return (n/1000).toFixed(n<10000?2:1)+" kb";}
function esc(s){return String(s).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));}
function slug(s){return String(s).replace(/[^A-Za-z0-9_.+~:@-]+/g,"-").replace(/^-+|-+$/g,"").slice(0,60)||"item";}

function renderSeq(segs, perLine, cap){
  perLine=perLine||60; cap=cap||6000;
  const cls=[]; let all="";
  for(const s of segs){ for(let i=0;i<s.seq.length;i++) cls.push(s.cls||"sX"); all+=s.seq; }
  let trimmed=false;
  if(all.length>cap){ all=all.slice(0,cap); cls.length=cap; trimmed=true; }
  const width=String(all.length).length;
  let html="";
  for(let i=0;i<all.length;i+=perLine){
    const line=all.slice(i,i+perLine), lc=cls.slice(i,i+perLine);
    let body="",cur=null,buf="";
    const flush=()=>{ if(buf){ body+='<span class="'+cur+'">'+buf+"</span>"; buf=""; } };
    for(let j=0;j<line.length;j++){
      if(j&&j%10===0){ flush(); body+=" "; }
      if(lc[j]!==cur){ flush(); cur=lc[j]; }
      buf+=line[j];
    }
    flush();
    html+='<span class="n">'+String(i+1).padStart(width," ")+"  </span>"+body+"\n";
  }
  if(trimmed) html+='<span class="n">… truncated for display — download the file for the complete sequence</span>\n';
  return '<div class="seqwrap"><div class="seq">'+html+"</div></div>";
}

/* ============================ maps ============================ */
function linearMap(feats, total){
  const W=920, H=70, pad=14, y=30, h=16;
  let s='<svg class="svgmap" viewBox="0 0 '+W+' '+H+'" role="img" aria-label="Linear construct map">';
  s+='<line x1="'+pad+'" y1="'+(y+h/2)+'" x2="'+(W-pad)+'" y2="'+(y+h/2)+'" stroke="var(--line)" stroke-width="2"></line>';
  for(const f of feats){
    const x1=pad+(W-2*pad)*(f.start/total), x2=pad+(W-2*pad)*(f.end/total);
    const w=Math.max(x2-x1,2.5);
    if(f.arrow&&w>16){
      const t=Math.min(11,w*0.28);
      s+='<path d="M'+x1+' '+y+' H'+(x1+w-t)+' L'+(x1+w)+' '+(y+h/2)+' L'+(x1+w-t)+' '+(y+h)+' H'+x1+' Z" fill="'+f.fill+'" stroke="'+f.stroke+'" stroke-width="1"></path>';
    } else {
      s+='<rect x="'+x1+'" y="'+y+'" width="'+w+'" height="'+h+'" rx="3" fill="'+f.fill+'" stroke="'+f.stroke+'" stroke-width="1"></rect>';
    }
    const cx=x1+w/2;
    if(w>44) s+='<text x="'+cx+'" y="'+(y+h/2+4)+'" text-anchor="middle" font-family="IBM Plex Mono, monospace" font-size="10.5" font-weight="500" fill="'+f.stroke+'">'+esc(f.label)+"</text>";
    else s+='<text x="'+cx+'" y="'+(y-7)+'" text-anchor="middle" font-family="IBM Plex Mono, monospace" font-size="10" fill="'+f.stroke+'">'+esc(f.label)+"</text>";
    s+='<text x="'+cx+'" y="'+(y+h+14)+'" text-anchor="middle" font-family="IBM Plex Mono, monospace" font-size="9.5" fill="var(--ink3)">'+f.len+"</text>";
  }
  return s+"</svg>";
}
function arcPath(cx,cy,r,a1,a2,w){
  const T=(rr,a)=>[cx+rr*Math.cos(a*Math.PI/180), cy+rr*Math.sin(a*Math.PI/180)];
  const ro=r+w/2, ri=r-w/2, large=(a2-a1)>180?1:0;
  const p1=T(ro,a1),p2=T(ro,a2),p3=T(ri,a2),p4=T(ri,a1);
  return "M"+p1[0].toFixed(2)+" "+p1[1].toFixed(2)+"A"+ro+" "+ro+" 0 "+large+" 1 "+p2[0].toFixed(2)+" "+p2[1].toFixed(2)+
         "L"+p3[0].toFixed(2)+" "+p3[1].toFixed(2)+"A"+ri+" "+ri+" 0 "+large+" 0 "+p4[0].toFixed(2)+" "+p4[1].toFixed(2)+"Z";
}
function circularMap(name, total, feats){
  const S=420, cx=S/2, cy=S/2, r=124, w=17;
  let s='<svg class="svgmap" viewBox="0 0 '+S+' '+S+'" role="img" aria-label="Plasmid map of '+esc(name)+'">';
  s+='<circle cx="'+cx+'" cy="'+cy+'" r="'+r+'" fill="none" stroke="var(--line)" stroke-width="'+w+'" opacity=".5"></circle>';
  const step=Math.max(500,Math.round(total/12/500)*500);
  for(let i=0;i<total;i+=step){
    const a=(i/total)*360-90, rad=a*Math.PI/180;
    s+='<line x1="'+(cx+(r+w/2)*Math.cos(rad)).toFixed(1)+'" y1="'+(cy+(r+w/2)*Math.sin(rad)).toFixed(1)+
       '" x2="'+(cx+(r+w/2+4)*Math.cos(rad)).toFixed(1)+'" y2="'+(cy+(r+w/2+4)*Math.sin(rad)).toFixed(1)+
       '" stroke="var(--ink3)" stroke-width="1"></line>';
  }
  const placed=[];
  for(const f of feats){
    const a1=(f.start/total)*360-90, a2=(f.end/total)*360-90;
    s+='<path d="'+arcPath(cx,cy,r,a1,Math.max(a2,a1+0.6),w)+'" fill="'+f.fill+'" stroke="'+f.stroke+'" stroke-width="1"></path>';
    let am=(a1+a2)/2;
    while(placed.some(p=>Math.abs(p-am)<7)) am+=7;
    placed.push(am);
    const rad=am*Math.PI/180;
    const lx=cx+(r+w/2+9)*Math.cos(rad), ly=cy+(r+w/2+9)*Math.sin(rad);
    const right=Math.cos(rad)>=-0.02;
    s+='<text x="'+lx.toFixed(1)+'" y="'+(ly+3.5).toFixed(1)+'" text-anchor="'+(right?"start":"end")+
       '" font-family="IBM Plex Mono, monospace" font-size="10" fill="'+f.stroke+'">'+esc(f.label)+"</text>";
  }
  const short=name.length>22?name.slice(0,21)+"…":name;
  s+='<text x="'+cx+'" y="'+(cy-4)+'" text-anchor="middle" font-family="IBM Plex Sans Condensed, sans-serif" font-size="15" font-weight="600" fill="var(--ink)">'+esc(short)+"</text>";
  s+='<text x="'+cx+'" y="'+(cy+15)+'" text-anchor="middle" font-family="IBM Plex Mono, monospace" font-size="11.5" fill="var(--ink3)">'+bp(total)+"</text>";
  return s+"</svg>";
}
const FC={
  B:{fill:"var(--attB-bg)",stroke:"var(--attB)"}, P:{fill:"var(--attP-bg)",stroke:"var(--attP)"},
  L:{fill:"var(--attL-bg)",stroke:"var(--attL)"}, R:{fill:"var(--attR-bg)",stroke:"var(--attR)"},
  O:{fill:"var(--orf-bg)",stroke:"var(--orf)"},   K:{fill:"var(--rrs-bg)",stroke:"var(--rrs)"},
  N:{fill:"var(--surface2)",stroke:"var(--ink3)"}
};
function feat(label,start,end,c,arrow){return {label,start,end,len:bp(end-start),fill:FC[c].fill,stroke:FC[c].stroke,arrow:!!arrow};}
function featColorFor(label){
  if(/^attB/.test(label))return "B"; if(/^attP/.test(label))return "P";
  if(/^attL/.test(label))return "L"; if(/^attR/.test(label))return "R";
  return "N";
}
