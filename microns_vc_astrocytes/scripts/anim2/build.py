"""Generate the animated soma-vs-endfoot schematic (1920x1080, scene-sequenced)."""
import math, json, os

W, H = 1920, 1080
PX = 9.0                      # px per micrometre
VES_HW, LUM_HW = 20, 13
SR = 45                       # soma radius (10 um cell)
VY0, VY1 = 215, 880
OA, OB = 70, 970              # panel origins
V1, V2 = 130, 690             # vessel x within a panel
SY = 545
GAP = int(9.4 * PX)           # 9.4 um median gap
SXA = OA + V1 + VES_HW + GAP + SR          # typical soma centre
FACEB = OB + V1 + VES_HW                   # vessel face in panel B
SXB = FACEB + 40                           # vc soma centre (clipped at the face)
SXB_START = OB + V1 + VES_HW + GAP + SR    # where it starts before sliding

def ves(cx):
    return (f'<rect class="vw" x="{cx-VES_HW}" y="{VY0}" width="{2*VES_HW}" height="{VY1-VY0}" rx="{VES_HW}"/>'
            f'<rect class="vl" x="{cx-LUM_HW}" y="{VY0+7}" width="{2*LUM_HW}" height="{VY1-VY0-14}" rx="{LUM_HW}"/>')

def cap(fx, y0, y1, bulge, side):
    b = bulge * side; ym = (y0 + y1) / 2
    return (f'M{fx},{y0} C{fx+b*0.9},{y0+5} {fx+b},{y0+(ym-y0)*0.45} {fx+b},{ym:.0f} '
            f'C{fx+b},{y1-(y1-ym)*0.45} {fx+b*0.9},{y1-5} {fx},{y1} Z')

def twigs(cx, cy, angles, lo=120, hi=185):
    prim, sec = [], []
    for i, a in enumerate(angles):
        r = math.radians(a); j = ((i * 37) % 11 - 5) / 5.0
        L = lo + (hi - lo) * (((i * 53) % 13) / 12.0)
        x0, y0 = cx + math.cos(r) * (SR - 3), cy + math.sin(r) * (SR - 3)
        bd = math.radians(a + j * 13)
        xm, ym = cx + math.cos(bd) * (SR + L * .55), cy + math.sin(bd) * (SR + L * .55)
        x1, y1 = cx + math.cos(r) * (SR + L), cy + math.sin(r) * (SR + L)
        prim.append(f'M{x0:.0f},{y0:.0f} Q{xm:.0f},{ym:.0f} {x1:.0f},{y1:.0f}')
        for k, da in enumerate((-26, 24)):
            r2 = math.radians(a + da + j * 5); L2 = 32 + ((i * 29 + k * 17) % 28)
            sec.append(f'M{x1:.0f},{y1:.0f} Q{x1+math.cos(r2)*L2*.6:.0f},{y1+math.sin(r2)*L2*.6:.0f} '
                       f'{x1+math.cos(r2)*L2:.0f},{y1+math.sin(r2)*L2:.0f}')
    return prim, sec

def proc(cx, cy, tx, ty):
    d = math.hypot(tx - cx, ty - cy); ux, uy = (tx - cx) / d, (ty - cy) / d
    x0, y0 = cx + ux * (SR - 5), cy + uy * (SR - 5)
    mx, my = (x0 + tx) / 2 - uy * 34, (y0 + ty) / 2 + ux * 34
    return f'M{x0:.0f},{y0:.0f} Q{mx:.0f},{my:.0f} {tx:.0f},{ty:.0f}'

S = []
A = S.append
A(f'<svg id="stage" viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg">')
A('<defs>'
  f'<clipPath id="clipB"><rect x="{FACEB}" y="0" width="{W}" height="{H}"/></clipPath>'
  '<clipPath id="magA"><circle cx="520" cy="330" r="150"/></clipPath>'
  '<clipPath id="magB"><circle cx="1452" cy="330" r="150"/></clipPath>'
  '<marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
  'orient="auto-start-reverse"><path d="M0,1 L9,5 L0,9 z" fill="currentColor" opacity=".6"/></marker>'
  '</defs>')

# ---------------- panel headers ----------------
A(f'<g id="hdrA" opacity="0"><text class="kick kA" x="{OA}" y="192">Typical astrocyte</text></g>')
A(f'<g id="hdrB" opacity="0"><text class="kick kB" x="{OB}" y="192">Vascular-coupled astrocyte</text></g>')
A(f'<line id="divider" x1="960" y1="190" x2="960" y2="905" class="dv" opacity="0"/>')

# ---------------- vessels ----------------
A(f'<g id="vesA" opacity="0">{ves(OA+V1)}{ves(OA+V2)}</g>')
A(f'<g id="vesB" opacity="0">{ves(OB+V1)}{ves(OB+V2)}</g>')

# ---------------- panel A cell ----------------
tp, ts = twigs(SXA, SY, [15, 48, 80, 112, 145, 170, 198, 232, 265, 298, 330, 352])
A('<g id="twA" opacity="0"><g class="tw2">' + ''.join(f'<path d="{d}"/>' for d in tp) +
  '</g><g class="tw1">' + ''.join(f'<path d="{d}"/>' for d in ts) + '</g></g>')
A('<g id="prA" opacity="0">'
  f'<path class="pr" d="{proc(SXA,SY,OA+V1+VES_HW+12,370)}"/>'
  f'<path class="pr" d="{proc(SXA,SY,OA+V2-VES_HW-12,400)}"/>'
  f'<path class="pr" d="{proc(SXA,SY,OA+V2-VES_HW-12,740)}"/></g>')
A('<g id="ftA" opacity="0">'
  f'<path class="ft" d="{cap(OA+V1+VES_HW,330,410,17,1)}"/>'
  f'<path class="ft" d="{cap(OA+V2-VES_HW,360,440,17,-1)}"/>'
  f'<path class="ft" d="{cap(OA+V2-VES_HW,700,780,17,-1)}"/></g>')
A(f'<g id="smA" opacity="0"><circle class="sm" cx="{SXA}" cy="{SY}" r="{SR}"/>'
  f'<circle class="nu" cx="{SXA}" cy="{SY}" r="24"/></g>')
# gap dimension
A(f'<g id="gapA" opacity="0"><line class="dim" x1="{OA+V1+VES_HW}" y1="690" x2="{SXA-SR}" y2="690" '
  'marker-start="url(#ar)" marker-end="url(#ar)"/>'
  f'<line class="dim" x1="{OA+V1+VES_HW}" y1="660" x2="{OA+V1+VES_HW}" y2="702"/>'
  f'<line class="dim" x1="{SXA-SR}" y1="600" x2="{SXA-SR}" y2="702"/>'
  f'<text class="lb" x="{(OA+V1+VES_HW+SXA-SR)//2}" y="680" text-anchor="middle">9.4 µm</text></g>')

# ---------------- panel B cell ----------------
tpB, tsB = twigs(SXB, SY, [-68, -40, -14, 14, 40, 68, 92, 268, 292, 320])
A('<g id="cellB" opacity="0">')
A('<g id="twB"><g class="tw2">' + ''.join(f'<path d="{d}"/>' for d in tpB) +
  '</g><g class="tw1">' + ''.join(f'<path d="{d}"/>' for d in tsB) + '</g></g>')
A('<g id="prB">'
  f'<path class="pr" d="{proc(SXB,SY,OB+V2-VES_HW-12,400)}"/>'
  f'<path class="pr" d="{proc(SXB,SY,OB+V2-VES_HW-12,740)}"/></g>')
A('</g>')
A('<g id="ftB" opacity="0">'
  f'<path class="ft" d="{cap(OB+V2-VES_HW,360,440,17,-1)}"/>'
  f'<path class="ft" d="{cap(OB+V2-VES_HW,700,780,17,-1)}"/></g>')
A(f'<g id="shB" opacity="0"><path id="shPath" class="sh" d="{cap(FACEB,390,700,24,1)}"/></g>')
A(f'<g id="smB" opacity="0" clip-path="url(#clipB)"><circle class="sm" cx="{SXB}" cy="{SY}" r="{SR}"/>'
  f'<circle class="nu" cx="{SXB+6}" cy="{SY}" r="24"/></g>')

# ---------------- callout labels ----------------
A(f'<g id="labFt" opacity="0"><line class="dim" x1="{OA+V1+VES_HW+18}" y1="352" x2="{OA+V1+VES_HW+70}" y2="300"/>'
  f'<text class="lb kA" x="{OA+V1+VES_HW+76}" y="294">endfoot</text>'
  f'<text class="lbs" x="{OA+V1+VES_HW+76}" y="316">88 µm² — one patch</text></g>')
A(f'<g id="labSh" opacity="0"><line class="dim" x1="{FACEB+24}" y1="420" x2="{FACEB+86}" y2="352"/>'
  f'<text class="lb kB" x="{FACEB+92}" y="346">somatic sheath</text>'
  f'<text class="lbs" x="{FACEB+92}" y="368">525 µm² — about 10× an endfoot</text></g>')
A(f'<g id="labFtB" opacity="0"><line class="dim" x1="{OB+V2+VES_HW+10}" y1="396" x2="{OB+V2+VES_HW+44}" y2="352"/>'
  f'<text class="lb kA" x="{OB+V2+VES_HW+50}" y="346">endfeet kept</text>'
  f'<text class="lbs" x="{OB+V2+VES_HW+50}" y="368">86% on a different vessel</text></g>')

# ---------------- magnifiers ----------------
def mag(idn, cx, cy, clip, body, labels):
    return (f'<g id="{idn}" opacity="0"><circle class="magbg" cx="{cx}" cy="{cy}" r="150"/>'
            f'<g clip-path="url(#{clip})">{body}</g>'
            f'<circle class="magring" cx="{cx}" cy="{cy}" r="150"/>{labels}</g>')

# -- endfoot magnifier: vessel wall at bottom, thin endfoot layer, mitochondria + ribosomes
eb = []
eb.append('<rect class="mag-lum" x="370" y="408" width="300" height="80"/>')
eb.append('<rect class="mag-endo" x="370" y="386" width="300" height="24"/>')
eb.append('<rect class="mag-bl" x="370" y="376" width="300" height="11"/>')
eb.append('<rect class="mag-foot" x="370" y="330" width="300" height="47"/>')
eb.append('<rect class="mag-neuro" x="370" y="180" width="300" height="151"/>')
for i, (mx, my, rot) in enumerate([(425, 353, -6), (505, 350, 4), (585, 355, -3)]):
    eb.append(f'<g class="mito" transform="translate({mx},{my}) rotate({rot})">'
              '<rect x="-27" y="-13" width="54" height="26" rx="13"/>'
              '<path class="crista" d="M-16,-11 L-16,11 M-6,-11 L-6,11 M4,-11 L4,11 M14,-11 L14,11"/></g>')
for i in range(16):
    rx = 388 + (i * 37) % 268; ry = 338 + ((i * 53) % 3) * 9
    eb.append(f'<circle class="ribo" cx="{rx}" cy="{ry}" r="3.6"/>')
el = ('<text class="maglb" x="520" y="206" text-anchor="middle">ENDFOOT</text>'
      '<text class="maglb2" x="520" y="286" text-anchor="middle">electron transport chain</text>'
      '<text class="maglb2" x="520" y="314" text-anchor="middle">and its own ribosomes</text>'
      '<text class="maglb3" x="520" y="403" text-anchor="middle">vessel wall</text>')
A('<g id="leadA" opacity="0"><line class="dim" x1="371" y1="350" x2="232" y2="368"/></g>')
A(mag('magAg', 520, 330, 'magA', ''.join(eb), el))

# -- soma magnifier: vessel at left, soma cytoplasm with nucleus, Golgi, glycogen
sb = []
sb.append('<rect class="mag-lum" x="1302" y="180" width="62" height="300"/>')
sb.append('<rect class="mag-endo" x="1364" y="180" width="22" height="300"/>')
sb.append('<rect class="mag-bl" x="1386" y="180" width="10" height="300"/>')
sb.append('<rect class="mag-soma" x="1396" y="180" width="210" height="300"/>')
sb.append('<circle class="mag-nuc" cx="1510" cy="330" r="62"/>')
sb.append('<circle class="mag-nucl" cx="1522" cy="318" r="17"/>')
sb.append('<g class="golgi"><path d="M1416,246 q26,-16 52,-4 M1414,258 q26,-16 52,-4 M1412,270 q26,-16 52,-4"/></g>')
sb.append('<g class="er"><path d="M1410,392 q22,12 44,0 q22,-12 44,0 M1410,408 q22,12 44,0 q22,-12 44,0"/></g>')
for gx, gy in [(1430, 340), (1448, 356), (1424, 364), (1442, 330), (1458, 344)]:
    for k in range(5):
        a = k * 72 * math.pi / 180
        sb.append(f'<circle class="glyc" cx="{gx+math.cos(a)*6:.0f}" cy="{gy+math.sin(a)*6:.0f}" r="2.6"/>')
sb.append('<g class="mito" transform="translate(1560,420) rotate(70)">'
          '<rect x="-24" y="-12" width="48" height="24" rx="12"/>'
          '<path class="crista" d="M-14,-10 L-14,10 M-5,-10 L-5,10 M4,-10 L4,10 M13,-10 L13,10"/></g>')
sl = ('<text class="maglb" x="1452" y="216" text-anchor="middle">SOMA CONTACT</text>'
      '<text class="maglb2" x="1510" y="334" text-anchor="middle">nucleus</text>'
      '<text class="maglb3" x="1444" y="232" text-anchor="start">Golgi</text>'
      '<text class="maglb3" x="1398" y="438" text-anchor="start">ER</text>'
      '<text class="maglb3" x="1414" y="392" text-anchor="start">glycogen</text>')
A('<g id="leadB" opacity="0"><line class="dim" x1="1326" y1="412" x2="1135" y2="528"/></g>')
A(mag('magBg', 1452, 330, 'magB', ''.join(sb), sl))

# -- disease overlays inside the magnifiers
A('<g id="dzA" opacity="0" clip-path="url(#magA)">')
for i in range(9):
    x = 392 + i * 32
    A(f'<rect class="aqp" id="aqp{i}" x="{x}" y="368" width="15" height="9" rx="2.5"/>')
A('</g>')
A('<g id="dzB" opacity="0" clip-path="url(#magB)">'
  '<circle class="mag-nuc dz" cx="1480" cy="296" r="46"/><circle class="mag-nucl dz" cx="1490" cy="286" r="12"/>'
  '<circle class="mag-nuc dz" cx="1548" cy="390" r="46"/><circle class="mag-nucl dz" cx="1558" cy="380" r="12"/>'
  '<path class="spindle" d="M1498,330 L1530,356"/>'
  '<text class="maglb2" x="1452" y="452" text-anchor="middle">the cell divides</text></g>')

# ---------------- stat bars ----------------
A('<g id="stats" opacity="0">'
  f'<text class="statlb" x="{OA}" y="952">Vessel-facing membrane, whole cell</text>'
  f'<text class="statlb" x="{OB}" y="952">Vessel-facing membrane, whole cell</text>'
  f'<rect class="barbg" x="{OA}" y="968" width="560" height="26" rx="13"/>'
  f'<rect class="barbg" x="{OB}" y="968" width="560" height="26" rx="13"/>'
  f'<rect id="barA" class="bar bA" x="{OA}" y="968" width="0" height="26" rx="13"/>'
  f'<rect id="barB" class="bar bB" x="{OB}" y="968" width="0" height="26" rx="13"/>'
  f'<text id="numA" class="statnum kA" x="{OA+576}" y="990">224 µm²</text>'
  f'<text id="numB" class="statnum kB" x="{OB+576}" y="990">715 µm²</text>'
  f'<text class="statsub" x="{OA}" y="1026">0.53% of the cell surface · median of 116 cells</text>'
  f'<text class="statsub" x="{OB}" y="1026">1.70% · the extra is almost all cell body · 82 cells</text>'
  '</g>')
A('</svg>')

svg = '\n'.join(S)

page = f'''<!doctype html><html><head><meta charset="utf-8">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600&family=Newsreader:opsz,wght@6..72,400;6..72,600;6..72,700&display=swap">
<style>
:root{{--bg:#0e1116;--panel:#161a21;--fg:#eef1f5;--muted:#99a3b0;--rule:#2a3039;
--vessel:#d1434d;--lumen:#3b2327;--proc:#6fc9de;--soma:#a98cf5;--contact:#ff8a3d;--foot:#f5c542;
--endo:#8e2730;--bl:#c9b06a;--neuro:#1b2026;}}
html,body{{margin:0;width:{W}px;height:{H}px;overflow:hidden;background:var(--bg);color:var(--fg);
font-family:"IBM Plex Sans",system-ui,sans-serif}}
svg{{position:absolute;left:0;top:0;width:{W}px;height:{H}px}}
.vw{{fill:var(--vessel)}} .vl{{fill:var(--lumen)}}
.sm{{fill:var(--soma)}} .nu{{fill:#fff;opacity:.3}}
.ft{{fill:var(--foot)}} .sh{{fill:var(--contact)}}
.pr{{fill:none;stroke:var(--proc);stroke-width:11;stroke-linecap:round}}
.tw2 path{{fill:none;stroke:var(--proc);stroke-width:3;opacity:.6;stroke-linecap:round}}
.tw1 path{{fill:none;stroke:var(--proc);stroke-width:1.8;opacity:.4;stroke-linecap:round}}
.dv{{stroke:var(--rule);stroke-width:2}}
.dim{{stroke:var(--fg);stroke-width:1.6;opacity:.55;stroke-dasharray:5 4}}
.kick{{font:600 19px/1 "IBM Plex Mono",monospace;letter-spacing:.13em;text-transform:uppercase}}
.kA{{fill:var(--foot)}} .kB{{fill:var(--contact)}}
.lb{{fill:var(--fg);font-size:21px;font-weight:500;paint-order:stroke;stroke:var(--bg);stroke-width:5px;stroke-linejoin:round}}
.lbs{{fill:var(--muted);font-size:17px;paint-order:stroke;stroke:var(--bg);stroke-width:5px;stroke-linejoin:round}}
.magbg{{fill:var(--panel);stroke:none}}
.magring{{fill:none;stroke:var(--rule);stroke-width:3}}
.mag-lum{{fill:var(--lumen)}} .mag-endo{{fill:var(--endo)}} .mag-bl{{fill:var(--bl)}}
.mag-foot{{fill:var(--foot);opacity:.5}} .mag-soma{{fill:var(--soma);opacity:.42}}
.mag-neuro{{fill:var(--neuro)}}
.mag-nuc{{fill:#fff;opacity:.26}} .mag-nucl{{fill:#fff;opacity:.4}}
.mito rect{{fill:#2f3a44;stroke:#7f8c99;stroke-width:2}}
.crista{{stroke:#7f8c99;stroke-width:1.6;fill:none}}
.ribo{{fill:#4fd1c5}}
.golgi path{{fill:none;stroke:#ffd98a;stroke-width:4;stroke-linecap:round}}
.er path{{fill:none;stroke:#b9c6d4;stroke-width:3;opacity:.8}}
.glyc{{fill:#8fe3b0}}
.aqp{{fill:#6fc9de;stroke:#0e1116;stroke-width:1.5}}
.spindle{{stroke:#fff;stroke-width:3;opacity:.6;stroke-dasharray:5 5}}
.maglb{{fill:var(--fg);font:600 17px/1 "IBM Plex Mono",monospace;letter-spacing:.1em}}
.maglb2{{fill:var(--fg);font-size:17px;paint-order:stroke;stroke:var(--panel);stroke-width:4px}}
.maglb3{{fill:var(--muted);font-size:14px;paint-order:stroke;stroke:var(--panel);stroke-width:4px}}
.statlb{{fill:var(--muted);font-size:17px}}
.statsub{{fill:var(--muted);font-size:16px}}
.statnum{{font:700 28px/1 "Newsreader",Georgia,serif}}
.barbg{{fill:var(--panel)}} .bA{{fill:var(--foot)}} .bB{{fill:var(--contact)}}
#cap{{position:absolute;left:70px;top:42px;width:1780px}}
#cap .ct{{font:700 44px/1.15 "Newsreader",Georgia,serif;margin:0;letter-spacing:-.01em}}
#cap .cs{{font:400 23px/1.4 "IBM Plex Sans",sans-serif;color:var(--muted);margin:10px 0 0}}
#src{{position:absolute;left:70px;bottom:14px;font:400 14px/1.4 "IBM Plex Mono",monospace;color:#5d6773}}
</style></head><body>
{svg}
<div id="cap"><p class="ct" id="ctitle"></p><p class="cs" id="csub"></p></div>
<div id="src">MICrONS minnie65 · geometry from 240 reconstructed astrocytes · compartment contents after Stokum 2021 and Boulay 2017 · schematic, not to scale inside the magnifiers</div>
<script>
const $=id=>document.getElementById(id);
const clamp=(x,a,b)=>Math.min(b,Math.max(a,x));
const ramp=(t,a,b)=>clamp((t-a)/(b-a),0,1);
const ease=x=>x<.5?2*x*x:1-Math.pow(-2*x+2,2)/2;
const show=(id,v)=>{{const e=$(id); if(e) e.setAttribute('opacity',v.toFixed(3));}};
const SXB_START={SXB_START}, SXB={SXB}, FACEB={FACEB}, SY={SY};

// measure path lengths once so processes can draw themselves on
const PATHS=[...document.querySelectorAll('#prA path,#prB path,.tw2 path,.tw1 path')];
PATHS.forEach(p=>{{const L=p.getTotalLength(); p.dataset.len=L; p.style.strokeDasharray=L; p.style.strokeDashoffset=L;}});
function drawPaths(sel,f){{document.querySelectorAll(sel).forEach((p,i)=>{{
  const L=+p.dataset.len, d=clamp((f*1.35)-(i%7)*0.05,0,1);
  p.style.strokeDashoffset=(L*(1-ease(d))).toFixed(1);}});}}

const CAPS=[
 [0.0,'Every cortical astrocyte touches a blood vessel','Most of them reach it the classical way — with thin endfeet at the ends of long processes'],
 [6.5,'The endfoot is a small, remote patch','88 µm² of membrane, about 23 µm from the nucleus, at the end of a process'],
 [12.0,'In 5.6% of astrocytes the cell body itself lies on the capillary','2.5× more often than chance placement predicts — a real population, not an accident of packing'],
 [18.5,'The sheath is roughly ten times the area of one endfoot','And the cell keeps its endfeet: same number as a typical astrocyte, 86% of them on a different vessel'],
 [24.0,'Same cell, three times the vessel-facing membrane','Total cell surface is unchanged. Nearly all the extra contact is at the cell body'],
 [28.5,'The two interfaces are built from different parts','An endfoot carries the electron transport chain and its own ribosomes — it is metabolically and translationally self-supplying'],
 [34.0,'A soma contact puts the nucleus on the basal lamina','With the bulk Golgi, the ER and the glycogen store — under 1 µm from blood instead of 23'],
 [39.5,'So they fail in different ways','The endfoot loses a surface: AQP4 leaves the membrane and water follows. The cell body changes what it is: at the vessel, astrocytes divide'],
 [45.5,'What does a cell body on a vessel mean?','A short line from blood to genome, a secretory face on the perivascular space, and a cell positioned to proliferate']
];
function caption(t){{let c=CAPS[0]; for(const k of CAPS) if(t>=k[0]) c=k;
  $('ctitle').textContent=c[1]; $('csub').textContent=c[2];
  const i=CAPS.indexOf(c), fade=ramp(t,c[0],c[0]+0.45);
  $('cap').style.opacity=fade.toFixed(3);}}

window.renderAt=function(t){{
  caption(t);
  show('hdrA',ramp(t,1.2,2.0)); show('vesA',ramp(t,0.3,1.4));
  show('divider',ramp(t,11.0,12.0));
  show('hdrB',ramp(t,11.4,12.2)); show('vesB',ramp(t,11.0,12.0));

  // panel A build
  show('smA',ramp(t,2.0,2.8));
  const gA=ramp(t,2.6,5.2); show('twA',ramp(t,2.6,3.2)); show('prA',ramp(t,2.6,3.2));
  drawPaths('#prA path',gA); drawPaths('.tw2 path',ramp(t,2.6,4.2)); drawPaths('.tw1 path',ramp(t,3.6,5.2));
  show('ftA',ramp(t,5.6,6.6));
  show('labFt',ramp(t,7.0,7.8)*(1-ramp(t,27.5,28.3)));
  show('gapA',ramp(t,8.4,9.2)*(1-ramp(t,27.5,28.3)));

  // panel B build, soma slides onto the vessel
  show('cellB',ramp(t,12.4,13.2)); show('ftB',ramp(t,12.4,13.2));
  drawPaths('#prB path',ramp(t,12.4,14.0));
  show('smB',ramp(t,12.4,13.2));
  const slide=ease(ramp(t,14.2,16.4));
  const cx=SXB_START+(SXB-SXB_START)*slide;
  $('smB').setAttribute('transform','translate('+(cx-SXB).toFixed(1)+',0)');
  $('cellB').setAttribute('transform','translate('+(cx-SXB).toFixed(1)+',0)');
  // sheath grows along the vessel from the contact point outward
  const g=ease(ramp(t,16.4,18.6));
  const y0=SY-(SY-390)*g, y1=SY+(700-SY)*g, b=24*Math.min(1,g*2.2);
  const ym=(y0+y1)/2;
  $('shPath').setAttribute('d','M'+FACEB+','+y0.toFixed(0)+' C'+(FACEB+b*0.9)+','+(y0+5).toFixed(0)+' '+(FACEB+b)+','+(y0+(ym-y0)*0.45).toFixed(0)+' '+(FACEB+b)+','+ym.toFixed(0)+' C'+(FACEB+b)+','+(y1-(y1-ym)*0.45).toFixed(0)+' '+(FACEB+b*0.9)+','+(y1-5).toFixed(0)+' '+FACEB+','+y1.toFixed(0)+' Z');
  show('shB',ramp(t,16.4,16.9));
  show('labSh',ramp(t,19.0,19.8)*(1-ramp(t,27.5,28.3)));
  show('labFtB',ramp(t,20.6,21.4)*(1-ramp(t,27.5,28.3)));

  // stat bars
  show('stats',ramp(t,24.2,25.0));
  const bg=ease(ramp(t,24.6,26.4));
  $('barA').setAttribute('width',(560*(224/715)*bg).toFixed(1));
  $('barB').setAttribute('width',(560*bg).toFixed(1));
  $('numA').setAttribute('opacity',ramp(t,25.6,26.2).toFixed(2));
  $('numB').setAttribute('opacity',ramp(t,25.6,26.2).toFixed(2));

  // magnifiers
  const mA=ramp(t,28.6,29.6)*(1-ramp(t,45.2,46.0));
  const mB=ramp(t,34.2,35.2)*(1-ramp(t,45.2,46.0));
  show('magAg',mA); show('magBg',mB); show('leadA',mA); show('leadB',mB);
  $('magAg').setAttribute('transform','translate(520,330) scale('+(0.25+0.75*ease(mA)).toFixed(3)+') translate(-520,-330)');
  $('magBg').setAttribute('transform','translate(1452,330) scale('+(0.25+0.75*ease(mB)).toFixed(3)+') translate(-1452,-330)');

  // disease
  const d1=ramp(t,40.0,43.0);
  show('dzA',mA*ramp(t,39.8,40.4));
  for(let i=0;i<9;i++){{const e=$('aqp'+i); if(!e) continue;
    const ph=clamp(d1*1.5-(i%5)*0.12,0,1);
    e.setAttribute('transform','translate('+((i%2?1:-1)*ph*26).toFixed(1)+','+(-ph*118).toFixed(1)+')');
    e.setAttribute('opacity',(1-0.45*ph).toFixed(2));}}
  show('dzB',mB*ramp(t,41.0,42.2));
  return true;
}};
window.READY=true;
</script></body></html>'''

os.makedirs(os.path.dirname(os.path.abspath(__file__)), exist_ok=True)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'anim.html'), 'w').write(page)
print('written', len(page), 'chars; SXB_START', SXB_START, 'SXB', SXB, 'FACEB', FACEB)
