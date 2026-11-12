"""Generate the hero schematic SVG geometrically so processes, endfeet and
labels line up exactly. Writes the <svg> block and splices it into the page."""
import math, re

PX = 7.0                      # px per micrometre
VW, VH = 1200, 700
SOMA_R = 35                   # 10 um soma
VES_HW, LUM_HW = 15, 9        # vessel outer / lumen half-widths
VY0, VY1 = 74, 600

def ves(cx):
    return (f'<rect class="vessel-wall" x="{cx-VES_HW}" y="{VY0}" width="{2*VES_HW}" '
            f'height="{VH-VY0-100}" rx="{VES_HW}"/>'
            f'<rect class="vessel-lumen" x="{cx-LUM_HW}" y="{VY0+6}" width="{2*LUM_HW}" '
            f'height="{VH-VY0-112}" rx="{LUM_HW}"/>')

def cap(face_x, y0, y1, bulge, side):
    """Membrane patch hugging a vessel face. side=+1 bulges right."""
    b = bulge * side
    ym = (y0 + y1) / 2
    return (f'M{face_x},{y0} C{face_x+b*0.9},{y0+4} {face_x+b},{y0+(ym-y0)*0.45} {face_x+b},{ym} '
            f'C{face_x+b},{y1-(y1-ym)*0.45} {face_x+b*0.9},{y1-4} {face_x},{y1} Z')

def twigs(cx, cy, angles, rng=(96, 150), seed=0):
    """Fine process bush: primary twigs from the soma edge plus secondaries."""
    prim, sec = [], []
    for i, a in enumerate(angles):
        r = math.radians(a)
        jitter = ((i * 37) % 11 - 5) / 5.0
        L = rng[0] + (rng[1] - rng[0]) * (((i * 53) % 13) / 12.0)
        x0, y0 = cx + math.cos(r) * (SOMA_R - 3), cy + math.sin(r) * (SOMA_R - 3)
        bend = math.radians(a + jitter * 13)
        xm, ym = cx + math.cos(bend) * (SOMA_R + L * 0.55), cy + math.sin(bend) * (SOMA_R + L * 0.55)
        x1, y1 = cx + math.cos(r) * (SOMA_R + L), cy + math.sin(r) * (SOMA_R + L)
        prim.append(f'M{x0:.0f},{y0:.0f} Q{xm:.0f},{ym:.0f} {x1:.0f},{y1:.0f}')
        for k, da in enumerate((-26, 24)):
            r2 = math.radians(a + da + jitter * 5)
            L2 = 26 + ((i * 29 + k * 17) % 22)
            sec.append(f'M{x1:.0f},{y1:.0f} Q{x1+math.cos(r2)*L2*0.6:.0f},{y1+math.sin(r2)*L2*0.6:.0f} '
                       f'{x1+math.cos(r2)*L2:.0f},{y1+math.sin(r2)*L2:.0f}')
    return prim, sec

def proc(cx, cy, tx, ty, w):
    """Thick process from soma edge to an endfoot."""
    d = math.hypot(tx - cx, ty - cy)
    ux, uy = (tx - cx) / d, (ty - cy) / d
    x0, y0 = cx + ux * (SOMA_R - 4), cy + uy * (SOMA_R - 4)
    mx, my = (x0 + tx) / 2 - uy * 26, (y0 + ty) / 2 + ux * 26
    return f'<path class="proc" stroke-width="{w}" d="M{x0:.0f},{y0:.0f} Q{mx:.0f},{my:.0f} {tx:.0f},{ty:.0f}"/>'

out = []
A = out.append
A(f'<svg class="d" viewBox="0 0 {VW} {VH}" role="img" aria-label="Side-by-side schematic of two '
  'astrocytes. In the typical astrocyte the cell body sits about nine micrometres from a capillary and '
  'reaches it with a thin endfoot. In the vascular-coupled astrocyte the cell body itself lies along the '
  'capillary as a large sheath, while the same number of endfeet still reach a second, separate vessel.">')
A('<defs><marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
  'orient="auto-start-reverse"><path d="M0,1 L9,5 L0,9 z" fill="currentColor" opacity=".6"/></marker></defs>')

# ---------------- PANEL A : typical ----------------
V1, V2 = 95, 420
SX, SY = V1 + VES_HW + int(9.4 * PX) + SOMA_R, 345     # 9.4 um gap, nucleus to lumen edge
A(f'<text class="lbl-k" x="22" y="26" fill="var(--endfoot)">Typical astrocyte</text>')
A(f'<text class="ttl" x="22" y="52">Reaches vessels with endfeet</text>')
A(ves(V1)); A(ves(V2))
ang = [15, 48, 80, 112, 145, 170, 198, 232, 265, 298, 330, 352]
p, s = twigs(SX, SY, ang)
A('<g class="twig" stroke-width="2.3">' + ''.join(f'<path d="{d}"/>' for d in p) + '</g>')
A('<g class="twig" stroke-width="1.4">' + ''.join(f'<path d="{d}"/>' for d in s) + '</g>')
A(proc(SX, SY, V1 + VES_HW + 10, 196, 8))
A(proc(SX, SY, V2 - VES_HW - 10, 252, 7))
A(proc(SX, SY, V2 - VES_HW - 10, 452, 8))
A(f'<path class="foot" d="{cap(V1+VES_HW,166,226,13,1)}"/>')
A(f'<path class="foot" d="{cap(V2-VES_HW,222,282,13,-1)}"/>')
A(f'<path class="foot" d="{cap(V2-VES_HW,422,482,13,-1)}"/>')
A(f'<circle class="soma-body" cx="{SX}" cy="{SY}" r="{SOMA_R}"/>')
A(f'<circle class="nucleus" cx="{SX}" cy="{SY}" r="19"/>')
# gap dimension, below the soma; labels carry a halo so they survive the twigs
A(f'<line class="dim" x1="{V1+VES_HW}" y1="446" x2="{SX-SOMA_R}" y2="446" '
  'marker-start="url(#ar)" marker-end="url(#ar)"/>')
A(f'<line class="dim" x1="{V1+VES_HW}" y1="424" x2="{V1+VES_HW}" y2="456"/>')
A(f'<line class="dim" x1="{SX-SOMA_R}" y1="386" x2="{SX-SOMA_R}" y2="456"/>')
A(f'<text class="lbl" x="{(V1+VES_HW+SX-SOMA_R)//2}" y="438" text-anchor="middle">9.4 µm</text>')
A(f'<text class="lbl-s" x="{SX-SOMA_R+8}" y="470">median nucleus-to-lumen gap,</text>')
A(f'<text class="lbl-s" x="{SX-SOMA_R+8}" y="484">astrocytes overall</text>')
A(f'<line class="dim" x1="{V1+VES_HW+12}" y1="182" x2="{V1+VES_HW+44}" y2="138"/>')
A(f'<text class="lbl" x="{V1+VES_HW+48}" y="128" fill="var(--endfoot)">endfoot</text>')
A(f'<text class="lbl-s" x="{V1+VES_HW+48}" y="144">88 µm² median — one patch</text>')
A('<text class="lbl-s" x="22" y="640">Vessel-facing membrane, whole cell</text>')
A('<text class="ttl" x="22" y="670" fill="var(--endfoot)">224 µm²</text>')
A('<text class="lbl-s" x="128" y="670">0.53% of the cell surface</text>')

A(f'<line x1="590" y1="18" x2="590" y2="604" class="axis"/>')

# ---------------- PANEL B : vascular-coupled ----------------
O = 590
V1b, V2b = V1 + O, 1010
SXb, SYb = V1b + VES_HW + 30, 345                       # soma centre, pressed on the vessel
A(f'<text class="lbl-k" x="{O+22}" y="26" fill="var(--contact)">Vascular-coupled astrocyte</text>')
A(f'<text class="ttl" x="{O+22}" y="52">Adds a cell-body sheath, keeps the endfeet</text>')
A(ves(V1b)); A(ves(V2b))
# the sheath: a long collar on the vessel, drawn before the soma so it reads as continuous
A(f'<path class="sheath" d="{cap(V1b+VES_HW,228,462,19,1)}"/>')
angb = [-68, -40, -14, 14, 40, 68, 92, 268, 292, 320]
p, s = twigs(SXb, SYb, angb)
A('<g class="twig" stroke-width="2.3">' + ''.join(f'<path d="{d}"/>' for d in p) + '</g>')
A('<g class="twig" stroke-width="1.4">' + ''.join(f'<path d="{d}"/>' for d in s) + '</g>')
A(proc(SXb, SYb, V2b - VES_HW - 10, 252, 7))
A(proc(SXb, SYb, V2b - VES_HW - 10, 452, 8))
A(f'<path class="foot" d="{cap(V2b-VES_HW,222,282,13,-1)}"/>')
A(f'<path class="foot" d="{cap(V2b-VES_HW,422,482,13,-1)}"/>')
# soma: a full-size circle with its left face flattened against the vessel
fx = V1b + VES_HW
CXb = fx + 30                                  # centre 30 px right of the face
half = math.sqrt(SOMA_R**2 - 30**2)            # half-chord where the circle meets the face
A(f'<path class="soma-body" d="M{fx},{SYb-half:.0f} '
  f'A{SOMA_R},{SOMA_R} 0 1 1 {fx},{SYb+half:.0f} Z"/>')
A(f'<ellipse class="nucleus" cx="{CXb}" cy="{SYb}" rx="18" ry="19"/>')
A(f'<line class="dim" x1="{fx}" y1="500" x2="{fx+16}" y2="516"/>')
A(f'<text class="lbl" x="{fx+20}" y="520">&lt; 1 µm</text>')
A(f'<text class="lbl-s" x="{fx+20}" y="536">nucleus to lumen</text>')
A(f'<line class="dim" x1="{fx+19}" y1="246" x2="{fx+54}" y2="206"/>')
A(f'<text class="lbl" x="{fx+58}" y="196" fill="var(--contact)">somatic sheath</text>')
A(f'<text class="lbl-s" x="{fx+58}" y="212">525 µm² median — about 10×</text>')
A(f'<text class="lbl-s" x="{fx+58}" y="226">the area of one endfoot</text>')
A(f'<line class="dim" x1="{V2b-VES_HW-10}" y1="246" x2="{V2b+VES_HW+10}" y2="214"/>')
A(f'<text class="lbl" x="{V2b+VES_HW+14}" y="206" fill="var(--endfoot)">endfeet</text>')
A(f'<text class="lbl-s" x="{V2b+VES_HW+14}" y="222">same number as a typical</text>')
A(f'<text class="lbl-s" x="{V2b+VES_HW+14}" y="236">cell — and 86% are on a</text>')
A(f'<text class="lbl-s" x="{V2b+VES_HW+14}" y="250">different vessel</text>')
A(f'<text class="lbl-s" x="{O+22}" y="640">Vessel-facing membrane, whole cell</text>')
A(f'<text class="ttl" x="{O+22}" y="670" fill="var(--contact)">715 µm²</text>')
A(f'<text class="lbl-s" x="{O+128}" y="670">1.70% — the extra is almost all cell body</text>')
# scale bar
A(f'<line x1="{V2b+VES_HW+14}" y1="588" x2="{V2b+VES_HW+14+10*PX:.0f}" y2="588" stroke="currentColor" stroke-width="3"/>')
A(f'<text class="lbl-s" x="{V2b+VES_HW+14}" y="606">10 µm</text>')
A('</svg>')

svg = '\n'.join(out)
page = open('soma_vs_endfoot.html').read()
new = re.sub(r'<svg class="d" viewBox="0 0 1200 \d+".*?</svg>', svg, page, count=1, flags=re.S)
assert new != page, 'hero svg not replaced'
open('soma_vs_endfoot.html', 'w').write(new)
print('replaced; svg chars', len(svg))
