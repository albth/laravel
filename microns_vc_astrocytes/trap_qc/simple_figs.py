"""Two simplified, presentation-grade figures for the vcAst TRAP QC.

Figure 1  "The two pull-downs differ in purity"      (replaces the 4-panel QC figure)
Figure 2  "What survives once that is corrected"     (replaces the 4-panel deconvolution figure)

Everything plotted from the striatal workbook except the two items marked
UPSTREAM, which come from the snRNA-seq deconvolution run elsewhere.
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PAN, VC, GREY, INK, MUT = '#2a78d6', '#eb6834', '#9aa3ad', '#15181c', '#5c6572'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.labelsize': 11,
    'axes.titlesize': 12.5, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': '#aab0b7', 'xtick.color': MUT, 'ytick.color': MUT,
    'axes.labelcolor': INK, 'text.color': INK, 'figure.dpi': 200})

# ---------------------------------------------------------------- data
s = pd.read_excel('str.xlsx', sheet_name='|FC| >2')
pan = [f'JD130_Str_{i}.norm' for i in range(1, 5)]
vc = [f'JD1884_Str_{i}.norm' for i in range(1, 5)]
X = s.set_index('Gene')[pan + vc]
X = X[~X.index.duplicated()]
lfc = -s.set_index('Gene')['condition-JD130_StrvsJD1884_Str.log2FC']
lfc = lfc[~lfc.index.duplicated()]          # positive = higher in VC

PANELS = {
 'Astrocyte': ['Aldh1l1','Slc1a3','Slc1a2','Gja1','Aqp4','Gfap','S100b','Sox9','Acsbg1','Slc6a11',
   'Atp1b2','Mlc1','Glul','Ntsr2','Fgfr3','Slc7a10','Gjb6','Aldoc','Cpe','Htra1','Mt1','Mt2','Id3',
   'Slc39a12','Cldn10','Prdx6','Atp1a2','Plcd4','Paqr6','Ppp1r3c','Dio2','Bcan','Ttyh1','Slc4a4',
   'Agt','Grin2c','Cbs','Gpr37l1','Sparcl1','Entpd2'],
 'Endothelial': ['Cldn5','Pecam1','Flt1','Slco1a4','Ly6c1','Vwf','Cdh5','Esam','Tie1','Podxl','Cav1','Emcn'],
 'Mural': ['Pdgfrb','Kcnj8','Acta2','Myh11','Des','Vtn','Abcc9','Rgs5','Higd1b','Ndufa4l2','Art3','Atp13a5'],
 'Microglia': ['Cx3cr1','P2ry12','Csf1r','Aif1','Tmem119','Hexb','C1qa','C1qb','C1qc','Ctss','Fcrls','Olfml3','Trem2','Laptm5'],
 'Neuron': ['Rbfox3','Snap25','Syt1','Syn1','Meg3','Camk2a','Gad1','Gad2','Nefl','Eno2','Tubb3','Stmn2',
   'Drd1','Drd2','Pde10a','Nrgn','Sncb','Atp1a3','Syp','Rab3a','Chga','Scn2a','Kcnc1','Slc17a7'],
 'Oligodendrocyte': ['Plp1','Mbp','Sox10','Mog','Mag','Cldn11','Mal','Cnp','Ugt8a','Trf','Apod','Opalin','Aspa','St18'],
}
med = {k: float(lfc[lfc.index.isin(v)].dropna().median()) for k, v in PANELS.items()}
npan = {k: int(lfc[lfc.index.isin(v)].dropna().shape[0]) for k, v in PANELS.items()}

# per-sample astrocyte-signal index: median-of-ratios size factors on astrocyte-specific genes
A = X.loc[X.index.isin(PANELS['Astrocyte'])]
A = A[(A > 50).all(axis=1)]
sf = (A.div(np.exp(np.log(A).mean(axis=1)), axis=0)).median()
sf = sf / sf[pan].mean() * 100          # express as % of the mean pan sample

# =================================================================== FIGURE 1
fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.3),
                       gridspec_kw={'width_ratios': [1.05, 1.5, 0.75], 'wspace': 0.42})

# -- a: per-sample astrocyte signal
a = ax[0]
vals = list(sf[pan]) + list(sf[vc])
cols = [PAN]*4 + [VC]*4
lab = [f'pan {i}' for i in range(1, 5)] + [f'VC {i}' for i in range(1, 5)]
a.barh(range(8), vals, color=cols, height=.66)
a.set_yticks(range(8)); a.set_yticklabels(lab, fontsize=10)
a.invert_yaxis(); a.axvline(100, color=MUT, ls='--', lw=1)
for i, v in enumerate(vals):
    a.text(v + 2.5, i, f'{v:.0f}', va='center', fontsize=9.5, color=MUT)
a.set_xlim(0, 135); a.set_xlabel('astrocyte signal, % of mean pan sample')
a.set_title('a   Every VC sample carries less\n     astrocyte signal', loc='left', fontweight='bold')

# -- b: which way each cell class moves
b = ax[1]
order = ['Astrocyte', 'Endothelial', 'Mural', 'Microglia', 'Neuron', 'Oligodendrocyte']
y = np.arange(len(order))
for i, k in enumerate(order):
    col = PAN if med[k] < 0 else (GREY if abs(med[k]) < 0.5 else VC)
    b.plot([0, med[k]], [i, i], color=col, lw=2.6, solid_capstyle='round', zorder=2)
    b.scatter(med[k], i, s=90, color=col, zorder=3)
    b.text(med[k] + (0.16 if med[k] > 0 else -0.16), i, f'{med[k]:+.2f}',
           va='center', ha='left' if med[k] > 0 else 'right', fontsize=9.5, color=MUT)
b.axvline(0, color=MUT, ls='--', lw=1)
b.set_yticks(y); b.set_yticklabels([f'{k}  (n={npan[k]})' for k in order], fontsize=10)
b.invert_yaxis(); b.set_xlim(-1.6, 3.6)
b.set_xlabel('median log$_2$ fold change,  VC ÷ pan')
b.set_title('b   Astrocyte genes fall, every other cell type rises —\n     a purity gradient, not a pathway',
            loc='left', fontweight='bold')

# -- c: deconvolution summary (upstream)
c = ax[2]
c.bar([0, 1], [96, 61], color=[PAN, VC], width=.6)
c.plot([1, 1], [30, 85], color=INK, lw=1.6)
c.plot([0.88, 1.12], [30, 30], color=INK, lw=1.6)
c.plot([0.88, 1.12], [85, 85], color=INK, lw=1.6)
for x, v in [(0, 96), (1, 61)]:
    c.text(x, v + 3, f'{v}%', ha='center', fontsize=11, fontweight='bold')
c.set_xticks([0, 1]); c.set_xticklabels(['pan\nJD130', 'VC\nJD1884'], fontsize=10)
c.set_ylim(0, 108); c.set_ylabel('% astrocyte in the pull-down')
c.set_title('c   snRNA-seq\n     deconvolution', loc='left', fontweight='bold')
c.text(1.17, 30, '30', va='center', fontsize=9, color=MUT)
c.text(1.17, 85, '85', va='center', fontsize=9, color=MUT)
c.set_xlim(-0.6, 1.62)

fig.suptitle('The two pull-downs differ in purity, not in biology',
             x=0.012, y=0.995, ha='left', fontsize=16, fontweight='bold')
fig.text(0.012, 0.015,
         'Striatum, 4 vs 4 replicates. a, b computed from the supplied counts; the cortical workbook is pre-filtered to '
         '|log$_2$FC| > 1 and cannot support these panels. c from the snRNA-seq deconvolution run upstream; '
         'its whisker spans the individual VC samples, 30 to 85%.',
         fontsize=8.6, color=MUT)
fig.subplots_adjust(top=0.78, bottom=0.21, left=0.085, right=0.985)
fig.savefig('Fig1_purity_simple.png', bbox_inches='tight')
fig.savefig('Fig1_purity_simple.pdf', bbox_inches='tight')
print('Fig1 written. per-sample astrocyte signal:',
      {k: round(v) for k, v in zip(lab, vals)})

# =================================================================== FIGURE 2
STEPS = [('Nominally significant', 395, 'moderated test, FDR < 0.05'),
         ('Effect big enough', 147, 'dilution-corrected |Δ| ≥ 1'),
         ('Seen in both regions', 121, 'same direction, cortex and striatum'),
         ('Not driven by one animal', 121, 'stable across 16 leave-one-out fits'),
         ('Actually expressed', 120, '≥ 50 counts'),
         ('Astrocyte-specific', 5, 'interpretable after dilution: 3 up, 2 down')]

OFF = -0.66                      # dilution offset used upstream
HITS = [('St6galnac5', 0.76, 0.0023, 'up'), ('Gabrb1', 0.44, 0.033, 'up'),
        ('Ncan', 0.24, 0.57, 'up'), ('Plscr4', -1.90, 3.5e-13, 'down'),
        ('Gfap', -2.77, 3.3e-21, 'down')]

fig, ax = plt.subplots(1, 2, figsize=(13.6, 4.6),
                       gridspec_kw={'width_ratios': [1.25, 1], 'wspace': 0.3})
a = ax[0]
ys = np.arange(len(STEPS))
for i, (nm, n, rule) in enumerate(STEPS):
    col = VC if i == len(STEPS) - 1 else GREY
    a.barh(i, n, color=col, height=.6)
    a.text(n + 8, i, str(n), va='center', fontsize=12, fontweight='bold',
           color=VC if i == len(STEPS) - 1 else INK)
a.set_yticks(ys); a.set_yticklabels([s[0] for s in STEPS], fontsize=10.5)
a.invert_yaxis(); a.set_xlim(0, 450); a.set_xlabel('genes remaining')
a.set_title('a   395 candidates, 5 survivors', loc='left', fontweight='bold')

b = ax[1]
for i, (g, raw, fdr, d) in enumerate(HITS):
    corr = raw - OFF
    col = VC if d == 'up' else PAN
    b.plot([raw, corr], [i, i], color=GREY, lw=1.4, ls=':', zorder=1)
    b.scatter(raw, i, s=52, facecolor='white', edgecolor=GREY, zorder=2, lw=1.4)
    b.scatter(corr, i, s=115, color=col, zorder=3)
    outer = max(raw, corr) + 0.26 if corr > 0 else min(raw, corr) - 0.26
    b.text(outer, i, f'{corr:+.2f}', va='center',
           ha='left' if corr > 0 else 'right', fontsize=10, color=col, fontweight='bold')
b.axvline(0, color=MUT, ls='--', lw=1)
b.set_yticks(range(len(HITS)))
b.set_yticklabels([f'{g}' for g, *_ in HITS], fontsize=11, style='italic')
b.invert_yaxis(); b.set_xlim(-4.3, 2.6)
b.set_xlabel('log$_2$ fold change, VC ÷ pan')
b.set_title('b   The five, before and after correcting for dilution', loc='left', fontweight='bold')
b.scatter([], [], s=52, facecolor='white', edgecolor=GREY, lw=1.4, label='as measured')
b.scatter([], [], s=115, color=GREY, label='after dilution correction')
b.legend(frameon=False, fontsize=9, loc='upper left', handletextpad=.4, borderaxespad=.1)

fig.suptitle('Correcting for dilution leaves five interpretable genes',
             x=0.012, y=0.995, ha='left', fontsize=16, fontweight='bold')
fig.text(0.012, 0.015,
         'Filters, in order: moderated test FDR < 0.05 · dilution-corrected |Δ| ≥ 1 · same direction in cortex and striatum · '
         'stable across 16 leave-one-out fits · ≥ 50 counts · astrocyte-specific in the snRNA-seq reference.\n'
         'a: counts from the upstream analysis. b: fold changes from the supplied striatal table, corrected by subtracting the '
         '−0.66 astrocyte-dilution offset. Ncan does not reach significance on its own (FDR 0.57).',
         fontsize=8.6, color=MUT)
fig.subplots_adjust(top=0.82, bottom=0.2, left=0.13, right=0.97)
fig.savefig('Fig2_survivors_simple.png', bbox_inches='tight')
fig.savefig('Fig2_survivors_simple.pdf', bbox_inches='tight')
print('Fig2 written. corrected values:', {g: round(r - OFF, 2) for g, r, _, _ in HITS})
