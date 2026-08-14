# Interleukin ligand / receptor signature of vcAst — TRAP-seq

Figures for the JD130 (pan-astrocyte) vs JD1884 (vascular-contacting astrocyte, vcAst)
TRAP-seq comparison, restricted to interleukin **ligands and receptors**, in **cortex**
and **striatum**.

Run with `python3 il_figures.py` (needs `pandas`, `openpyxl`, `matplotlib`, `numpy`).

## Sign convention — important

The source workbooks report `JD130 vs JD1884`, so a **negative** log2FC there means
vcAst-enriched. Every figure here plots

```
log2FC(plotted) = -1 x log2FC(source)
```

so **positive / warm = enriched in JD1884 (vcAst)** and **negative / cool = enriched in
JD130 (pan-astrocyte)**. This flip is applied once, in `load_region()`.

## Which sheet is which

The two workbooks are not laid out the same way:

| Workbook | Sheet used | Why |
|---|---|---|
| `JD130_CtxvsJD1884_Ctx.xlsx` | `\|FC\| >1` | 3,307 genes, all with \|log2FC\| > 1 — the cortical DE table |
| `JD130_StrvsJD1884_Str.xlsx` | `\|FC\| >2` | 55,536 genes — the **complete, unfiltered** striatal table |

The striatal workbook's `|FC| >1` sheet is a byte-for-byte duplicate of the cortical
data (same `JD130_Ctx…` column headers), so it is not used. The sheet labelled
`|FC| >2` in that file is in fact the full striatal result table, and the
|log2FC| > 1 cut is applied here in code.

**Consequence:** the cortical workbook is pre-filtered, so for cortex we only have
values for genes that already pass |log2FC| > 1. A grey cell in the Cortex column
therefore means "not present in the cortical |log2FC| > 1 table"; for striatum, grey
means the gene was measured and did not reach the cut.

## Gene selection

* Symbol matches `^Il<digit>` — this captures interleukin ligands (`Il1a`, `Il34`, …)
  and interleukin receptors / accessory subunits (`Il1rap`, `Il10ra`, `Il17re`, …).
* Excluded as non-ligand/non-receptor: `Il4i1` (L-amino-acid oxidase), `Il1bos`
  (lncRNA), `Il18bp` (soluble binding protein).
* `Il6st` (gp130) is classed as a receptor; `Il1rn` (IL-1RA) is classed as a
  secreted ligand despite the `r` in its symbol.
* Detection floor: **group-mean normalised counts ≥ 10** in at least one genotype,
  applied in both selections below. It removes genes driven by 1–3 stray reads
  (e.g. `Il9r`, max 2.1 counts). Set `EXPR_FLOOR = 0` in the script to disable it.

## Two selections, same figure design

`python3 il_figures.py` emits **both** sets in one run.

| Selection | Criterion (in ≥ 1 region) | Genes | File suffix |
|---|---|---|---|
| Effect size | \|log2FC\| > 1 | 23 (9 ligands, 14 receptors) | *(none)* |
| Significance | FDR < 0.05, **no fold-change cut** | 21 (11 ligands, 10 receptors) | `_significant` |

**Effect-size set** — Cortex n = 10, all vcAst-enriched; Striatum n = 22,
17 vcAst / 5 pan-astrocyte.

**Significance set** — Cortex n = 10, all vcAst-enriched; Striatum n = 17,
11 vcAst / 6 pan-astrocyte. Relative to the effect-size set it *adds* five striatal
genes that are significant with sub-2-fold effects — `Il33` (log2FC −0.80,
FDR 3×10⁻⁵), `Il18` (−0.69), `Il17rc` (−0.68), `Il13ra1` (+0.95), `Il2` (+0.99) —
and *drops* seven that clear the fold-change bar on unstable statistics
(`Il23a`, `Il1rl2`, `Il1rap`, `Il12rb1`, `Il6ra`, `Il1r2`, `Il17re`). `Il33` and
`Il18` are the most notable additions: both are highly expressed astrocyte cytokines
(mean 18,700 and 10,700 normalised counts) and both are pan-astrocyte-enriched.

Because the cortical workbook is pre-filtered at |log2FC| > 1, the cortical column of
the significance set is unchanged — a genuinely significant cortical IL gene with a
sub-2-fold effect could not appear in either figure. Recovering those needs the
unfiltered cortical table.

## Figures

All outputs are written to `figures/` as 300 dpi PNG **and** vector PDF with Type-42
(editable) fonts, sized for a single journal column.

**`Fig1_IL_logFC_heatmap[_significant]`** — Heatmap of interleukin ligand and receptor
log2FC in cortex and striatum. Colour, log2FC of vcAst (JD1884) over pan-astrocyte
(JD130); warm = vcAst-enriched, cool = pan-astrocyte-enriched, scale saturated at ±6
in both figures so the two are directly comparable. Rows are split into ligands and
receptors and ordered by mean vcAst enrichment. Plain grey, gene was measured in that
region but did not meet the selection criterion; hatched grey, gene is absent from the
pre-filtered cortical table and was therefore not tested there. \*FDR < 0.05,
\*\*FDR < 0.001, \*\*\*FDR < 10⁻⁴.

**`Fig2_IL_dotplot[_significant]`** — Dot plot of the same genes and contrasts. Dot
colour, log2FC (vcAst / pan-astrocyte); dot area, −log₁₀FDR. In the effect-size
figure, dots with FDR ≥ 0.05 carry a grey ring and × marks an effect for which FDR
was not estimable (zero counts in one genotype); neither case can arise in the
significance figure, so those keys are omitted there.

**`Fig3_IL_expression_zscore[_significant]`** — Per-replicate expression of the same
genes, n = 4 per genotype per region, plotted as the row z-score of log₂(normalised
counts + 1) computed within each region. Grey rows in the cortical panel are genes
absent from the pre-filtered cortical table.

**`source_data_IL_genes[_significant].csv`** — per-gene, per-region source data behind
each set: sign-flipped log2FC, FDR, group-mean expression, and the pass/fail flag.

## Notes for interpretation

* In cortex, **every** interleukin gene passing the cut is vcAst-enriched — the
  cortical IL signature is unidirectional.
* `Il1rapl1`, `Il1rapl2`, `Il10ra`, `Il34`, `Il4` and `Il16` are vcAst-enriched in
  **both** regions with FDR < 0.05, and are the most defensible cross-regional calls.
* `Il20rb` and `Il1rl2` show the largest striatal effects but sit at 0 counts in all
  four JD130 samples, which inflates log2FC and makes the FDR unstable — treat as
  qualitative presence/absence, and confirm by qPCR or RNAscope before quoting a
  fold change. (`Il1rl2` has no estimable FDR and so drops out of the significance
  set; `Il20rb` survives at FDR 0.042 but rests on the same zero-count denominator.)
* The significance set is the more conservative figure for a main panel; the
  effect-size set is the better supplementary companion, since between them they show
  which calls are driven by magnitude and which by consistency.
