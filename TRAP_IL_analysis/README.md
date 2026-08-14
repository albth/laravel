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
* Cut-offs: **|log2FC| > 1** (as requested) **and** a detection floor of
  **group-mean normalised counts ≥ 10** in at least one genotype. The floor removes
  genes driven by 1–3 stray reads (e.g. `Il9r`, max 2.1 counts). It is set by
  `EXPR_FLOOR` in the script — set it to `0` to reproduce the unfiltered list.

23 genes pass (9 ligands, 14 receptors): **Cortex** n = 10, all vcAst-enriched;
**Striatum** n = 22, 17 vcAst-enriched / 5 pan-astrocyte-enriched.

## Figures

All outputs are written to `figures/` as 300 dpi PNG **and** vector PDF with Type-42
(editable) fonts, sized for a single journal column.

**`Fig1_IL_logFC_heatmap`** — Heatmap of interleukin ligand and receptor log2FC in
cortex and striatum. Colour, log2FC of vcAst (JD1884) over pan-astrocyte (JD130);
warm = vcAst-enriched, cool = pan-astrocyte-enriched, scale saturated at ±6. Rows are
split into ligands and receptors and ordered by mean vcAst enrichment. Grey, gene did
not reach |log2FC| > 1 in that region. \*FDR < 0.05, \*\*FDR < 0.001,
\*\*\*FDR < 10⁻⁴; absence of a star means FDR ≥ 0.05 or not estimable.

**`Fig2_IL_dotplot`** — Dot plot of the same genes and contrasts. Dot colour, log2FC
(vcAst / pan-astrocyte); dot area, −log₁₀FDR. Dots with FDR ≥ 0.05 carry a grey ring;
× marks an effect that passes the fold-change cut but for which FDR was not estimable
(zero counts in one genotype).

**`Fig3_IL_expression_zscore`** — Per-replicate expression of the same genes,
n = 4 per genotype per region, plotted as the row z-score of log₂(normalised
counts + 1) computed within each region. Grey rows in the cortical panel are genes
absent from the pre-filtered cortical table.

**`source_data_IL_genes.csv`** — per-gene, per-region source data behind all three
figures: sign-flipped log2FC, FDR, group-mean expression, and the pass/fail flag.

## Notes for interpretation

* In cortex, **every** interleukin gene passing the cut is vcAst-enriched — the
  cortical IL signature is unidirectional.
* `Il1rapl1`, `Il1rapl2`, `Il10ra`, `Il34`, `Il4` and `Il16` are vcAst-enriched in
  **both** regions with FDR < 0.05, and are the most defensible cross-regional calls.
* `Il20rb` and `Il1rl2` show the largest striatal effects but sit at 0 counts in all
  four JD130 samples, which inflates log2FC and makes the FDR unstable — treat as
  qualitative presence/absence, and confirm by qPCR or RNAscope before quoting a
  fold change.
