# vcAst TRAP-seq: composition QC

QC of the JD130 (pan-astrocyte) vs JD1884 (vcAst) TRAP-seq comparison, source
workbooks on branch `claude/trap-rnaseq-il-heatmap-dotplot-eg89ik`
(`TRAP_IL_analysis/JD130_*vsJD1884_*.xlsx`). Sign convention follows that
analysis: positive log2FC = enriched in JD1884 (vcAst).

## Finding

The comparison is dominated by immunoprecipitation purity, not by biology.

| Marker panel | n | median log2FC (vcAst ÷ pan), striatum |
|---|---|---|
| Astrocyte | 40 | **−0.92** (bootstrap 95% CI −1.05 to −0.74) |
| Endothelial | 10 | +0.11 |
| Mural / pericyte | 7 | +0.96 |
| Microglia | 13 | +1.90 |
| Neuron | 24 | +2.12 |
| Oligodendrocyte | 13 | +2.85 |

39 of 40 astrocyte-specific genes are lower in the vcAst IP; every other cell
type is higher. A two-component mixing model
(`pan = b·Ast + (1−b)·Bg`, `vc = a·Ast + (1−a)·Bg`) puts the JD130 IP at
85–93% astrocyte and the JD1884 IP at 45–49%, consistently across the neuron,
oligodendrocyte and microglia panels.

Endothelial markers are the one panel that does **not** rise. If the vcAst line
were pulling vessel-associated material into the IP, endothelium is what should
have come with it, so physical vascular proximity does not explain the pattern.

## Consequence for interpretation

In a mixture, cell-type-specific genes are diluted and broadly-expressed genes
are buffered, so a gene's apparent fold change mostly reports how
astrocyte-restricted it is. Any pathway built from the raw table will recover
"broadly expressed" (mitochondrial, ER/Golgi, lipid) rather than vcAst biology.
A flat offset correction does **not** fix this; it is valid only for genes with
negligible expression in the contaminating cell types.

## What is still interpretable

Restricting to astrocyte-specific genes and recomputing DESeq2-style
median-of-ratios size factors on those genes only (`astro_normalised_DE.csv`,
58 genes, 4 vs 4 replicates, striatum):

* **Nothing reaches FDR < 0.05.**
* Strongest nominal signals: `Utrn` +2.55 (p = 0.001), `Dag1` +1.03 (p = 0.009),
  `Kcnj10` +0.76 (p = 0.041), `Slc1a2` +0.59 (p = 0.048), `Gfap` −1.91
  (p = 0.004), `Cbs` −0.79 (p = 0.002).
* The direction is coherent — two components of the dystrophin-associated
  complex that anchors AQP4 move up together while GFAP moves down — but this
  is hypothesis-generating only.
* If either line uses a `Gfap`-based driver, the GFAP result is circular and
  should be dropped. Worth checking before use.

## To make the dataset answer the question

1. An **input / total-RNA sample per animal**, so composition is estimable
   rather than assumed.
2. A **purity-matched comparison** — a 2-fold difference in IP specificity
   cannot be normalised away post hoc.
3. **More replicates**: the effects of interest are smaller than the artefact
   currently on top of them.

## Files

* `load.py` — loads both workbooks, flips the sign so positive = vcAst.
* `astro_specific_adj.csv` — astrocyte-specific genes with the flat offset applied.
* `astro_normalised_DE.csv` — astrocyte-content-normalised DE with BH FDR.

## Simplified presentation figures

`simple_figs.py` rebuilds the two QC figures as one-argument-each slides from
`ctx.xlsx` / `str.xlsx` (striatal counts; the cortical workbook is pre-filtered
to |log2FC| > 1 and cannot support the marker panels).

* `figures/Fig1_purity_simple.*` — **"the two pull-downs differ in purity, not
  in biology"**: per-sample astrocyte signal (pan 78/87/128/107 vs VC
  54/61/53/43, % of mean pan sample), median log2FC per cell class
  (Astrocyte −0.92, Endothelial +0.11, Mural +1.12, Microglia +1.99,
  Neuron +2.12, Oligodendrocyte +2.91), and the deconvolution estimate
  (96% vs 61%). Endothelium being flat is the diagnostic: this is general
  background, not vascular carry-over.
* `figures/Fig2_survivors_simple.*` — **"correcting for dilution leaves five
  interpretable genes"**: the attrition funnel 395 → 147 → 121 → 121 → 120 → 5
  and the five survivors before and after subtracting the −0.66 offset.

Verified against the raw striatal table:

| gene | log2FC | FDR | after −0.66 |
|---|---|---|---|
| St6galnac5 | +0.76 | 0.0023 | +1.42 |
| Gabrb1 | +0.44 | 0.033 | +1.10 |
| Ncan | +0.24 | 0.57 | +0.90 |
| Plscr4 | −1.90 | 3.5e-13 | −1.24 |
| Gfap | −2.77 | 3.3e-21 | −2.11 |

*Ncan* does not reach significance on its own; it clears the bar only after
correction. The offset itself is panel-dependent: −0.66 (their panel, 63%
purity), −0.74 (core 8 markers, 60%), −0.92 (the 40-gene panel here, 53%).

Presentation spec, scripts and the deconvolution explainer:
https://claude.ai/artifact/MTNuRiZ8x4VGte14o61FQu
