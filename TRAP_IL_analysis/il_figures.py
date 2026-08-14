"""
TRAP-seq: interleukin ligand / receptor signature of vascular-contacting astrocytes
==================================================================================

Comparison in the source workbooks is  JD130 (pan-astrocyte)  vs  JD1884 (vcAst),
so a NEGATIVE log2FC means vcAst-enriched.  Every figure below reports

        log2FC(plotted) = -1 x log2FC(source)

so that POSITIVE = enriched in JD1884 (vcAst) and NEGATIVE = enriched in
JD130 (pan-astrocyte).

Genes shown: interleukin ligands and interleukin receptors (symbol ^Il<digit>),
passing |log2FC| > 1 with a detectable expression floor (see EXPR_FLOOR).

Outputs (300 dpi PNG + vector PDF, Type-42 fonts) into ./figures.
"""

import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

HERE = Path(__file__).resolve().parent
FIGDIR = HERE / "figures"
FIGDIR.mkdir(exist_ok=True)

CTX_XLSX = HERE / "JD130_CtxvsJD1884_Ctx.xlsx"
STR_XLSX = HERE / "JD130_StrvsJD1884_Str.xlsx"

LFC_CUT = 1.0     # |log2FC| threshold requested
EXPR_FLOOR = 10.0  # min. group-mean normalised counts in at least one genotype
CLIP = 6.0        # colour-scale saturation for log2FC

# --------------------------------------------------------------------------- #
# publication rcParams
# --------------------------------------------------------------------------- #
mpl.rcParams.update({
    "pdf.fonttype": 42,          # editable text in Illustrator
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

INK = "#1b1b1b"
INK_MUTED = "#6b6b6b"
SURFACE = "#ffffff"
MISSING = "#e9e9ec"   # gene did not pass the threshold in that region

# Diverging ramp: cool = pan-astrocyte, neutral grey midpoint, warm = vcAst.
CMAP = LinearSegmentedColormap.from_list(
    "vcast_div",
    ["#1a4f8f", "#3f7cc0", "#9dc0e0", "#eeeeee", "#f2b48a", "#d9612c", "#8c2d04"],
)
NORM = TwoSlopeNorm(vmin=-CLIP, vcenter=0.0, vmax=CLIP)

# --------------------------------------------------------------------------- #
# gene curation
# --------------------------------------------------------------------------- #
IL_PATTERN = re.compile(r"^Il\d", re.IGNORECASE)

# ^Il<digit>… symbols that are NOT interleukin ligands or receptors
NOT_IL = {
    "Il4i1",    # interleukin-4 induced 1 — an L-amino-acid oxidase
    "Il4i1b",
    "Il1bos",   # Il1b opposite strand lncRNA
    "Il18bp",   # soluble binding protein, not a signalling receptor
}
# receptor subunits whose symbol does not match ^Il\d+r
EXTRA_RECEPTORS = {"Il6st"}          # gp130
# ^Il\d+r symbols that are secreted ligands, not receptors
LIGAND_EXCEPTIONS = {"Il1rn"}        # IL-1 receptor antagonist (secreted)


def classify(symbol: str) -> str:
    """Ligand vs receptor for an interleukin gene symbol."""
    if symbol in EXTRA_RECEPTORS:
        return "Receptor"
    if symbol in LIGAND_EXCEPTIONS:
        return "Ligand"
    return "Receptor" if re.match(r"^Il\d+r", symbol, re.IGNORECASE) else "Ligand"


# --------------------------------------------------------------------------- #
# data loading
# --------------------------------------------------------------------------- #
def load_region(xlsx: Path, sheet: str, region: str) -> pd.DataFrame:
    """Return tidy per-gene table with the log2FC sign flipped to vcAst-positive."""
    df = pd.read_excel(xlsx, sheet_name=sheet)
    lfc_col = next(c for c in df.columns if c.endswith(".log2FC"))
    fdr_col = next(c for c in df.columns if c.endswith(".FDR"))
    pan_cols = [c for c in df.columns if c.startswith("JD130_")]
    vc_cols = [c for c in df.columns if c.startswith("JD1884_")]

    il = df[df["Gene"].astype(str).str.match(IL_PATTERN, na=False)].copy()
    il = il[~il["Gene"].isin(NOT_IL)]

    out = pd.DataFrame({
        "Gene": il["Gene"].astype(str),
        # sign flip: positive = vcAst (JD1884) enriched
        "log2FC": -il[lfc_col].astype(float),
        "FDR": il[fdr_col].astype(float),
        "mean_pan": il[pan_cols].mean(axis=1).astype(float),
        "mean_vc": il[vc_cols].mean(axis=1).astype(float),
    })
    out["Region"] = region
    out["Class"] = out["Gene"].map(classify)
    out["max_expr"] = out[["mean_pan", "mean_vc"]].max(axis=1)
    # collapse duplicated symbols (multiple Ensembl IDs) on the largest effect
    out = (out.sort_values("log2FC", key=lambda s: s.abs(), ascending=False)
              .drop_duplicates("Gene")
              .reset_index(drop=True))
    return out


# The Str workbook stores the complete (unfiltered) striatal table on the sheet
# labelled "|FC| >2"; its "|FC| >1" sheet is a duplicate of the cortical data.
ctx = load_region(CTX_XLSX, "|FC| >1", "Cortex")
stri = load_region(STR_XLSX, "|FC| >2", "Striatum")

tidy = pd.concat([ctx, stri], ignore_index=True)
tidy["passes"] = (tidy["log2FC"].abs() > LFC_CUT) & (tidy["max_expr"] >= EXPR_FLOOR)

kept = tidy[tidy["passes"]].copy()

# gene order: ligands then receptors, each sorted by mean vcAst enrichment
mean_lfc = kept.groupby("Gene")["log2FC"].mean()
cls = kept.drop_duplicates("Gene").set_index("Gene")["Class"]
order = []
for grp in ("Ligand", "Receptor"):
    genes = [g for g in mean_lfc.index if cls[g] == grp]
    order += sorted(genes, key=lambda g: -mean_lfc[g])
GENES = order
CLASS_OF = {g: cls[g] for g in GENES}
N_LIG = sum(1 for g in GENES if CLASS_OF[g] == "Ligand")

REGIONS = ["Cortex", "Striatum"]
lookup = {(r["Gene"], r["Region"]): r for _, r in tidy.iterrows()}


def cell(gene, region):
    return lookup.get((gene, region))


def stars(fdr):
    if fdr is None or not np.isfinite(fdr):
        return ""
    if fdr < 1e-4:
        return "***"
    if fdr < 1e-3:
        return "**"
    if fdr < 0.05:
        return "*"
    return ""


# --------------------------------------------------------------------------- #
# shared drawing helpers
# --------------------------------------------------------------------------- #
def style_gene_axis(ax, n):
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_yticklabels(GENES, fontstyle="italic", color=INK)
    ax.set_ylim(n, 0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0, pad=2)


def class_bracket(fig, ax, n):
    """Ligand / receptor bracket placed just outside the gene labels.

    Positioned in figure coordinates after a draw pass so it clears the
    y tick labels whatever their width.
    """
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    left = min(t.get_window_extent().transformed(inv).x0
               for t in ax.get_yticklabels())
    x = left - 0.014

    def yfig(row):
        return inv.transform(ax.transData.transform((0, row)))[1]

    for label, lo, hi in (("Ligands", 0, N_LIG), ("Receptors", N_LIG, n)):
        if hi <= lo:
            continue
        y0, y1 = yfig(lo + 0.15), yfig(hi - 0.15)
        fig.add_artist(Line2D([x, x], [y0, y1], color=INK_MUTED, lw=1.0,
                              solid_capstyle="butt"))
        fig.text(x - 0.012, (y0 + y1) / 2, label, rotation=90, va="center",
                 ha="center", color=INK_MUTED, fontsize=7)


def colorbar(fig, rect, label="log$_2$FC   (vcAst / pan-astrocyte)"):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP),
                      cax=cax, orientation="vertical",
                      ticks=[-CLIP, -3, 0, 3, CLIP])
    cb.set_label(label, fontsize=6.5, color=INK, labelpad=4)
    cb.ax.set_yticklabels([f"≤ -{CLIP:.0f}", "-3", "0", "3", f"≥ {CLIP:.0f}"])
    cb.ax.tick_params(labelsize=6, length=2, width=0.5, color=INK_MUTED)
    cb.outline.set_linewidth(0.5)
    cb.outline.set_edgecolor(INK_MUTED)
    cb.ax.text(0.5, 1.06, "vcAst ↑", transform=cb.ax.transAxes,
               ha="center", va="bottom", fontsize=6, color="#8c2d04")
    cb.ax.text(0.5, -0.06, "pan-astro ↑", transform=cb.ax.transAxes,
               ha="center", va="top", fontsize=6, color="#1a4f8f")
    return cb


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(FIGDIR / f"{stem}.{ext}")
    plt.close(fig)
    print(f"  wrote figures/{stem}.png / .pdf")


# --------------------------------------------------------------------------- #
# Figure 1 — heatmap of log2FC
# --------------------------------------------------------------------------- #
def figure_heatmap():
    n = len(GENES)
    fig = plt.figure(figsize=(3.15, 0.155 * n + 1.15))
    ax = fig.add_axes([0.30, 0.075, 0.30, 0.855])

    for j, region in enumerate(REGIONS):
        for i, gene in enumerate(GENES):
            rec = cell(gene, region)
            passes = rec is not None and bool(rec["passes"])
            face = CMAP(NORM(rec["log2FC"])) if passes else MISSING
            ax.add_patch(Rectangle((j + 0.03, i + 0.03), 0.94, 0.94,
                                   facecolor=face, edgecolor=SURFACE,
                                   linewidth=0.8, zorder=2))
            if passes:
                s = stars(rec["FDR"])
                if s:
                    lum = float(np.mean(mpl.colors.to_rgb(face)[:3]))
                    ax.text(j + 0.5, i + 0.5, s, ha="center", va="center",
                            fontsize=6, color=INK if lum > 0.55 else "white",
                            zorder=3)

    ax.set_xlim(0, len(REGIONS))
    ax.set_xticks(np.arange(len(REGIONS)) + 0.5)
    ax.set_xticklabels(REGIONS, color=INK)
    style_gene_axis(ax, n)
    ax.xaxis.set_ticks_position("top")
    ax.tick_params(axis="x", pad=3)
    class_bracket(fig, ax, n)

    # ligand / receptor divider
    ax.plot([0, len(REGIONS)], [N_LIG, N_LIG], color=SURFACE, lw=2.2, zorder=4)

    ax.set_title("Interleukin ligands and receptors\nvcAst (JD1884) vs pan-astrocyte (JD130)",
                 pad=22, fontsize=8, color=INK, loc="left", x=-0.02)

    colorbar(fig, [0.70, 0.50, 0.035, 0.30])

    legend = fig.legend(
        handles=[Patch(facecolor=MISSING, edgecolor="none",
                       label=f"|log$_2$FC| ≤ {LFC_CUT:.0f}")],
        loc="lower left", bbox_to_anchor=(0.655, 0.395), frameon=False,
        handlelength=1.0, handletextpad=0.5,
    )
    for t in legend.get_texts():
        t.set_color(INK)
    fig.text(0.655, 0.335, "FDR", fontsize=6.5, color=INK)
    for k, txt in enumerate(("*    < 0.05", "**   < 0.001", "***  < 10$^{-4}$")):
        fig.text(0.655, 0.288 - 0.040 * k, txt, fontsize=6.5, color=INK,
                 family="monospace")
    save(fig, "Fig1_IL_logFC_heatmap")


# --------------------------------------------------------------------------- #
# Figure 2 — dot plot (size = significance, colour = log2FC)
# --------------------------------------------------------------------------- #
SIZE_BREAKS = [1.3, 3, 6, 12]          # -log10(FDR)


def size_from_fdr(fdr):
    if fdr is None or not np.isfinite(fdr) or fdr <= 0:
        return None
    return float(np.clip(-np.log10(fdr), 0, 14))


def area(v):
    return 6 + 13.5 * v


def figure_dotplot():
    n = len(GENES)
    fig = plt.figure(figsize=(3.35, 0.155 * n + 1.15))
    ax = fig.add_axes([0.30, 0.075, 0.26, 0.855])

    # significant and non-significant dots are drawn separately so that
    # FDR ≥ 0.05 can carry a grey ring (identity never by colour alone)
    sig = {"x": [], "y": [], "c": [], "s": []}
    ns = {"x": [], "y": [], "c": [], "s": []}
    for j, region in enumerate(REGIONS):
        for i, gene in enumerate(GENES):
            rec = cell(gene, region)
            if rec is None or not rec["passes"]:
                continue
            v = size_from_fdr(rec["FDR"])
            if v is None:      # effect passes but FDR not estimable
                ax.plot(j + 0.5, i + 0.5, marker="x", ms=2.6, mew=0.7,
                        color=INK_MUTED, zorder=3)
                continue
            tgt = sig if rec["FDR"] < 0.05 else ns
            tgt["x"].append(j + 0.5)
            tgt["y"].append(i + 0.5)
            tgt["c"].append(rec["log2FC"])
            tgt["s"].append(area(v))

    ax.set_axisbelow(True)
    for i in range(n):
        ax.axhline(i + 0.5, color="#f0f0f2", lw=0.5, zorder=0)
    for j in range(len(REGIONS)):
        ax.axvline(j + 0.5, color="#f0f0f2", lw=0.5, zorder=0)

    ax.scatter(ns["x"], ns["y"], s=ns["s"], c=ns["c"], cmap=CMAP, norm=NORM,
               linewidths=0.6, edgecolors=INK_MUTED, zorder=3)
    ax.scatter(sig["x"], sig["y"], s=sig["s"], c=sig["c"], cmap=CMAP, norm=NORM,
               linewidths=0.5, edgecolors=SURFACE, zorder=4)

    ax.set_xlim(0, len(REGIONS))
    ax.set_xticks(np.arange(len(REGIONS)) + 0.5)
    ax.set_xticklabels(REGIONS, color=INK)
    style_gene_axis(ax, n)
    ax.xaxis.set_ticks_position("top")
    ax.tick_params(axis="x", pad=3)
    class_bracket(fig, ax, n)
    ax.axhline(N_LIG, color=INK_MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)

    ax.set_title("Interleukin ligands and receptors\nvcAst (JD1884) vs pan-astrocyte (JD130)",
                 pad=22, fontsize=8, color=INK, loc="left", x=-0.02)

    colorbar(fig, [0.66, 0.53, 0.033, 0.28])

    handles = [Line2D([], [], linestyle="none", marker="o", color=INK_MUTED,
                      markerfacecolor="#c9c9cf", markeredgecolor=SURFACE,
                      markeredgewidth=0.5, markersize=np.sqrt(area(b)),
                      label=f"{b:g}" + ("  (FDR 0.05)" if b == 1.3 else ""))
               for b in SIZE_BREAKS]
    handles += [
        Line2D([], [], linestyle="none", marker="o", color=INK_MUTED,
               markerfacecolor="#c9c9cf", markeredgecolor=INK_MUTED,
               markeredgewidth=0.6, markersize=np.sqrt(area(0.8)),
               label="n.s. (grey ring)"),
        Line2D([], [], linestyle="none", marker="x", color=INK_MUTED,
               markersize=2.6, markeredgewidth=0.7, label="FDR not estimable"),
    ]
    leg = fig.legend(handles=handles, title="-log$_{10}$FDR",
                     loc="upper left", bbox_to_anchor=(0.63, 0.47),
                     frameon=False, handlelength=1.0, handletextpad=0.6,
                     labelspacing=0.75, borderpad=0.2)
    leg.get_title().set_fontsize(6.5)
    leg.get_title().set_color(INK)
    for t in leg.get_texts():
        t.set_color(INK)
    save(fig, "Fig2_IL_dotplot")


# --------------------------------------------------------------------------- #
# Figure 3 — per-sample expression (row z-score), one panel per region
# --------------------------------------------------------------------------- #
EXPR_CMAP = LinearSegmentedColormap.from_list(
    "expr", ["#2c3d8f", "#7d92c9", "#f1f1f4", "#e9a05c", "#9a3a06"])


def load_samples(xlsx, sheet):
    df = pd.read_excel(xlsx, sheet_name=sheet)
    cols = [c for c in df.columns if c.startswith("JD130_") or c.startswith("JD1884_")]
    df = df[df["Gene"].astype(str).str.match(IL_PATTERN, na=False)]
    df = df[df["Gene"].isin(GENES)]
    mat = df.groupby(df["Gene"].astype(str))[cols].mean()
    return mat.reindex([g for g in GENES if g in mat.index]), cols


def figure_expression():
    n = len(GENES)
    fig = plt.figure(figsize=(4.9, 0.155 * n + 1.3))
    panels = [
        ("Cortex", load_samples(CTX_XLSX, "|FC| >1"), [0.235, 0.075, 0.235, 0.83]),
        ("Striatum", load_samples(STR_XLSX, "|FC| >2"), [0.505, 0.075, 0.235, 0.83]),
    ]
    vlim = 1.8
    for k, (region, (mat, cols), rect) in enumerate(panels):
        ax = fig.add_axes(rect)
        z = np.full((n, len(cols)), np.nan)
        for i, gene in enumerate(GENES):
            if gene not in mat.index:
                continue
            v = np.log2(mat.loc[gene].to_numpy(dtype=float) + 1)
            sd = v.std(ddof=0)
            if sd > 0:
                z[i] = (v - v.mean()) / sd
        ax.imshow(np.ma.masked_invalid(z), aspect="auto", cmap=EXPR_CMAP,
                  vmin=-vlim, vmax=vlim, interpolation="nearest",
                  extent=(0, len(cols), n, 0))
        ax.set_facecolor(MISSING)

        ax.set_xticks([])
        ax.axvline(len(cols) / 2, color=SURFACE, lw=1.6)
        ax.axhline(N_LIG, color=SURFACE, lw=1.6)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_ylim(n, 0)
        ax.set_yticks([])
        if k == 0:
            style_gene_axis(ax, n)
            class_bracket(fig, ax, n)
        # genotype bar above each panel
        for x0, x1, col, lab in [
                (0, len(cols) / 2, "#5b6b7f", "JD130\npan-astro"),
                (len(cols) / 2, len(cols), "#c26a33", "JD1884\nvcAst")]:
            ax.add_patch(Rectangle((x0 + 0.06, -0.022 * n), x1 - x0 - 0.12, 0.013 * n,
                                   facecolor=col, edgecolor="none",
                                   clip_on=False, transform=ax.transData))
            ax.text((x0 + x1) / 2, -0.030 * n, lab, ha="center", va="bottom",
                    fontsize=6, color=col, clip_on=False, linespacing=1.3)
        ax.text(len(cols) / 2, -0.085 * n, region, ha="center", va="bottom",
                fontsize=8, color=INK, clip_on=False)

    cax = fig.add_axes([0.79, 0.45, 0.022, 0.26])
    cb = fig.colorbar(mpl.cm.ScalarMappable(
        norm=mpl.colors.Normalize(-vlim, vlim), cmap=EXPR_CMAP), cax=cax,
        ticks=[-vlim, 0, vlim])
    cb.set_label("row z-score\nlog$_2$(norm. counts + 1)", fontsize=6.5, color=INK)
    cb.ax.set_yticklabels([f"≤ -{vlim:g}", "0", f"≥ {vlim:g}"])
    cb.ax.tick_params(labelsize=6, length=2, width=0.5, color=INK_MUTED)
    cb.outline.set_linewidth(0.5)
    cb.outline.set_edgecolor(INK_MUTED)

    leg = fig.legend(
        handles=[Patch(facecolor=MISSING, edgecolor="none",
                       label="not in the |log$_2$FC| > 1\ncortical table")],
        loc="upper left", bbox_to_anchor=(0.775, 0.40), frameon=False,
        handlelength=1.0, handletextpad=0.5)
    for t in leg.get_texts():
        t.set_color(INK)
        t.set_fontsize(6.5)

    fig.text(0.235, 1.055, "Interleukin ligand / receptor expression per replicate",
             fontsize=8, color=INK, ha="left", va="bottom")
    save(fig, "Fig3_IL_expression_zscore")


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    print(f"IL genes passing |log2FC| > {LFC_CUT:g} and expression ≥ {EXPR_FLOOR:g}: "
          f"{len(GENES)}  ({N_LIG} ligands, {len(GENES) - N_LIG} receptors)")
    for region in REGIONS:
        sub = kept[kept["Region"] == region]
        print(f"  {region:9s} n = {len(sub):2d}  "
              f"vcAst-enriched {int((sub['log2FC'] > 0).sum()):2d} / "
              f"pan-astro-enriched {int((sub['log2FC'] < 0).sum()):2d}")

    figure_heatmap()
    figure_dotplot()
    figure_expression()

    src = tidy[tidy["Gene"].isin(GENES)].copy()
    src = src.assign(order=src["Gene"].map({g: i for i, g in enumerate(GENES)}))
    src = src.sort_values(["order", "Region"]).drop(columns="order")
    src = src.rename(columns={
        "log2FC": "log2FC_vcAst_over_panAstro",
        "mean_pan": "mean_JD130_panAstro",
        "mean_vc": "mean_JD1884_vcAst",
        "passes": "passes_filters",
    })
    src.to_csv(FIGDIR / "source_data_IL_genes.csv", index=False)
    print("  wrote figures/source_data_IL_genes.csv")
