# Vascular-coupled (juxtavascular) astrocytes in MICrONS minnie65

Analysis of astrocytes whose soma contacts a blood vessel, compared with the rest of the astrocytes and with other cell types. Built from the public data because the MICrONS Explorer and CAVE API were not reachable from the analysis environment.

## Data sources
- EM, `seg_m1300` segmentation and meshes, and nucleus segmentation: `gs://iarpa_microns/minnie/minnie65/`
- Nucleus table and neuron/non-neuron SVM: BossDB `iarpa_microns/minnie/minnie65/`
- Cell-type predictions and hand-verified identities: embedded in S. Grubb's `microns-neuroglancer-tool` (CC BY 4.0)
- Arteriole, venule and capillary-order tracings: S. Grubb's `micronsvasculature` (CC BY 4.0)

## Pipeline (`scripts/`)
1. `decode_nuc.py`, `parse_vessels.py`, `vessels.py`, `dist.py`: decode nuclei and cell types; compute distances to the traced vessels.
2. `lumen_scan.py`: detect vessel lumens across the whole volume at 512x512x640 nm (per-slice z-score, hysteresis threshold 4/6.5, 2D hole fill).
3. `nucclean.py`: remove lumen detections that overlap the nucleus segmentation.
4. `assemble.py`, `globalmask.py`, `lumenlib.py`: build the lumen surface KD-tree and a global bitmap.
5. `dist_lumen.py`, `localnull.py`: sphere-approximation distances and a first chance model.
6. `nucgap.py`, `gapstats.py`: nucleus-surface-to-lumen gap from the real nucleus shape; shape-preserving chance model (20 random rigid shifts of 10-40 um per nucleus, rejecting any overlap with lumen).
7. `astromesh.py`, `meshbatch.py`, `samevessel.py`, `mergecheck.py`: m1300 mesh analysis of 120 vc-astrocytes and 120 distant astrocytes (endfoot sites, vessel type, same-vessel geodesic test, merge check).
8. `overlay.py`, `emview.py`: EM rendering used for QC and the mitochondria survey.

Scripts expect to run from one working directory holding the intermediate files. The intermediates (lumen blocks, about 1 GB) are not committed.

## Results (`results/`)
- `vc_astrocytes_minnie65.csv`: the 413 MICrONS-predicted astrocytes whose nucleus lies within 1 um of a lumen. Coordinates are Neuroglancer voxels (4, 4, 40 nm); root IDs are m1300.
- `gap_enrichment.csv`: observed vs chance contact fraction per cell group, with bootstrap 95% CIs.
- `gapcurves.json`: observed and chance gap histograms (0.5 um bins).
- `astro_mesh_summary.csv`, `endfoot_sites.csv`, `meshbatch.json`, `mergecheck.json`: mesh analysis outputs.

## Headline numbers
- 5.6% of astrocytes have a nucleus within 1 um of a lumen vs 2.2% by chance: 2.5x enrichment (95% CI 2.2-2.8). The excess is confined to gaps under 1.5 um.
- Other glia: OPC 2.9x, oligodendrocyte 1.9x, microglia 6.8x.
- vc-astrocytes have the same total membrane area as other astrocytes, a somatic vessel sheath, and slightly fewer distal endfeet (median 1.5 vs 2). 83% have at least one distal endfoot, and 84% of those are on a different vessel segment than the soma.
- No prismatic cristae were seen in about 60 endfoot and perivascular fields at 8x8x40 nm; this is not conclusive at this resolution.

## Refined endfoot analysis and movie (November 2026)

`scripts/refined/` reruns the mesh comparison of the same 240 astrocytes with a local vessel map at 256 nm (4x finer than the volume-wide map):
- `lumen256.py`: local lumen detection from EM at 256x256x320 nm, with nucleus voxels removed.
- `refine.py`: per astrocyte (m1300 mesh, 45 um radius), membrane 0.15-1.0 um from the detected lumen counts as perivascular. That is about 0.7 um of the true lumen wall, because the detected lumen sits ~0.35 um inside the EM edge. Membrane closer than 0.15 um is flagged as a possible merge with endothelium; none of the 240 had any. Perivascular membrane is grouped into contact sites. Sites touching the cell-body zone (nucleus radius + 3 um) are the soma contact; the rest are endfeet. A lumen-geodesic test says whether an endfoot lies on the same vessel segment as the soma contact.
- `summarize.py`: per-cell table `results/refined/refined_astro_summary.csv` and per-site table `results/refined/refined_sites.csv`.

Refined results:
- Membrane-level soma contact was confirmed in 82 of 120 nucleus-defined vc candidates (68%) and in 4 of 120 typical astrocytes (3%).
- Distal endfeet: confirmed vc-astrocytes median 2 per cell (mean 2.05), typical median 2 (mean 2.00), p = 0.64. 82% of vc vs 66% of typical have at least one (Fisher p = 0.023).
- 86% of vc-astrocyte endfeet are on a different vessel segment than the soma contact. The median distance from the nucleus is 23 um (typical 22.5 um).
- This supersedes the earlier "fewer distal endfeet" result: once soma contact is confirmed at membrane level, vc-astrocytes keep a full endfoot repertoire and add a somatic contact.
- 22-33% of astrocytes had no detected endfoot within 45 um; thin endfeet are lost in the lod-2 meshes, so treat absolute counts as lower bounds.

`scripts/movie/` renders the split-screen movie: `export_cell.py` (mesh classes and vessel surface), `em_inset.py` (EM sections painted with the movie's colours), and `movie.html` + `frames.js` (three.js in headless Chromium, frames encoded with ffmpeg). The exemplars are typical astrocyte nucleus 292597 and vc-astrocyte nucleus 493402 (`results/refined/*_meta.json`).
