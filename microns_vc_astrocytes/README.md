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
