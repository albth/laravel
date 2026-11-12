# Are astrocyte and microglial soma contacts exclusive along capillaries?

Test of hypothesis 5 from the soma-vs-endfoot schematic: a 525 µm² astrocyte
somatic sheath occupies a lot of capillary surface, and microglia are known to
reach the vasculature through gaps between astrocyte endfeet, so the two
contact types might avoid each other.

**They do not. The two are positioned independently along shared capillaries.**

## Method

1. `contacts.py` — nucleus-surface-to-lumen gap for every MICrONS-predicted and
   hand-verified astrocyte and microglion (10,488 cells screened, 7,356 with a
   nucleus centre within 15 µm of a lumen, 7,305 measured against the nucleus
   segmentation). Contact rates reproduce the earlier analysis:

   | population | n | gap < 1 µm | gap < 2 µm |
   |---|---|---|---|
   | astrocyte (predicted) | 7,390 | 5.7% | 10.6% |
   | astrocyte (verified) | 1,497 | 9.6% | 22.4% |
   | microglia (predicted) | 2,044 | 13.0% | 19.6% |
   | microglia (verified) | 863 | 11.4% | 19.8% |

2. `exclusivity.py` — for each contact, isolate the lumen connected component it
   sits on inside a 120 µm box, compute the geodesic (along-lumen) distance
   field from it with `MCP_Geometric`, and locate every other contact on that
   same component.

3. **Exact null, no simulation.** Because the geodesic field is known at every
   lumen voxel of the component, a partner at along-vessel distance *g* is
   scored as

   ```
   u = (lumen voxels of that component with geodesic distance < g) / (all its lumen voxels)
   ```

   Under random placement along the same vessel, `u ~ Uniform(0,1)`. The null
   therefore conditions on the vessel's own shape, length and contact count.
   `u > 0.5` means the partner is further away than a random position;
   `u < 0.5` means closer. A physical exclusion zone empties the lowest decile.

## Result (gap < 1 µm; 563 astrocyte and 363 microglial contacts, 1,411 pairs)

| comparison | pairs | mean u (95% CI) | nearest decile | p |
|---|---|---|---|---|
| astrocyte → microglia | 298 | 0.502 (0.466–0.538) | 13.1% | 0.08 |
| microglia → astrocyte | 301 | 0.519 (0.480–0.558) | 10.3% | 0.85 |
| astrocyte → astrocyte *(control)* | 554 | 0.468 (0.438–0.499) | **13.5%** | **0.007** |
| microglia → microglia *(control)* | 258 | 0.523 (0.484–0.561) | 7.4% | 0.18 |

Nearest decile = % of partners in the closest tenth of the shared vessel
(about 30 µm), against 10% expected; p from an exact binomial test.

* **No exclusion, at any scale tested.** The nearest decile is never depleted
  for the heterotypic comparison — if anything it runs slightly above chance.
  The ten closest astrocyte/microglia pairs are 1.7–10.6 µm apart along the
  vessel.
* **Power.** The CI on mean u is ±0.036, which on a median 203 µm component
  excludes a systematic displacement larger than about 7 µm. A sheath spanning
  tens of µm of vessel would have been detected easily.
* **Robustness.** The same conclusion holds at a 2 µm contact threshold
  (`pairs_gap2um.csv`; astrocyte → microglia nearest decile 11.5%, p = 0.14)
  and with the hand-verified sets dropped, in case curation clustered them
  (13.8%, p = 0.09).

## Unexpected: astrocyte soma contacts cluster with each other

The control comparison is the positive finding. Astrocyte contacts sit closer
together along a shared vessel than random placement allows:

| threshold | cells | nearest decile | p |
|---|---|---|---|
| gap < 1 µm, all | 554 pairs | 13.5% | 0.007 |
| gap < 2 µm, all | 1,920 pairs | 15.5% | 8.7 × 10⁻¹⁴ |
| gap < 1 µm, predicted only | 325 pairs | 16.0% | 0.0008 |
| gap < 2 µm, predicted only | 972 pairs | 14.8% | 2.3 × 10⁻⁶ |

The null already conditions on how many contacts that capillary carries, so
this is a within-vessel effect, not "some vessels attract astrocytes". Two
readings to separate: particular segments are attractive (branch point,
pericyte soma, thinner basal lamina), or contacts are laid down together,
which would point at a developmental or clonal origin given that astrocytes
migrate along vessels before settling. The first is testable in this volume.

Microglial contacts show a weak opposite trend (nearest decile 7.4% / 6.7%,
p = 0.18 / 0.04) — possible mild self-avoidance, not firm.

## Caveats

* The 120 µm box truncates each vessel component. Observed and null use the
  same truncated component, so the comparison is fair, but partners further
  than the box are not counted — the test is about relative position given
  both contacts are present on a shared segment.
* Sampling the null by lumen voxel is uniform in vessel *volume*, which equals
  uniform in arc length only for constant-calibre vessels. Measured calibre at
  these contacts is 2.5–2.8 µm, so the approximation is close.
* Cell types are MICrONS model predictions except where marked verified; the
  predicted-only control above addresses annotation bias, not classifier error.

## Files

* `contacts.py`, `exclusivity.py` — the pipeline (needs the rebuilt 512 nm
  lumen map from `../scripts/`).
* `soma_contacts.csv` — all 1,693 contacts with gap < 2 µm: nucleus id, cell
  type, coordinates, gap, nucleus-surface and lumen points.
* `pairs_gap1um.csv`, `pairs_gap2um.csv` — every scored pair with its geodesic
  distance, straight-line distance and `u`.
