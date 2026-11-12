"""Find soma-vessel contacts for astrocytes and microglia.

Stage 1: cheap screen on nucleus-centre distance to the lumen surface.
Stage 2: for survivors, the true nucleus-surface-to-lumen gap from the
nucleus segmentation (same measure used in the earlier analysis).
"""
import numpy as np, pandas as pd, pickle, time, os, sys, json
from scipy.spatial import cKDTree
from scipy import ndimage as ndi
from cloudvolume import CloudVolume
from concurrent.futures import ThreadPoolExecutor

RES = np.array([0.512, 0.512, 0.64])          # mip6 voxel size, um
NR = np.array([0.256, 0.256, 0.32])           # nucleus segmentation voxel size

df = pd.read_pickle('nuclei.pkl')
df['grp'] = [('V:' + v) if v else c for v, c in zip(df.verified, df.ctype)]
# cells of interest: MICrONS-predicted astrocytes and microglia, plus the
# hand-verified sets, kept separate so we can check they agree
WANT = {'astrocyte': 'astrocyte', 'microglia': 'microglia',
        'V:Astrocyte': 'astrocyte_verified', 'V:Microglia': 'microglia_verified'}
sub = df[df.grp.isin(WANT)].copy()
sub['kind'] = sub.grp.map(WANT)
print('candidate nuclei:', sub.kind.value_counts().to_dict(), flush=True)

S = np.load('lumen_surface_vox6.npy')
P = (S + 0.5) * RES
T = cKDTree(P, leafsize=32, balanced_tree=False, compact_nodes=False)
X = np.c_[sub.x * 0.004, sub.y * 0.004, sub.z * 0.04]
d, j = T.query(X, workers=4)
sub['d_centre'] = d
sub['lum_pt'] = list(P[j])
print('with centre within 15 um of a lumen:', (d < 15).sum(), flush=True)

cand = sub[sub.d_centre < 15].copy()
NV = CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/nuclei',
                 mip=[256, 256, 320], use_https=True, progress=False, fill_missing=True)

def gap(row):
    c = (np.array([row.x * 0.004, row.y * 0.004, row.z * 0.04]) / NR).astype(int)
    h = (12 / NR).astype(int)
    for att in range(4):
        try:
            a = np.asarray(NV[c[0]-h[0]:c[0]+h[0], c[1]-h[1]:c[1]+h[1], c[2]-h[2]:c[2]+h[2]])[..., 0]
            break
        except Exception:
            time.sleep(2 ** att)
    else:
        return None
    m = a == row.nid
    if m.sum() < 10:
        return None
    er = ndi.binary_erosion(m)
    pts = (np.argwhere(m & ~er) + c - h + 0.5) * NR
    if len(pts) > 400:
        pts = pts[np.random.default_rng(0).choice(len(pts), 400, replace=False)]
    dd, jj = T.query(pts)
    k = int(np.argmin(dd))
    return float(dd[k]), pts[k].tolist(), P[jj[k]].tolist()

rows = []
t0 = time.time()
with ThreadPoolExecutor(6) as ex:
    for n, (res, (i, r)) in enumerate(zip(ex.map(gap, [r for _, r in cand.iterrows()]), cand.iterrows())):
        if res is None:
            continue
        g, surf, lum = res
        rows.append(dict(idx=i, nid=int(r.nid), kind=r.kind, grp=r.grp,
                         x=int(r.x), y=int(r.y), z=int(r.z), gap=g,
                         sx=surf[0], sy=surf[1], sz=surf[2],
                         lx=lum[0], ly=lum[1], lz=lum[2]))
        if n % 500 == 0:
            print(n, len(cand), round(time.time() - t0), flush=True)
R = pd.DataFrame(rows)
R.to_pickle('contacts_all.pkl')
print('\nmeasured:', len(R), flush=True)
for k, g in R.groupby('kind'):
    tot = (sub.kind == k).sum()
    print(f'{k:<22} n={tot:5d}  gap<1um {(g.gap<1).sum():4d} ({100*(g.gap<1).sum()/tot:.1f}%)'
          f'  gap<2um {(g.gap<2).sum():4d} ({100*(g.gap<2).sum()/tot:.1f}%)', flush=True)
print('DONE', flush=True)
