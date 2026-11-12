"""Are astrocyte soma contacts and microglial soma contacts spatially exclusive
along individual capillaries?

For each focal contact we take the lumen component it sits on inside a local
box, compute the geodesic (along-lumen) distance field from it, and ask where
the other cell type's contacts fall within that field.

Null: the other contacts are placed uniformly at random on the SAME lumen
component. Because the geodesic field is already computed over every lumen
voxel of that component, the null is exact and needs no simulation: for a
partner at along-vessel distance g,

    u = fraction of that component's lumen voxels closer than g

is Uniform(0,1) under the null. u > 0.5 means the partner sits further away
than a random position would; u < 0.5 means closer.
"""
import numpy as np, pandas as pd, time, sys, json, os
from skimage.graph import MCP_Geometric
from scipy import ndimage as ndi
from scipy import stats
from concurrent.futures import ThreadPoolExecutor

RES = np.array([0.512, 0.512, 0.64])
X0 = np.array([384, 512, 926])
HALF_UM = 60.0                      # local box half-width
MIN_COMP_VOX = 400                  # ignore tiny lumen fragments

LB = np.load('global_lumen_bits.npy', mmap_mode='r')
SH = np.array([LB.shape[0], LB.shape[1], LB.shape[2] * 8])

def get_box(lo, hi):
    zl, zh = lo[2] // 8, (hi[2] + 7) // 8
    sub = np.asarray(LB[lo[0]:hi[0], lo[1]:hi[1], zl:zh])
    un = np.unpackbits(sub, axis=2)
    o = lo[2] - zl * 8
    return un[:, :, o:o + (hi[2] - lo[2])].astype(bool)

C = pd.read_pickle('contacts_all.pkl')
C = C[C.gap < float(sys.argv[1]) if len(sys.argv) > 1 else C.gap < 1.0].copy()
C['base'] = C.kind.str.replace('_verified', '', regex=False)
# one row per nucleus: prefer the verified call when a cell appears in both sets
C = C.sort_values('kind').drop_duplicates('nid', keep='last')
print('contacts used:', C.base.value_counts().to_dict(), flush=True)
LUM = np.c_[C.lx, C.ly, C.lz]
VOX = np.floor(LUM / RES).astype(int) - X0

def one(k):
    r = C.iloc[k]
    c = VOX[k]
    h = (HALF_UM / RES).astype(int)
    lo = np.maximum(c - h, 0); hi = np.minimum(c + h, SH)
    if np.any(hi - lo < 20):
        return None
    L = get_box(lo, hi)
    if not L.any():
        return None
    lab, n = ndi.label(L)
    src = tuple(c - lo)
    if not L[src]:                                  # snap to the nearest lumen voxel
        idx = np.argwhere(L)
        src = tuple(idx[np.argmin(np.linalg.norm((idx - np.array(src)) * RES, axis=1))])
    comp = lab[src]
    M = lab == comp
    nvox = int(M.sum())
    if nvox < MIN_COMP_VOX:
        return None
    cum, _ = MCP_Geometric(np.where(M, 1.0, np.inf), sampling=tuple(RES)).find_costs([src])
    cum = np.where(M, cum, np.inf)
    fin = cum[np.isfinite(cum)]
    if len(fin) < MIN_COMP_VOX:
        return None
    fin_sorted = np.sort(fin)
    # every other contact inside this box that lies on the SAME component
    out = []
    near = np.where((np.abs(VOX[:, 0] - c[0]) < h[0]) & (np.abs(VOX[:, 1] - c[1]) < h[1])
                    & (np.abs(VOX[:, 2] - c[2]) < h[2]))[0]
    for q in near:
        if q == k:
            continue
        v = VOX[q] - lo
        if np.any(v < 0) or np.any(v >= np.array(M.shape)):
            continue
        if not M[tuple(v)]:
            continue                                 # different vessel
        g = cum[tuple(v)]
        if not np.isfinite(g):
            continue
        u = float(np.searchsorted(fin_sorted, g) / len(fin_sorted))
        out.append(dict(focal=int(r.nid), focal_kind=r.base, partner=int(C.iloc[q].nid),
                        partner_kind=C.iloc[q].base, geo=float(g),
                        euc=float(np.linalg.norm((VOX[q] - c) * RES)), u=u,
                        comp_len_vox=nvox, comp_max_geo=float(fin.max())))
    return out

t0 = time.time()
rows = []
with ThreadPoolExecutor(4) as ex:
    for i, res in enumerate(ex.map(one, range(len(C)))):
        if res:
            rows.extend(res)
        if i % 100 == 0:
            print(i, len(C), round(time.time() - t0), 's  pairs', len(rows), flush=True)
P = pd.DataFrame(rows)
P.to_csv('exclusivity_pairs.csv', index=False)
print('\npairs found:', len(P), flush=True)
if len(P):
    print(P.groupby(['focal_kind', 'partner_kind']).size().to_string(), flush=True)
print('DONE', flush=True)
