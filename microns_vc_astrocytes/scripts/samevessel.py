import json,numpy as np,pandas as pd
from skimage.graph import MCP_Geometric
from scipy.spatial import cKDTree
from lumenlib import RES,get_region
R=[r for r in json.load(open('meshbatch.json')) if 'err' not in r]
VP=np.load('vessel_pts.npy',allow_pickle=True).item()
trees={k:cKDTree(v) for k,v in VP.items()}
def vtype(p):
    best=(1e9,None)
    for k,t in trees.items():
        d,_=t.query(p); best=min(best,(d,k))
    return best[1] if best[0]<6 else 'untraced'
df=pd.read_pickle('nuclei_L.pkl')
out=[]
for r in R:
    s=df[df.nid==r['nid']].iloc[0]; cen=np.array([s.x*0.004,s.y*0.004,s.z*0.04])
    P=[p for p in r['patches'] if p['area']>=5]
    for p in P: out.append(dict(nid=r['nid'],cls=r['cls_vc'],soma=p['soma'],area=p['area'],dist=p['dmin_soma'],caliber=p['caliber'],vtype=vtype(np.array(p['lumen_pt']))))
    if r['cls_vc']!='vc': continue
    sp=[p for p in P if p['soma']]; dp=[p for p in P if not p['soma']]
    if not sp or not dp: continue
    src=np.array(sp[0]['lumen_pt']); lo=np.floor((cen-60)/RES).astype(int); hi=np.ceil((cen+60)/RES).astype(int)
    L,_=get_region(lo,hi)
    cost=np.where(L,1.0,np.inf)
    m=MCP_Geometric(cost,sampling=tuple(RES))
    si=tuple(np.clip(np.floor(src/RES).astype(int)-lo,0,np.array(L.shape)-1))
    if not L[si]:
        # snap to nearest lumen voxel
        idx=np.argwhere(L); si=tuple(idx[np.argmin(np.linalg.norm((idx+lo)*RES-src,axis=1))])
    cum,_=m.find_costs([si])
    idx=np.argwhere(L)
    for p in dp:
        t=np.array(p['lumen_pt']); j=np.argmin(np.linalg.norm((idx+lo+0.5)*RES-t,axis=1)); g=cum[tuple(idx[j])]
        eu=np.linalg.norm(t-src)
        out.append(dict(nid=r['nid'],cls='vc_pair',geo=float(g),eu=float(eu),dist=p['dmin_soma']))
O=pd.DataFrame(out); O.to_pickle('sites.pkl')
S=O[O.cls!='vc_pair']
print(pd.crosstab([S.cls,S.soma],S.vtype))
pr=O[O.cls=='vc_pair']; pr['ratio']=pr.geo/pr.eu
print(len(pr), 'geo finite',np.isfinite(pr.geo).mean())
print(pr[['geo','eu','ratio']].replace(np.inf,np.nan).describe())
print('same vessel segment (geo<1.5*eu and finite, within 60um box):',((pr.ratio<1.5)&np.isfinite(pr.geo)).mean())
