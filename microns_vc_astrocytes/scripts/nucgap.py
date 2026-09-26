import numpy as np,pandas as pd,pickle,time,sys
from cloudvolume import CloudVolume
from scipy import ndimage as ndi
from concurrent.futures import ThreadPoolExecutor
from lumenlib import RES,load_block,B
T=pickle.load(open('lumen_tree.pkl','rb'))
df=pd.read_pickle('nuclei_L.pkl')
rng=np.random.default_rng(3)
pick=[]
for g,n in [('astrocyte',None),('V:Astrocyte',None),('microglia',None),('V:Microglia',None),('OPC',None),('oligo',2500),('23P',2000),('4P',1500),('5P-IT',1000),('6P-IT',1000),('BC',1500),('pericyte',500),('V:Endothelial cell',500)]:
    s=df[df.grp==g]
    if n and len(s)>n: s=s.sample(n,random_state=4)
    pick.append(s)
sub=pd.concat(pick); print(len(sub),flush=True)
nv=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/nuclei',mip=[256,256,320],use_https=True,progress=False,fill_missing=True)
NR=np.array([0.256,0.256,0.32])
def fetch(r):
    c=(np.array([r.x*4,r.y*4,r.z*40])/1000/NR).astype(int); h=(10/NR).astype(int)
    for att in range(4):
        try:
            a=np.asarray(nv[c[0]-h[0]:c[0]+h[0],c[1]-h[1]:c[1]+h[1],c[2]-h[2]:c[2]+h[2]])[...,0]; break
        except Exception: time.sleep(2**att)
    else: return None
    m=a==r.nid
    if m.sum()<10: return None
    er=ndi.binary_erosion(m); S=np.argwhere(m&~er)
    pts=(S+c-h+0.5)*NR
    allp=(np.argwhere(m)+c-h+0.5)*NR
    cen=allp.mean(0); cov=np.cov((allp-cen).T); ev,evec=np.linalg.eigh(cov)
    return pts,cen,np.sqrt(np.maximum(ev,0)),evec
X0=np.array([384,512,926])
def inlumen_old(P):
    vox=np.floor(P/RES).astype(int); bidx=(vox-X0)//B; res=np.zeros(len(P),bool); ok=np.zeros(len(P),bool)
    for key in np.unique(bidx,axis=0):
        sel=np.where(np.all(bidx==key,axis=1))[0]; o=X0+B*key
        try: _,m,vv=load_block(f'blocks/{o[0]}_{o[1]}_{o[2]}.npz')
        except FileNotFoundError: continue
        l=vox[sel]-o; g=np.all((l>=0)&(l<np.array(m.shape)),axis=1); s=sel[g]; l=l[g]
        res[s]=m[l[:,0],l[:,1],l[:,2]]; ok[s]=vv[l[:,0],l[:,1],l[:,2]]
    return res,ok
from lumenlib import lookup as inlumen
import os,glob as _g
K=20
t=time.time()
CH=2000
for c0 in range(0,len(sub),CH):
  fnc=f'nucgap_parts/part_{c0:06d}.pkl'
  if os.path.exists(fnc): continue
  part=sub.iloc[c0:c0+CH]; rows=[]
  with ThreadPoolExecutor(8) as ex:
    for k,(res_,(_,r)) in enumerate(zip(ex.map(fetch,[r for _,r in part.iterrows()]),part.iterrows())):
          if res_ is None: rows.append(dict(idx=r.name,ok=False)); continue
          pts,cen,sd,evec=res_
          if len(pts)>300: pts=pts[rng.choice(len(pts),300,replace=False)]
          d,j=T.query(pts); gap=d.min()
          # nearest lumen direction vs nucleus long/short axes
          near=T.data[j[d.argmin()]]; u=(near-cen); u/=np.linalg.norm(u)+1e-9
          cos_short=abs(u@evec[:,0])
          # null: translate shape
          v=rng.normal(size=(K,3)); v/=np.linalg.norm(v,axis=1)[:,None]; rr=(10**3+(40**3-10**3)*rng.random(K))**(1/3)
          sh=v*rr[:,None]
          cins,cok=inlumen(cen+sh)
          ng=[]
          for q in range(K):
              if cins[q] or not cok[q]: ng.append(np.nan); continue
              il,_=inlumen(pts+sh[q])
              if il.any(): ng.append(np.nan); continue
              ng.append(T.query(pts+sh[q])[0].min())
          il0,_=inlumen(pts); rows.append(dict(idx=r.name,ok=True,gap=gap,overlap=float(il0.mean()),sd1=sd[0],sd2=sd[1],sd3=sd[2],cos_short=cos_short,null=ng))
    pd.DataFrame(rows).to_pickle(fnc); print(c0,round(time.time()-t),flush=True)
R=pd.concat([pd.read_pickle(f) for f in sorted(_g.glob('nucgap_parts/*.pkl'))]).set_index('idx')
out=sub.join(R); out.to_pickle('nucgap.pkl'); print('DONE',flush=True)
