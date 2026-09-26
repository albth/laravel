import numpy as np,pandas as pd,pickle
from lumenlib import RES,load_block,B
T=pickle.load(open('lumen_tree.pkl','rb'))
df=pd.read_pickle('nuclei_L.pkl')
groups=['astrocyte','V:Astrocyte','microglia','V:Microglia','oligo','OPC','23P','4P','5P-IT','5P-ET','6P-IT','6P-CT','BC','MC','BPC','NGC','pericyte','V:Endothelial cell']
sub=df[df.grp.isin(groups)].copy()
X=np.c_[sub.x*0.004,sub.y*0.004,sub.z*0.04]
K=30; rng=np.random.default_rng(1)
idx=np.repeat(np.arange(len(sub)),K)
v=rng.normal(size=(len(idx),3)); v/=np.linalg.norm(v,axis=1)[:,None]
rr=(10**3+(40**3-10**3)*rng.random(len(idx)))**(1/3)   # uniform in 10-40um shell
R=X[idx]+v*rr[:,None]
vox=np.floor(R/RES).astype(int); X0=np.array([384,512,926]); bidx=(vox-X0)//B
ok=np.zeros(len(R),bool)
keys,inv=np.unique(bidx,axis=0,return_inverse=True)
for k,key in enumerate(keys):
    sel=np.where(inv.ravel()==k)[0]; o=X0+B*key
    try: _,m,vv=load_block(f'blocks/{o[0]}_{o[1]}_{o[2]}.npz')
    except FileNotFoundError: continue
    l=vox[sel]-o; good=np.all((l>=0)&(l<np.array(m.shape)),axis=1)
    s=sel[good]; l=l[good]
    ok[s]=vv[l[:,0],l[:,1],l[:,2]]&~m[l[:,0],l[:,1],l[:,2]]
d,_=T.query(R,workers=4); d[~ok]=np.nan
D=d.reshape(len(sub),K)
np.save('localnull_D.npy',D); sub.to_pickle('localnull_sub.pkl')
r=sub.r_eq.values[:,None]
rows=[]
for G in (1.0,1.5,2.0,3.0):
    valid=(D>r)
    e=np.nansum(((D>r)&(D<r+G)),axis=1)/np.maximum(valid.sum(1),1)
    sub[f'exp{G}']=e; sub[f'obs{G}']=(sub.gap_eq<G).astype(float)
g=sub.groupby('grp')
out=pd.DataFrame({'n':g.size()})
for G in (1.0,1.5,2.0,3.0):
    out[f'obs{G}']=g[f'obs{G}'].mean(); out[f'exp{G}']=g[f'exp{G}'].mean(); out[f'enr{G}']=out[f'obs{G}']/out[f'exp{G}']
print(out.round(3).sort_values('enr1.5').to_string()); out.to_csv('enrichment_local.csv'); sub.to_pickle('localnull_sub.pkl')
