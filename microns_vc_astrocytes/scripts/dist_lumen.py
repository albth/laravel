import numpy as np,pandas as pd,pickle
from scipy.spatial import cKDTree
from lumenlib import RES
S=np.load('lumen_surface_vox6.npy')
P=(S+0.5)*RES
T=cKDTree(P,leafsize=32,balanced_tree=False,compact_nodes=False)
pickle.dump(T,open('lumen_tree.pkl','wb'),protocol=4)
df=pd.read_pickle('nuclei_d.pkl')
X=np.c_[df.x*0.004,df.y*0.004,df.z*0.04]
d,i=T.query(X,workers=4)
df['dl_center']=d
df['r_eq']=(3*df.vol.clip(lower=1)/(4*np.pi))**(1/3)
df['gap_eq']=d-df.r_eq
# sign: is nucleus center inside lumen? check via nearest surface point direction is ambiguous; flag centers inside mask later
df['lum_nx'],df['lum_ny'],df['lum_nz']=P[i].T
df['grp']=[('V:'+v) if v else c for v,c in zip(df.verified,df.ctype)]
df.to_pickle('nuclei_L.pkl')
g=df.groupby('grp')
print(pd.DataFrame({'n':g.size(),'med_dl':g.dl_center.median(),'med_gap':g.gap_eq.median(),'p_gap<2':g.gap_eq.apply(lambda s:(s<2).mean()),'p_gap<1':g.gap_eq.apply(lambda s:(s<1).mean()),'r_eq':g.r_eq.median()}).sort_values('p_gap<2').round(3).to_string())
