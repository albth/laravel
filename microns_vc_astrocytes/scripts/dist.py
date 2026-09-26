import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from vessels import sample,use
df=pd.read_pickle('nuclei.pkl')
P=np.c_[df.x*0.004,df.y*0.004,df.z*0.04]
maj=['pial_art','pen_art','pial_ven','asc_ven']; cap=['acap1','acap2','acap3','vcap1','vcap2','vcap3']
allpts={}
for k in use:
    pts,_=sample(k); allpts[k]=pts
    d,_=cKDTree(pts).query(P); df['d_'+k]=d
df['d_major']=df[['d_'+k for k in maj]].min(axis=1)
df['d_cap']=df[['d_'+k for k in cap]].min(axis=1)
df['d_any']=df[['d_'+k for k in use]].min(axis=1)
df['near_type']=df[['d_'+k for k in use]].idxmin(axis=1).str[2:]
np.save('vessel_pts.npy',allpts,allow_pickle=True)
df.to_pickle('nuclei_d.pkl')
def grp(r):
    if r.verified: return 'V:'+r.verified
    return r.ctype
df['grp']=[('V:'+v) if v else c for v,c in zip(df.verified,df.ctype)]
g=df.groupby('grp')
print(pd.DataFrame({'n':g.size(),'med_dany':g.d_any.median(),'p_lt5':g.d_any.apply(lambda s:(s<5).mean()),'p_lt10':g.d_any.apply(lambda s:(s<10).mean()),'tool_med':g.vdist_tool_um.median(),'corr_tool':g.apply(lambda x:np.corrcoef(x.d_any,x.vdist_tool_um)[0,1])}).sort_values('p_lt5').to_string())
