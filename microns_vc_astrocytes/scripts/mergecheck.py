import json,numpy as np,pandas as pd
from concurrent.futures import ThreadPoolExecutor
from scipy.spatial import cKDTree
from astromesh import root_at
R=[r for r in json.load(open('meshbatch.json')) if 'err' not in r]
df=pd.read_pickle('nuclei_L.pkl')
X=np.c_[df.x*0.004,df.y*0.004,df.z*0.04]; T=cKDTree(X)
def check(r):
    s=df[df.nid==r['nid']].iloc[0]; c=np.array([s.x*0.004,s.y*0.004,s.z*0.04])
    idx=[i for i in T.query_ball_point(c,60) if df.nid.iat[i]!=r['nid'] and df.dl_center.iat[i]<3 and df.z.iat[i]<27850]
    same=[]
    for i in idx:
        q=df.iloc[i]
        try: rr=root_at(q.x,q.y,q.z)
        except Exception: continue
        if rr==r['root']: same.append((int(q.nid),q.grp,float(q.dl_center)))
    return r['nid'],same
with ThreadPoolExecutor(8) as ex: res=list(ex.map(check,R))
json.dump(res,open('mergecheck.json','w'))
n_any=sum(1 for _,s in res if s); n_vasc=sum(1 for _,s in res if any(d<3 for _,_,d in s))
print('astros',len(res),'share root with another nucleus',n_any,'with a vessel-adjacent nucleus',n_vasc)
from collections import Counter
print(Counter(g for _,s in res for _,g,_ in s).most_common(12))
