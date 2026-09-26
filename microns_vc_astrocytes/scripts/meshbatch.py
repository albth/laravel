import numpy as np,pandas as pd,glob,json,os,sys,time
from concurrent.futures import ThreadPoolExecutor
from astromesh import analyze
parts=pd.concat([pd.read_pickle(f) for f in sorted(glob.glob('nucgap_parts/*.pkl'))]).set_index('idx')
df=pd.read_pickle('nuclei_L.pkl').join(parts,how='inner')
svm=pd.read_csv('../mic/nucleus_neuron_svm.csv',header=None,names=['id','valid','sv','root','x','y','z','tag','cls'])[['id','cls']]
df=df.merge(svm,left_on='nid',right_on='id',how='left').set_index(df.index)
a=df[(df.grp=='astrocyte')&(df.cls=='not-neuron')&df.ok.astype(bool)]
vc=a[a.gap<1.0].sample(min(120,(a.gap<1.0).sum()),random_state=5)
far=a[a.gap>6].sample(120,random_state=5)
sel=pd.concat([vc.assign(cls_vc='vc'),far.assign(cls_vc='nonvc')])
print(len(vc),len(far),flush=True)
def run(item):
    i,r=item
    try: out=analyze(r)
    except Exception as e: return dict(idx=i,err=str(e))
    if out is None: return dict(idx=i,err='noroot')
    out['idx']=i; out['cls_vc']=r.cls_vc; out['gap']=r.gap; out['nid']=int(r.nid)
    return out
res=[]
with ThreadPoolExecutor(4) as ex:
    for k,o in enumerate(ex.map(run,list(sel.iterrows()))):
        res.append(o)
        if k%20==0: print(k,flush=True)
json.dump(res,open('meshbatch.json','w'),default=float); print('DONE',flush=True)
