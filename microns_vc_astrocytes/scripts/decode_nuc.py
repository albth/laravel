import re,json,base64,numpy as np,pandas as pd
s=open('../microns-neuroglancer-tool/ujump.html').read()
def js(id_):
    m=re.search(r'<script type="application/json" id="'+id_+r'">(.*?)</script>',s,re.S); return json.loads(m.group(1))
def b(x,dt): return np.frombuffer(base64.b64decode(x),dtype=dt)
D=js('nucdata')
df=pd.DataFrame(dict(x=b(D['XB'],'<u4'),y=b(D['YB'],'<u4'),z=b(D['ZB'],'<u2'),t=b(D['TB'],'u1'),h=b(D['HB'],'u1'),lo=b(D['LB'],'<u4'),nid=b(D['NB'],'<u4'),vol=b(D['VB'],'<f4')))
ct=D['CT_NAMES']+['?']*300
df['ctype']=['none' if i==0 else D['CT_NAMES'][i-1] for i in df.t]
HI=D['HI']
df['root_id']=[ (HI[h]*2**32+int(l)) if h!=255 else 0 for h,l in zip(df.h,df.lo)]
V=js('vascdist12')
dist=b(V['OV_DIST_ALL'],'<u2').astype(float); dist[dist==65535]=np.nan
df['vdist_tool_um']=dist[:len(df)]*V['UNIT_NM']/1000
O=js('owndata'); ov=b(O['OVERRIDE']['TYPE'],'u1'); names=O['TYPE_NAMES']
df['verified']=[names[i] if i!=255 else '' for i in ov]
print(df.head()); print(df.ctype.value_counts()); print(df.verified.value_counts())
print(df[['x','y','z']].describe())
df.to_parquet('nuclei.parquet') if False else df.to_pickle('nuclei.pkl')
