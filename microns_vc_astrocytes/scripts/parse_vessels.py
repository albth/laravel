import re,glob,numpy as np,json,collections
num=r'(-?[\d.eE+-]+)'
pat=re.compile(r'"pointA":\s*\[\s*'+num+r',\s*'+num+r',\s*'+num+r'\s*\],\s*"pointB":\s*\[\s*'+num+r',\s*'+num+r',\s*'+num+r'\s*\]'
               r'|"point":\s*\[\s*'+num+r',\s*'+num+r',\s*'+num+r'\s*\]'
               r'|"name":\s*"([^"]*)"')
out={}
for f in sorted(glob.glob('*.txt')):
    s=open(f).read()
    pending_l=[];pending_p=[]
    for m in pat.finditer(s):
        g=m.groups()
        if g[0] is not None: pending_l.append([float(x) for x in g[0:6]])
        elif g[6] is not None: pending_p.append([float(x) for x in g[6:9]])
        else:
            name=g[9]
            if pending_l or pending_p:
                key=f+' :: '+name
                out.setdefault(key,{'lines':[],'points':[]})
                out[key]['lines']+=pending_l; out[key]['points']+=pending_p
            pending_l=[];pending_p=[]
    if pending_l or pending_p:
        out[f+' :: <unnamed>']={'lines':pending_l,'points':pending_p}
res={}
for k,v in out.items():
    L=np.array(v['lines']).reshape(-1,6); P=np.array(v['points']).reshape(-1,3)
    print(f'{k:60s} lines={len(L):6d} points={len(P):5d}', 'x',L[:,0].min() if len(L) else '', L[:,0].max() if len(L) else '','y',L[:,1].min() if len(L) else '',L[:,1].max() if len(L) else '','z',L[:,2].min() if len(L) else '',L[:,2].max() if len(L) else '')
    res[k]=(L,P)
np.save('../work/vessel_layers.npy',res,allow_pickle=True)
