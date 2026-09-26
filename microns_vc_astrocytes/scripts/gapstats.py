import numpy as np,pandas as pd,glob,json
parts=pd.concat([pd.read_pickle(f) for f in sorted(glob.glob('nucgap_parts/*.pkl'))]).set_index('idx')
df=pd.read_pickle('nuclei_L.pkl').join(parts,how='inner'); df=df[df.ok.astype(bool)]
df.to_pickle('nucgap_partial.pkl')
rng=np.random.default_rng(0)
rows=[]
for g,s in df.groupby('grp'):
    if len(s)<200: continue
    N=np.array(s.null.tolist(),float)
    r={'group':g,'n':len(s),'overlap_frac':(s.overlap>0).mean()}
    for G in (0.5,1.0,1.5,2.0,3.0):
        obs=(s.gap<G).mean(); exp=np.nanmean(N<G)
        # bootstrap CI of enrichment over cells
        bs=[]
        for _ in range(300):
            ii=rng.integers(0,len(s),len(s)); bs.append((s.gap.values[ii]<G).mean()/np.nanmean(N[ii]<G))
        r[f'obs{G}']=obs; r[f'exp{G}']=exp; r[f'enr{G}']=obs/exp; r[f'lo{G}'],r[f'hi{G}']=np.percentile(bs,[2.5,97.5])
    rows.append(r)
R=pd.DataFrame(rows).set_index('group').sort_values('enr1.0')
print(R[['n','overlap_frac','obs1.0','exp1.0','enr1.0','lo1.0','hi1.0','obs2.0','exp2.0','enr2.0']].round(3).to_string())
R.to_csv('gap_enrichment.csv')
