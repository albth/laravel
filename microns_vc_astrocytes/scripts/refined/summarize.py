import json,glob,numpy as np,pandas as pd
rows=[];sites=[]
for f in glob.glob('refine_out/*.json'):
    r=json.load(open(f))
    if 'err' in r: print('ERR',r['nid'],r['err'][:120]); continue
    S=r['sites']; d=[s for s in S if not s['soma']]
    soma_comps={s['comp'] for s in S if s['soma']}
    for s in S: sites.append(dict(nid=r['nid'],cls=r['cls'],**{k:v for k,v in s.items() if k not in('pos','lumen_vox')}))
    rows.append(dict(nid=r['nid'],root=r['root'],cls=r['cls'],area_total=r['area_total'],area_peri=r['area_peri'],lum_frac=r['area_luminal']/max(r['area_peri']+r['area_luminal'],1e-9),
        soma_peri=r['soma_peri'],soma_lum=r['soma_luminal'],n_sites=len(S),n_distal=len(d),soma_site=any(s['soma'] for s in S),
        n_distal_samecomp=sum(s['comp'] in soma_comps for s in d),
        n_distal_sameseg=sum(('geo' in s) and np.isfinite(s['geo']) and s['geo']<1.5*s['eu'] for s in d),
        distal_dist=np.median([s['dmin_soma'] for s in d]) if d else np.nan, distal_area=np.median([s['area'] for s in d]) if d else np.nan))
D=pd.DataFrame(rows); X=pd.DataFrame(sites)
D.to_csv('refined_astro_summary.csv',index=False); X.to_csv('refined_sites.csv',index=False)
print(len(D)); print(D.groupby('cls')[['area_total','area_peri','lum_frac','soma_peri','n_sites','n_distal','distal_dist','distal_area']].median().round(2).T)
print(D.groupby('cls').apply(lambda x: pd.Series({'soma_site':x.soma_site.mean(),'no_site':(x.n_sites==0).mean(),'distal>=1':(x.n_distal>=1).mean(),'lumfrac>0.3':(x.lum_frac>0.3).mean()})).round(3))
