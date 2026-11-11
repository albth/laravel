import numpy as np,pandas as pd,json,time,os,sys
from cloudvolume import CloudVolume
from scipy import ndimage as ndi
from skimage.graph import MCP_Geometric
from concurrent.futures import ThreadPoolExecutor
from lumen256 import lumen_box,RES
SEG=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/seg_m1300',use_https=True,progress=False)
nd=pd.read_csv('../mic/nucleus_detection_v0.csv',header=None,names=['id','valid','sv','root','x','y','z','vol']).set_index('id')
old=[r for r in json.load(open('old_meshbatch.json')) if 'err' not in r]
os.makedirs('refine_out',exist_ok=True)
R=45.0
def run(o):
    fn=f"refine_out/{o['nid']}.json"
    if os.path.exists(fn): return
    try:
        q=nd.loc[o['nid']]; cen=np.array([q.x*0.004,q.y*0.004,q.z*0.04]); rn=(3*q.vol/(4*np.pi))**(1/3)
        m=SEG.mesh.get(o['root'],lod=2); m=list(m.values())[0] if isinstance(m,dict) else m
        V=m.vertices/1000.0; F=m.faces
        keep=np.linalg.norm(V-cen,axis=1)<R; F=F[keep[F].all(1)]
        a=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
        L,lo,valid,_=lumen_box(cen,R+4)
        dout=ndi.distance_transform_edt(~L,sampling=RES).astype(np.float32); din=ndi.distance_transform_edt(L,sampling=RES).astype(np.float32)
        def samp(f,P):
            ij=np.floor(P/RES).astype(int)-lo; ij=np.clip(ij,0,np.array(f.shape)-1); return f[ij[:,0],ij[:,1],ij[:,2]]
        vd=samp(dout,V); fd=vd[F].mean(1); fc=V[F].mean(1)
        lumi=fd<0.15; peri=(fd>=0.15)&(fd<=1.0)
        dc=np.linalg.norm(fc-cen,axis=1); somaz=dc<rn+3
        out=dict(nid=int(o['nid']),root=int(o['root']),cls=o['cls_vc'],gap_old=o['gap'],rn=float(rn),
                 area_total=float(a.sum()),area_peri=float(a[peri].sum()),area_luminal=float(a[lumi].sum()),
                 soma_peri=float(a[peri&somaz].sum()),soma_luminal=float(a[lumi&somaz].sum()),lumen_frac=float(L.mean()))
        lab,nl=ndi.label(L)
        sites=[]
        pf=np.where(peri)[0]
        if len(pf):
            g=np.floor((fc[pf]-cen+R+2)/1.0).astype(int); gs=g.max(0)+2
            grid=np.zeros(gs,bool); grid[g[:,0],g[:,1],g[:,2]]=True
            gl,nc=ndi.label(ndi.binary_dilation(grid,structure=np.ones((3,3,3))),structure=np.ones((3,3,3)))
            fl=gl[g[:,0],g[:,1],g[:,2]]
            Lidx=None
            for k in range(1,nc+1):
                s=pf[fl==k]; ar=float(a[s].sum())
                if ar<3.0: continue
                pc=fc[s].mean(0); dmin=float(dc[s].min())
                # nearest lumen voxel to the patch face closest to patch centroid
                ps=fc[s][np.argmin(np.linalg.norm(fc[s]-pc,axis=1))]
                ij=np.floor(ps/RES).astype(int)-lo; w=np.ceil(3/RES).astype(int)
                sl=tuple(slice(max(0,ij[t]-w[t]),min(L.shape[t],ij[t]+w[t]+1)) for t in range(3))
                sub=L[sl]
                if not sub.any(): continue
                idx=np.argwhere(sub)+np.array([x.start for x in sl]); P=(idx+lo+0.5)*RES
                j=np.argmin(np.linalg.norm(P-ps,axis=1)); lp=idx[j]
                near=np.linalg.norm(P-P[j],axis=1)<2.5
                cal=2*float(din[tuple(idx[near].T)].max())
                sites.append(dict(area=ar,dmin_soma=dmin,soma=bool(dmin<rn+3),pos=ps.tolist(),lumen_vox=lp.tolist(),comp=int(lab[tuple(lp)]),caliber=cal))
        # same-vessel geodesic from soma site
        ss=[s for s in sites if s['soma']]
        if ss and len(sites)>1:
            src=tuple(max(ss,key=lambda s:s['area'])['lumen_vox'])
            cost=np.where(L,1.0,np.inf); cum,_=MCP_Geometric(cost,sampling=tuple(RES)).find_costs([src])
            for s in sites:
                s['geo']=float(cum[tuple(s['lumen_vox'])]); s['eu']=float(np.linalg.norm((np.array(s['lumen_vox'])-np.array(src))*RES))
        out['sites']=sites
        json.dump(out,open(fn,'w'))
    except Exception as e:
        json.dump(dict(nid=int(o['nid']),err=repr(e)),open(fn,'w'))
t=time.time()
with ThreadPoolExecutor(2) as ex:
    for k,_ in enumerate(ex.map(run,old)):
        if k%10==0: print(k,round(time.time()-t),flush=True)
print('DONE',flush=True)
