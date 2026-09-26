import numpy as np,pandas as pd,time
from cloudvolume import CloudVolume
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from lumenlib import RES,get_region
SEG=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/seg_m1300',use_https=True,progress=False,fill_missing=True)
SEG0=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/seg_m1300',mip=0,use_https=True,progress=False,fill_missing=True)
def root_at(x4,y4,z40):
    x,y,z=int(x4//2),int(y4//2),int(z40)
    v=np.asarray(SEG0[x-3:x+3,y-3:y+3,z:z+1]).ravel(); v=v[v>0]
    if len(v)==0: return 0
    u,c=np.unique(v,return_counts=True); return int(u[c.argmax()])
def analyze(r,R=45.0,thr=1.0,lod=2):
    rid=root_at(r.x,r.y,r.z)
    if rid==0: return None
    m=SEG.mesh.get(rid,lod=lod); m=list(m.values())[0] if isinstance(m,dict) else m
    V=m.vertices/1000.0; F=m.faces
    cen=np.array([r.x*0.004,r.y*0.004,r.z*0.04])
    extent=np.ptp(V,axis=0)
    keepv=np.linalg.norm(V-cen,axis=1)<R
    fk=keepv[F].all(1); F=F[fk]
    # face areas
    a=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    # lumen region + EDT
    lo=np.floor((cen-R-3)/RES).astype(int); hi=np.ceil((cen+R+3)/RES).astype(int)
    L,Vd=get_region(lo,hi)
    dist_out=ndi.distance_transform_edt(~L,sampling=RES)
    dist_in=ndi.distance_transform_edt(L,sampling=RES)   # lumen radius field
    def samp(field,P):
        ij=np.floor(P/RES).astype(int)-lo; ij=np.clip(ij,0,np.array(field.shape)-1); return field[ij[:,0],ij[:,1],ij[:,2]]
    vd=samp(dist_out,V)
    fc=V[F].mean(1); fd=vd[F].min(1)
    peri=fd<=thr
    tot=a.sum(); pa=a[peri].sum()
    # soma zone: faces within nucleus-radius+3um
    rn=(3*max(r.vol,1)/(4*np.pi))**(1/3)
    dcen=np.linalg.norm(fc-cen,axis=1)
    soma=dcen<rn+3
    soma_peri=a[peri&soma].sum()
    # contact sites: perivascular faces grouped on a 1um grid with 26-connectivity after 1-voxel dilation (merges fragments <~2um apart)
    pf=np.where(peri)[0]; patches=[]
    if len(pf):
        g=np.floor((fc[pf]-cen+R+2)/1.0).astype(int); gs=g.max(0)+2
        grid=np.zeros(gs,bool); grid[g[:,0],g[:,1],g[:,2]]=True
        lab,nc=ndi.label(ndi.binary_dilation(grid,structure=np.ones((3,3,3))),structure=np.ones((3,3,3)))
        fl=lab[g[:,0],g[:,1],g[:,2]]
        for k in range(1,nc+1):
            s=pf[fl==k]; ar=float(a[s].sum())
            if ar<2.0: continue
            pc=fc[s].mean(0)
            ij=np.floor(pc/RES).astype(int)-lo; w=np.ceil(8/RES).astype(int)
            sl=tuple(slice(max(0,ij[q]-w[q]),min(L.shape[q],ij[q]+w[q]+1)) for q in range(3))
            # caliber: lumen diameter at the lumen point nearest the patch
            sub=dist_in[sl]
            if sub.size and sub.max()>0:
                idx=np.argwhere(sub>0); P=(idx+np.array([x.start for x in sl])+lo+0.5)*RES
                j=np.argmin(np.linalg.norm(P-pc,axis=1)); nearest=P[j]
                # local max of inscribed radius within 3um of nearest lumen point
                near=np.linalg.norm(P-nearest,axis=1)<3
                cal=2*float(sub[tuple(idx[near].T)].max())
            else: cal=0.0; nearest=pc
            dmin=float(np.linalg.norm(V[F[s]].reshape(-1,3)-cen,axis=1).min())
            patches.append(dict(area=ar,dist_soma=float(np.linalg.norm(pc-cen)),dmin_soma=dmin,caliber=cal,pos=pc.tolist(),lumen_pt=list(map(float,nearest)),soma=bool(dmin<rn+3)))
    return dict(root=rid,n_faces=len(F),area_total=tot,area_peri=pa,frac_peri=pa/tot if tot else np.nan,soma_peri=soma_peri,
                n_patches=len(patches),patches=patches,extent=extent.tolist(),mesh_nverts=len(m.vertices))
if __name__=='__main__':
    df=pd.read_pickle('nuclei_L.pkl'); r=df[df.nid==290882].iloc[0]
    t=time.time(); out=analyze(r); print(time.time()-t)
    print({k:v for k,v in out.items() if k!='patches'})
    for p in sorted(out['patches'],key=lambda p:-p['area']): print(p)
