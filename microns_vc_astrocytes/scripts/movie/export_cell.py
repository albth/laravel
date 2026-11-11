"""Export one astrocyte + local vessels for the three.js renderer.
Usage: python export_cell.py <nid> <root> <outname>"""
import numpy as np,pandas as pd,json,sys
from cloudvolume import CloudVolume
from scipy import ndimage as ndi
from skimage import measure
import trimesh,fast_simplification
from lumen256 import lumen_box,RES
nid,root,name=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
R=45.0
SEG=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/seg_m1300',use_https=True,progress=False)
nd=pd.read_csv('../mic/nucleus_detection_v0.csv',header=None,names=['id','valid','sv','root','x','y','z','vol']).set_index('id')
q=nd.loc[nid]; cen=np.array([q.x*0.004,q.y*0.004,q.z*0.04]); rn=(3*q.vol/(4*np.pi))**(1/3)
m=SEG.mesh.get(root,lod=2); m=list(m.values())[0] if isinstance(m,dict) else m
V=m.vertices/1000.0; F=m.faces
keep=np.linalg.norm(V-cen,axis=1)<R; F=F[keep[F].all(1)]
tm=trimesh.Trimesh(V,F,process=False); tm.remove_unreferenced_vertices()
# keep the largest connected piece containing the soma plus pieces near it
V2,F2=fast_simplification.simplify(tm.vertices.astype(np.float32),tm.faces.astype(np.int32),target_reduction=max(0,1-350000/len(tm.faces)))
tm=trimesh.Trimesh(V2,F2,process=True)
comps=tm.split(only_watertight=False)
keepc=[c for c in comps if np.linalg.norm(c.vertices-cen,axis=1).min()<R and len(c.faces)>200]
tm=trimesh.util.concatenate(keepc)
V,F=np.asarray(tm.vertices),np.asarray(tm.faces)
L,lo,valid,_=lumen_box(cen,R+4)
dout=ndi.distance_transform_edt(~L,sampling=RES)
ij=np.clip(np.floor(V/RES).astype(int)-lo,0,np.array(L.shape)-1); vd=dout[ij[:,0],ij[:,1],ij[:,2]]
fd=vd[F].mean(1); fc=V[F].mean(1); dc=np.linalg.norm(fc-cen,axis=1)
soma=dc<rn+3.0
peri=(fd>=0.15)&(fd<=1.0); lumi=fd<0.15
area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
cls=np.zeros(len(F),np.uint8)        # 0 process
cls[soma]=1                          # 1 soma
cls[lumi]=4                          # 4 luminal (merge suspect)
# cluster ALL perivascular faces; clusters touching the soma zone = soma contact + its sheath (3), others = endfeet (2)
pf=np.where(peri)[0]; clusters=[]
if len(pf):
    g=np.floor((fc[pf]-cen+R+2)/1.0).astype(int); gs=g.max(0)+2
    grid=np.zeros(gs,bool); grid[g[:,0],g[:,1],g[:,2]]=True
    gl,nc=ndi.label(ndi.binary_dilation(grid,structure=np.ones((3,3,3))),structure=np.ones((3,3,3)))
    fl=gl[g[:,0],g[:,1],g[:,2]]
    for k in range(1,nc+1):
        s_=pf[fl==k]; ar=float(area[s_].sum())
        if ar<3: continue
        touches=bool(soma[s_].any() and area[s_][soma[s_]].sum()>=2.0)
        cls[s_]=3 if touches else 2
        pc=fc[s_].mean(0); ps=fc[s_][np.argmin(np.linalg.norm(fc[s_]-pc,axis=1))]
        if touches:
            sz=s_[soma[s_]]; ps=fc[sz][np.argmin(np.linalg.norm(fc[sz]-fc[sz].mean(0),axis=1))]
        clusters.append(dict(area=ar,soma=touches,pos=(ps-cen).tolist(),dist=float(dc[s_].min())))
# vessel surface
Ls=ndi.binary_closing(dout<=0.35,iterations=1)  # display wall: detected lumen grown 0.35 um to the EM lumen edge
vv,ff,_,_=measure.marching_cubes(np.pad(Ls,1).astype(np.float32),0.5,spacing=tuple(RES))
vv=vv-RES+lo*RES
vm=trimesh.Trimesh(vv,ff,process=True)
trimesh.smoothing.filter_taubin(vm,iterations=15)
vm=vm.slice_plane(cen-R*0.98,[1,0,0]).slice_plane(cen+R*0.98,[-1,0,0]) if False else vm
inside=np.linalg.norm(vm.vertices-cen,axis=1)<R+2
fk=inside[vm.faces].all(1); vm=trimesh.Trimesh(vm.vertices,vm.faces[fk],process=True)
if len(vm.faces)>250000:
    a,b=fast_simplification.simplify(vm.vertices.astype(np.float32),vm.faces.astype(np.int32),target_reduction=1-250000/len(vm.faces)); vm=trimesh.Trimesh(a,b,process=True)
sites=[c for c in clusters if not c['soma']]
sc=[c for c in clusters if c['soma']]
soma_contact=max(sc,key=lambda c:c['area'])['pos'] if sc else None
# same-vessel test for endfeet: lumen geodesic from soma contact
if sc and sites:
    from skimage.graph import MCP_Geometric
    lab_,_=ndi.label(L)
    def nearest_lumen(p):
        ij=np.floor((np.array(p)+cen)/RES).astype(int)-lo; w=np.ceil(3/RES).astype(int)
        sl=tuple(slice(max(0,ij[t]-w[t]),min(L.shape[t],ij[t]+w[t]+1)) for t in range(3)); idx=np.argwhere(L[sl])
        if len(idx)==0: return None
        idx=idx+np.array([x.start for x in sl]); return tuple(idx[np.argmin(np.linalg.norm((idx+lo+0.5)*RES-(np.array(p)+cen),axis=1))])
    src=nearest_lumen(soma_contact)
    if src is not None:
        cum,_=MCP_Geometric(np.where(L,1.0,np.inf),sampling=tuple(RES)).find_costs([src])
        for s_ in sites:
            t_=nearest_lumen(s_['pos'])
            if t_ is None: continue
            geo=float(cum[t_]); eu=float(np.linalg.norm((np.array(t_)-np.array(src))*RES))
            s_['geo']=geo if np.isfinite(geo) else None; s_['eu']=eu; s_['same_vessel']=bool(np.isfinite(geo) and geo<1.5*eu)
# nearest lumen to soma surface
somav=np.where(np.linalg.norm(V-cen,axis=1)<rn+3)[0]
near_lumen=float(vd[somav].min()) if len(somav) else None
meta=dict(nid=nid,root=root,center_um=cen.tolist(),center_vox=[int(q.x),int(q.y),int(q.z)],rn=float(rn),sites=sorted(sites,key=lambda s:-s['area']),
          soma_contact=soma_contact,soma_min_dist=near_lumen,area_by_class={int(k):float(area[cls==k].sum()) for k in range(5)},
          nverts=len(V),nfaces=len(F),vessel_nverts=len(vm.vertices),vessel_nfaces=len(vm.faces))
out=f'web/{name}'
(V-cen).astype('<f4').tofile(out+'_astro_v.bin')
order=np.argsort(cls,kind='stable'); F=F[order]; cls=cls[order]
F.astype('<u4').tofile(out+'_astro_f.bin')
meta['class_ranges']={int(k):[int(np.searchsorted(cls,k)),int(np.searchsorted(cls,k,side='right'))] for k in range(5)}
(np.asarray(vm.vertices)-cen).astype('<f4').tofile(out+'_vessel_v.bin'); np.asarray(vm.faces).astype('<u4').tofile(out+'_vessel_f.bin')
json.dump(meta,open(out+'_meta.json','w'),indent=1)
print(json.dumps({k:v for k,v in meta.items() if k!='sites'},indent=0)[:1500]); print('sites',[(round(s['area']),round(s['dist'])) for s in meta['sites']])
