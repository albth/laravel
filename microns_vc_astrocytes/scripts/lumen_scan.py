import numpy as np,os,sys,time
from cloudvolume import CloudVolume
from scipy import ndimage as ndi
from concurrent.futures import ThreadPoolExecutor
em=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/em',mip=6,use_https=True,progress=False,fill_missing=True,parallel=1)
X0,X1=384,3328+216; Y0,Y1=512,2432; Z0,Z1=926,926+818
B=np.array([256,256,128]); M=np.array([8,8,4])
def proc(o):
    o=np.array(o); fn=f'blocks/{o[0]}_{o[1]}_{o[2]}.npz'
    if os.path.exists(fn): return 0
    lo=np.maximum(o-M,[X0,Y0,Z0]); hi=np.minimum(o+B+M,[X1,Y1,Z1])
    for attempt in range(4):
        try:
            e=np.asarray(em[lo[0]:hi[0],lo[1]:hi[1],lo[2]:hi[2]])[...,0].astype(np.float32); break
        except Exception as ex:
            time.sleep(2**attempt)
    else:
        print('FAIL',o,flush=True); return 0
    valid=e>0
    med=np.zeros(e.shape[2],np.float32); mad=np.ones(e.shape[2],np.float32)
    for k in range(e.shape[2]):
        v=e[:,:,k][valid[:,:,k]]
        if v.size>1000:
            m=np.median(v); med[k]=m; mad[k]=max(np.median(np.abs(v-m))*1.4826,1.0)
        else: med[k]=1e9
    zs=(e-med)/mad
    sm=ndi.uniform_filter(zs,size=(3,3,1))
    strong=sm>6.5; weak=sm>4.0
    lab,n=ndi.label(weak)
    keep=np.zeros(n+1,bool); keep[np.unique(lab[strong])]=True; keep[0]=False
    mask=keep[lab]
    # fill 2D holes (RBCs etc.)
    for k in range(mask.shape[2]): mask[:,:,k]=ndi.binary_fill_holes(mask[:,:,k])
    # crop margins
    s=o-lo; mask=mask[s[0]:s[0]+B[0],s[1]:s[1]+B[1],s[2]:s[2]+B[2]]
    vv=valid[s[0]:s[0]+B[0],s[1]:s[1]+B[1],s[2]:s[2]+B[2]]
    np.savez_compressed(fn,mask=np.packbits(mask,axis=None),shape=mask.shape,valid=np.packbits(vv,axis=None),origin=o)
    return int(mask.sum())
origins=[(x,y,z) for x in range(X0,X1,B[0]) for y in range(Y0,Y1,B[1]) for z in range(Z0,Z1,B[2])]
print('nblocks',len(origins),flush=True)
t=time.time(); done=0
with ThreadPoolExecutor(6) as ex:
    for r in ex.map(proc,origins):
        done+=1
        if done%25==0: print(done,len(origins),round(time.time()-t),flush=True)
print('DONE',flush=True)
