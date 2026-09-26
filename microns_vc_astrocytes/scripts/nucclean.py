import numpy as np,glob,os,time
from cloudvolume import CloudVolume
from scipy import ndimage as ndi
from concurrent.futures import ThreadPoolExecutor
from lumenlib import load_block
NV=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/nuclei',mip=[512,512,640],use_https=True,progress=False,fill_missing=True)
bb=NV.bounds
stats=[]
def proc(fn):
    out='blocks2/'+os.path.basename(fn)
    if os.path.exists(out): return None
    o,m,v=load_block(fn)
    removed=0
    if m.any():
        lo=np.maximum(o,bb.minpt); hi=np.minimum(o+np.array(m.shape),bb.maxpt)
        nuc=np.zeros(m.shape,bool)
        if (hi>lo).all():
            for att in range(4):
                try: a=np.asarray(NV[lo[0]:hi[0],lo[1]:hi[1],lo[2]:hi[2]])[...,0]; break
                except Exception: time.sleep(2**att)
            s=lo-o; nuc[s[0]:s[0]+a.shape[0],s[1]:s[1]+a.shape[1],s[2]:s[2]+a.shape[2]]=a>0
        lab,n=ndi.label(m)
        if n:
            tot=np.bincount(lab.ravel(),minlength=n+1); innuc=np.bincount(lab[nuc].ravel(),minlength=n+1)
            bad=(innuc>0.3*tot); bad[0]=False
            before=m.sum(); m=m&~bad[lab]
            m=m&~nuc   # also strip voxels inside nuclei
            removed=int(before-m.sum())
    np.savez_compressed(out,mask=np.packbits(m,axis=None),shape=m.shape,valid=np.packbits(v,axis=None),origin=o)
    return removed
fns=sorted(glob.glob('blocks/*.npz')); t=time.time(); tot=0
with ThreadPoolExecutor(6) as ex:
    for k,r in enumerate(ex.map(proc,fns)):
        tot+=r or 0
        if k%100==0: print(k,tot,round(time.time()-t),flush=True)
print('DONE removed',tot,flush=True)
