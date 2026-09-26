import numpy as np
from cloudvolume import CloudVolume
from PIL import Image
_em={}
def em(mip):
    if mip not in _em: _em[mip]=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/em',mip=mip,use_https=True,progress=False,fill_missing=True)
    return _em[mip]
def slice_nm(cnm,half_um,mip,out,scale=1,nz=1):
    v=em(mip); res=np.array(v.resolution); c=(np.array(cnm)/res).astype(int); h=int(half_um*1000/res[0])
    a=np.asarray(v[c[0]-h:c[0]+h,c[1]-h:c[1]+h,c[2]:c[2]+nz])[...,0].astype(float)
    sl=a.mean(axis=2).T
    lo,hi=np.percentile(sl[sl>0],[0.5,99.8]) if (sl>0).any() else (0,255)
    im=np.clip((sl-lo)/(hi-lo+1e-6)*255,0,255).astype(np.uint8)
    img=Image.fromarray(im)
    if scale!=1: img=img.resize((int(im.shape[1]*scale),int(im.shape[0]*scale)),Image.LANCZOS)
    img.save(out); return a
