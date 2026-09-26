import numpy as np
from cloudvolume import CloudVolume
from PIL import Image
from scipy import ndimage as ndi
EM0=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/em',mip=0,use_https=True,progress=False,fill_missing=True)
SG0=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/seg_m1300',mip=0,use_https=True,progress=False,fill_missing=True)
def render(cnm,half_um,root,out,scale=1.0,z_off=0,outline=True,label=None):
    c=(np.array(cnm)/np.array([8,8,40])).astype(int); c[2]+=z_off; h=int(half_um*1000/8)
    e=np.asarray(EM0[c[0]-h:c[0]+h,c[1]-h:c[1]+h,c[2]:c[2]+1])[...,0,0].astype(float).T
    lo,hi=np.percentile(e[e>0],[0.5,99.7]) if (e>0).any() else (0,255)
    g=np.clip((e-lo)/(hi-lo+1e-6)*255,0,255).astype(np.uint8)
    rgb=np.stack([g]*3,-1)
    if outline and root:
        s=np.asarray(SG0[c[0]-h:c[0]+h,c[1]-h:c[1]+h,c[2]:c[2]+1])[...,0,0].T==root
        edge=s&~ndi.binary_erosion(s,iterations=2)
        rgb[edge]=[255,60,200]
    img=Image.fromarray(rgb)
    if scale!=1: img=img.resize((int(rgb.shape[1]*scale),int(rgb.shape[0]*scale)),Image.LANCZOS)
    img.save(out); return c
