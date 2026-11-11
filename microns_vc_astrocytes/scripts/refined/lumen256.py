import numpy as np
from cloudvolume import CloudVolume
from scipy import ndimage as ndi
EM=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/em',mip=[256,256,320],use_https=True,progress=False,fill_missing=True,parallel=4)
NUC=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/nuclei',mip=[256,256,320],use_https=True,progress=False,fill_missing=True,parallel=4)
RES=np.array([0.256,0.256,0.32])
def lumen_box(cen_um,half_um):
    lo=np.floor((np.array(cen_um)-half_um)/RES).astype(int); hi=np.ceil((np.array(cen_um)+half_um)/RES).astype(int)
    b=EM.bounds; lo=np.maximum(lo,b.minpt); hi=np.minimum(hi,b.maxpt)
    e=np.asarray(EM[lo[0]:hi[0],lo[1]:hi[1],lo[2]:hi[2]])[...,0].astype(np.float32)
    valid=e>0
    zs=np.zeros_like(e)
    for k in range(e.shape[2]):
        v=e[:,:,k][valid[:,:,k]]
        if v.size<1000: continue
        m=np.median(v); mad=max(np.median(np.abs(v-m))*1.4826,1.0); zs[:,:,k]=(e[:,:,k]-m)/mad
    sm=ndi.uniform_filter(zs,size=(5,5,1))
    strong=sm>6.0; weak=sm>3.5
    lab,n=ndi.label(weak); keep=np.zeros(n+1,bool); keep[np.unique(lab[strong])]=True; keep[0]=False
    L=keep[lab]
    nb=NUC.bounds; nlo=np.maximum(lo,nb.minpt); nhi=np.minimum(hi,nb.maxpt)
    nuc=np.zeros(L.shape,bool)
    if (nhi>nlo).all():
        a=np.asarray(NUC[nlo[0]:nhi[0],nlo[1]:nhi[1],nlo[2]:nhi[2]])[...,0]; s=nlo-lo
        nuc[s[0]:s[0]+a.shape[0],s[1]:s[1]+a.shape[1],s[2]:s[2]+a.shape[2]]=a>0
    lab,n=ndi.label(L)
    if n:
        tot=np.bincount(lab.ravel(),minlength=n+1); inn=np.bincount(lab[nuc].ravel(),minlength=n+1)
        bad=(inn>0.3*tot)|(tot<400); bad[0]=False; L&=~bad[lab]
    L&=~nuc
    for k in range(L.shape[2]): L[:,:,k]=ndi.binary_fill_holes(L[:,:,k])
    L=ndi.binary_opening(L,iterations=1)
    return L,lo,valid,e
