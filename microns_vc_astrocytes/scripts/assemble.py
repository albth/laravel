import numpy as np,glob
from scipy import ndimage as ndi
from lumenlib import load_block,RES
surf=[];nl=0;nv=0;minsz=40
allc=[]
for fn in sorted(glob.glob('blocks/*.npz')):
    o,m,v=load_block(fn)
    nv+=v.sum()
    if not m.any(): continue
    lab,n=ndi.label(m); sz=np.bincount(lab.ravel()); keep=sz>=minsz; keep[0]=False; m=keep[lab]
    nl+=m.sum()
    er=ndi.binary_erosion(m,structure=ndi.generate_binary_structure(3,1),border_value=0)
    s=np.argwhere(m&~er)+o
    surf.append(s.astype(np.int32))
S=np.vstack(surf)
print('lumen voxels',nl,'valid voxels',nv,'lumen vol frac',nl/nv,'surface pts',len(S))
np.save('lumen_surface_vox6.npy',S)
