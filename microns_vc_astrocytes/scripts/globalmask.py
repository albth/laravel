import numpy as np,glob
from lumenlib import load_block,B
X0=np.array([384,512,926]); X1=np.array([3544,2432,1744])
sh=X1-X0; sh=((sh+B-1)//B)*B
L=np.zeros((sh[0],sh[1],sh[2]//8),np.uint8); V=np.zeros_like(L)
for fn in glob.glob('blocks/*.npz'):
    o,m,v=load_block(fn); s=o-X0
    mm=np.zeros(tuple(B),bool); vv=np.zeros(tuple(B),bool); mm[:m.shape[0],:m.shape[1],:m.shape[2]]=m; vv[:v.shape[0],:v.shape[1],:v.shape[2]]=v
    L[s[0]:s[0]+B[0],s[1]:s[1]+B[1],s[2]//8:(s[2]+B[2])//8]=np.packbits(mm,axis=2)
    V[s[0]:s[0]+B[0],s[1]:s[1]+B[1],s[2]//8:(s[2]+B[2])//8]=np.packbits(vv,axis=2)
np.save('global_lumen_bits.npy',L); np.save('global_valid_bits.npy',V); print(L.shape)
