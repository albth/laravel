import numpy as np
R=np.load('vessel_layers.npy',allow_pickle=True).item()
use={'pial_art':'004 - pial arterioles.txt :: Arterioles','pen_art':'004_2 - Pen arterioles.txt :: Pen Arterioles',
'pial_ven':'005 - pial Venules.txt :: Venules','asc_ven':'005_2 - Asc venules.txt :: Asc Ven 1',
'acap1':'008_1 - arterial 1 order cap.txt :: 1. order','acap2':'008_1 - arterial 2 order cap.txt :: 2. order',
'acap3':'arterial 3 order cap.txt :: 3. order1','vcap1':'Venule 1 order cap.txt :: 1 order  ','vcap2':'Venule 2 order cap.txt :: 2 order ','vcap3':'Venule 3 order cap.txt :: 3 order'}
S=np.array([4,4,40.])/1000  # um per voxel
def segs(k):
    key=[kk for kk in R if kk.strip()==use[k].strip()][0]
    L=R[key][0]; return L[:,:3]*S, L[:,3:]*S
def sample(k,step=0.25):
    A,B=segs(k); d=np.linalg.norm(B-A,axis=1); pts=[]
    for a,b,l in zip(A,B,d):
        n=max(1,int(np.ceil(l/step))); t=np.linspace(0,1,n+1)[:,None]; pts.append(a+(b-a)*t)
    return np.vstack(pts), d.sum()
if __name__=='__main__':
    for k in use:
        A,B=segs(k); d=np.linalg.norm(B-A,axis=1)
        zs=np.unique(np.round(A[:,2]/0.04))
        print(f'{k:9s} nseg={len(A):6d} totlen={d.sum()/1000:8.1f} mm  medseg={np.median(d):.2f}um  dz_frac={np.mean(np.abs(B[:,2]-A[:,2])>0.001):.2f} uniqueZsections={len(zs)}')
