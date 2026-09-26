import numpy as np,glob,os
from scipy import ndimage as ndi
RES=np.array([0.512,0.512,0.64])  # um per mip6 voxel
B=np.array([256,256,128])
def load_block(fn):
    d=np.load(fn); sh=tuple(d['shape']); n=np.prod(sh)
    m=np.unpackbits(d['mask'])[:n].reshape(sh).astype(bool)
    v=np.unpackbits(d['valid'])[:n].reshape(sh).astype(bool)
    return d['origin'],m,v
def get_region(lo,hi,bdir='blocks'):
    """lo,hi in mip6 voxel coords -> (mask,valid) arrays"""
    lo=np.array(lo);hi=np.array(hi)
    X0,Y0,Z0=384,512,926
    M=np.zeros(hi-lo,bool); V=np.zeros(hi-lo,bool)
    b0=(lo-[X0,Y0,Z0])//B; b1=(hi-1-np.array([X0,Y0,Z0]))//B
    for bx in range(b0[0],b1[0]+1):
        for by in range(b0[1],b1[1]+1):
            for bz in range(b0[2],b1[2]+1):
                o=np.array([X0,Y0,Z0])+B*[bx,by,bz]
                fn=f'{bdir}/{o[0]}_{o[1]}_{o[2]}.npz'
                if not os.path.exists(fn): continue
                _,m,v=load_block(fn)
                s0=np.maximum(lo,o); s1=np.minimum(hi,o+np.array(m.shape))
                if (s1<=s0).any(): continue
                M[tuple(slice(a-l,b-l) for a,b,l in zip(s0,s1,lo))]=m[tuple(slice(a-oo,b-oo) for a,b,oo in zip(s0,s1,o))]
                V[tuple(slice(a-l,b-l) for a,b,l in zip(s0,s1,lo))]=v[tuple(slice(a-oo,b-oo) for a,b,oo in zip(s0,s1,o))]
    return M,V
def nm_to_vox6(p_nm): return np.asarray(p_nm)/ (RES*1000)
_G={}
def lookup(P_um):
    """P in um -> (in_lumen, valid) using global packed arrays"""
    if not _G:
        _G['L']=np.load('global_lumen_bits.npy',mmap_mode='r'); _G['V']=np.load('global_valid_bits.npy',mmap_mode='r')
    vox=np.floor(np.asarray(P_um)/RES).astype(int)-np.array([384,512,926])
    Lb=_G['L']; sh=np.array([Lb.shape[0],Lb.shape[1],Lb.shape[2]*8])
    ok=np.all((vox>=0)&(vox<sh),axis=1); inl=np.zeros(len(vox),bool); val=np.zeros(len(vox),bool)
    v=vox[ok]; bit=7-(v[:,2]%8)
    inl[ok]=(Lb[v[:,0],v[:,1],v[:,2]//8]>>bit)&1; val[ok]=(_G['V'][v[:,0],v[:,1],v[:,2]//8]>>bit)&1
    return inl,val
