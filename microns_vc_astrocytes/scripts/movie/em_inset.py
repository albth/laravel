"""EM inset: 32 nm section through a point, astrocyte segment painted with the movie's colours.
usage: python em_inset.py <meta.json> <which: contact|site0|site1..> <out.png> [half_um]"""
import numpy as np,json,sys
from cloudvolume import CloudVolume
from scipy import ndimage as ndi
from PIL import Image,ImageDraw
meta=json.load(open(sys.argv[1])); which=sys.argv[2]; out=sys.argv[3]; half=float(sys.argv[4]) if len(sys.argv)>4 else 7
cen=np.array(meta['center_um'])
p=(np.zeros(3) if which=='soma' else np.array(meta['soma_contact'] if which=='contact' else meta['sites'][int(which[4:])]['pos']))+cen
EM=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/em',mip=[32,32,40],use_https=True,progress=False,fill_missing=True)
SG=CloudVolume('precomputed://gs://iarpa_microns/minnie/minnie65/seg_m1300',mip=[32,32,40],use_https=True,progress=False,fill_missing=True)
res=np.array([0.032,0.032,0.04]); c=(p/res).astype(int); h=int(half/0.032)
e=np.asarray(EM[c[0]-h:c[0]+h,c[1]-h:c[1]+h,c[2]:c[2]+1])[...,0,0].astype(np.float32).T
s=np.asarray(SG[c[0]-h:c[0]+h,c[1]-h:c[1]+h,c[2]:c[2]+1])[...,0,0].T
v=e[e>0]
from lumen256 import lumen_box,RES
L,lo3,_,_=lumen_box(p,max(half+2,30))
k=int(np.floor(p[2]/RES[2]))-lo3[2]
sl=L[:,:,k].astype(np.float32).T
# map 32nm grid -> 256nm grid
xs=(np.arange(2*h)+c[0]-h)*0.032; ys=(np.arange(2*h)+c[1]-h)*0.032
ix=np.clip((xs/RES[0]-lo3[0]),0,sl.shape[1]-1); iy=np.clip((ys/RES[1]-lo3[1]),0,sl.shape[0]-1)
up=ndi.map_coordinates(ndi.gaussian_filter(sl,0.7),np.meshgrid(iy,ix,indexing='ij'),order=1)
lum0=up>0.5
lum=ndi.distance_transform_edt(~lum0)*0.032<=0.35
astro=s==meta['root']
d=ndi.distance_transform_edt(~lum0)*0.032
yy,xx=np.mgrid[0:2*h,0:2*h]; pts=np.stack([xx+c[0]-h,yy+c[1]-h],-1)*0.032
dc=np.sqrt(((pts-cen[:2])**2).sum(-1)+(p[2]-cen[2])**2)
soma=astro&(dc<meta['rn']+3)
peri=astro&(d<=1.0)&(d>=0.15)
lo,hi=np.percentile(v,[0.5,99.7]); g=np.clip((e-lo)/(hi-lo)*255,0,255)
rgb=np.stack([g]*3,-1)
def paint(mask,col,a): rgb[mask]=rgb[mask]*(1-a)+np.array(col)*a
paint(astro&~soma&~peri,(127,211,230),0.45); paint(soma&~peri,(138,92,240),0.55)
paint(peri&soma,(255,138,61),0.8); paint(peri&~soma,(245,197,66),0.8)
edge=lum&~ndi.binary_erosion(lum,iterations=3); rgb[edge]=(195,58,69)
im=Image.fromarray(rgb.astype(np.uint8)).resize((520,520),Image.LANCZOS)
dr=ImageDraw.Draw(im); px=520/(2*half); dr.rectangle([20,494,20+5*px,500],fill=(255,255,255))
im.save(out); print('saved',out,'astro px',astro.sum(),'peri px',peri.sum(),'lumen px',lum.sum())
