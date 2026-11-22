import numpy as np
from pyscf import gto, symm
from scipy.spatial.transform import Rotation
L="SPDFGHIKLM"
def mk(geom, shells, cart=False):
    atoms=[[f"X{i}",c] for i,c in enumerate(geom)]
    basis={f"X{i}":gto.basis.parse("".join(f"X {L[l]}\n {e:.17g} 1.0\n" for l,e in sh)) for i,sh in enumerate(shells)}
    m=gto.Mole(atom=atoms,basis=basis,unit="bohr",verbose=0,cart=cart); m.build(dump_input=False,parse_arg=False); return m
# rotation matrix acting on real spherical harmonics of the molecule: build from AO overlap trick:
# D_l = S(rot ao, orig ao) for a single-centre shell -> exact numerically (overlap ints are accurate)
def Dmat(l, Rm, a=1.0):
    # overlap between shell at origin and same shell evaluated with rotated coordinates:
    # use PySCF's real spherical harmonics on random points to fit D by least squares (exact linear map)
    pts=np.random.default_rng(0).normal(size=(200,3))
    from pyscf.symm.sph import real_sph_vec
    Y=np.array(real_sph_vec(pts/np.linalg.norm(pts,axis=1)[:,None], l, True)[l]).T  # (np, 2l+1)
    Yr=np.array(real_sph_vec((pts@Rm.T)/np.linalg.norm(pts,axis=1)[:,None], l, True)[l]).T
    D,*_=np.linalg.lstsq(Y,Yr,rcond=None); return D  # Y(R r) = Y(r) D
for lmax,R,a,b in ((6,4.75,0.3,3.0),(7,5.5,0.3,3.0),(4,3.75,0.3,3.0)):
    g=np.array([[0,0,0],[0.3*R,-0.2*R,0.93*R]]); sh=[[(lmax,a)],[(lmax,b)]]
    m=mk(g,sh); n=m.nao; e=m.intor("int2e",aosym="s1").reshape(n,n,n,n)
    Rm=Rotation.random(random_state=11).as_matrix()
    e2=mk(g@Rm.T,sh).intor("int2e",aosym="s1").reshape(n,n,n,n)
    d=2*lmax+1; D=Dmat(lmax,Rm); 
    print("D orthogonal err", np.abs(D.T@D-np.eye(d)).max())
    B=np.zeros((n,n)); B[:d,:d]=D; B[d:,d:]=D
    for cand in (B, B.T):
        er=np.einsum('ijkl,ia,jb,kc,ld->abcd',e,cand,cand,cand,cand,optimize=True)
        diff=np.abs(er-e2)
        if diff.max()<0.5*np.abs(e2).max(): break
    i=np.unravel_index(diff.argmax(),diff.shape)
    big=np.abs(e2)>1e-3*np.abs(e2).max()
    relbig=(diff[big]/np.abs(e2[big])).max()
    print(f"l={lmax} R={R}: max |ERI|={np.abs(e2).max():.3e}  max abs dev={diff.max():.3e} at {i} (value {e2[i]:.3e})  max rel dev over elements >1e-3*max: {relbig:.2e}")
    # cartesian path for comparison
    mc=mk(g,sh,cart=True); c2s=mc.cart2sph_coeff(); nc=mc.nao
    ec=mc.intor("int2e",aosym="s1").reshape(nc,nc,nc,nc)
    es=np.einsum('ijkl,ia,jb,kc,ld->abcd',ec,c2s,c2s,c2s,c2s,optimize=True)
    print(f"      cart->sph vs sph: max abs dev {np.abs(es-e).max():.3e}")
