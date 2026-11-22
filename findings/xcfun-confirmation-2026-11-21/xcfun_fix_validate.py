"""XCFun (PySCF backend) vs libxc on 400 random meta-GGA points; run under LD_PRELOAD=<libxcfun.so> to test a build."""
import numpy as np, warnings; warnings.filterwarnings("ignore")
from pyscf.dft import libxc, xcfun
rng=np.random.default_rng(1); n=400
rho=10**rng.uniform(-3,1.5,n); s=rng.uniform(0.02,4,n); a=rng.uniform(0.02,4,n)
kf=(3*np.pi**2*rho)**(1/3); sig=(2*kf*rho*s)**2; tu=0.3*(3*np.pi**2)**(2/3)*rho**(5/3); tau=sig/(8*rho)+a*tu
R=np.vstack([rho,np.sqrt(sig/3)*np.ones((3,n)),np.zeros(n),tau])
for xf,lx in (("REVTPSSX","MGGA_X_REVTPSS"),("REVTPSSC","MGGA_C_REVTPSS"),("R4SCANX","MGGA_X_R4SCAN"),
              ("TPSSX","MGGA_X_TPSS"),("TPSSC","MGGA_C_TPSS"),("R2SCANX","MGGA_X_R2SCAN"),("SCANX","MGGA_X_SCAN")):
    try:
        e1,v1=xcfun.eval_xc(xf,R,spin=0,deriv=1)[:2]; e2,v2=libxc.eval_xc(lx,R,spin=0,deriv=1)[:2]
        de=np.max(np.abs(e1-e2)/np.abs(e2)); dv=max(np.max(np.abs(v1[i]-v2[i])/(np.abs(v2[i])+1e-12)) for i in (0,1,3))
        print("%-9s type=%-4s  eps rel %.1e   vrho/vsigma/vtau rel %.1e"%(xf,xcfun.xc_type(xf),de,dv))
    except Exception as ex: print("%-9s type=%-4s  ERROR %s"%(xf,xcfun.xc_type(xf),str(ex)[:60]))
for ra,rb in ((0.7,0.3),(0.9,0.1),(1.0,0.0)):
    pol=xcfun.eval_xc("TFK",np.array([[ra],[rb]]),spin=1,deriv=0)[0][0]*(ra+rb)
    u=lambda r: xcfun.eval_xc("TFK",np.array([r]),spin=0,deriv=0)[0][0]*r if r>0 else 0.0
    ex=(u(2*ra)+u(2*rb))/2; print("TFK zeta=%.1f  polarized vs exact spin scaling: rel %.1e"%(ra-rb,pol/ex-1))
