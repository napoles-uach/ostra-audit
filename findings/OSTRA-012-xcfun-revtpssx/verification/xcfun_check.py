import numpy as np, warnings
warnings.filterwarnings("ignore")
from pyscf.dft import libxc, xcfun
rng=np.random.default_rng(0); n=200
rho=10**rng.uniform(-3,1,n); s=rng.uniform(0.05,3,n); a=rng.uniform(0.05,3,n)
kf=(3*np.pi**2*rho)**(1/3); sig=(2*kf*rho*s)**2; tw=sig/(8*rho); tu=0.3*(3*np.pi**2)**(2/3)*rho**(5/3)
tau=tw+a*tu
R=np.vstack([rho,np.sqrt(sig/3)*np.ones((3,n)),np.zeros(n),tau])  # PySCF MGGA layout: rho, dx,dy,dz, lapl, tau
def e(mod,code): return mod.eval_xc(code,R,spin=0,deriv=0)[0]
for lx,xf in (("MGGA_X_REVTPSS","REVTPSSX"),("MGGA_C_REVTPSS","REVTPSSC"),("MGGA_X_R4SCAN","R4SCANX"),("MGGA_X_TPSS","TPSSX"),("MGGA_C_TPSS","TPSSC"),("MGGA_X_R2SCAN","R2SCANX")):
    try:
        a1=e(libxc,lx); a2=e(xcfun,xf)
        d=np.abs(a2-a1)/np.abs(a1); print("%-9s vs libxc %-15s max rel diff %.1e  nan_xcfun=%d"%(xf,lx,np.nanmax(d),np.isnan(a2).sum()))
    except Exception as ex: print(xf,"error",type(ex).__name__,str(ex)[:80])
# TFK spin scaling: T[ra,rb] must equal (T[2ra]+T[2rb])/2
ra=np.array([0.7]); rb=np.array([0.3])
pol=xcfun.eval_xc("TFK",(np.vstack([ra,np.zeros((3,1))]),np.vstack([rb,np.zeros((3,1))])),spin=1,deriv=0)[0]*(ra+rb)
u=lambda r: xcfun.eval_xc("TFK",np.vstack([r,np.zeros((3,1))]),spin=0,deriv=0)[0]*r
print("TFK zeta=0.4: polarized %.6f  exact spin-scaled %.6f  rel %.2f%%"%(pol[0],((u(2*ra)+u(2*rb))/2)[0],100*(pol[0]/((u(2*ra)+u(2*rb))/2)[0]-1)))
