import numpy as np, warnings; warnings.filterwarnings("ignore")
from pyscf.dft import libxc, xcfun
rng=np.random.default_rng(1); n=400
rho=10**rng.uniform(-3,1.5,n); s=rng.uniform(0.02,4,n); a=rng.uniform(0.02,4,n)
kf=(3*np.pi**2*rho)**(1/3); sig=(2*kf*rho*s)**2; tu=0.3*(3*np.pi**2)**(2/3)*rho**(5/3); tw=sig/(8*rho); tau=tw+a*tu
R=np.vstack([rho,np.sqrt(sig/3)*np.ones((3,n)),np.zeros(n),tau]); z=tw/tau
for xf,lx in (("REVTPSSX","MGGA_X_REVTPSS"),("REVTPSSC","MGGA_C_REVTPSS")):
    e1=xcfun.eval_xc(xf,R,spin=0,deriv=0)[0]; e2=libxc.eval_xc(lx,R,spin=0,deriv=0)[0]; d=np.abs(e1-e2)/np.abs(e2)
    print(xf,"median %.1e  90%% %.1e  max %.1e"%(np.median(d),np.percentile(d,90),d.max()))
    for lab,v in (("rho",rho),("s",s),("alpha",a),("z",z)): print("   corr(log d, %s) = %+.2f"%(lab,np.corrcoef(np.log10(d+1e-17),np.log10(v))[0,1]))
    i=np.argsort(d)[-3:]; print("   worst:",[(round(rho[k],4),round(s[k],2),round(a[k],2),"%.1e"%d[k]) for k in i])
