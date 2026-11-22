import sys, json, numpy as np, spinprobe as sp, lxc, collections
def pointwise(fid, low=True):
    fu, fp = lxc.Func(fid,1), lxc.Func(fid,2)
    if low: sp.lower_thresholds(fu); sp.lower_thresholds(fp)
    R,S,A,sigma,tau,lapl = sp.grid()
    xu={"rho":R[:,None],"sigma":sigma[:,None],"tau":tau[:,None],"lapl":lapl[:,None]}
    xp={"rho":np.stack([R/2,R/2],1),"sigma":np.stack([sigma/4]*3,1),"tau":np.stack([tau/2]*2,1),"lapl":np.stack([lapl/2]*2,1)}
    ou,op=sp.call(fu,xu),sp.call(fp,xp); fu.close(); fp.close()
    e=R*ou["zk"]; D={}
    with np.errstate(all="ignore"):
        D["zk"]=np.abs(op["zk"]-ou["zk"])/(np.abs(ou["zk"])+1e-300)
        if "vrho" in ou: D["vrho"]=np.abs(op["vrho"][:,0]-ou["vrho"][:,0])/(np.abs(ou["vrho"][:,0])+1e-300)
        if "vsigma" in ou:
            vs=ou["vsigma"][:,0]; sc=np.maximum(np.abs(vs)*sigma,np.abs(e))/np.maximum(sigma,1e-300)
            D["vsigma"]=np.abs(op["vsigma"].sum(1)/4-vs)/sc
        if "vtau" in ou:
            v=ou["vtau"][:,0]; sc=np.maximum(np.abs(v)*tau,np.abs(e))/np.maximum(tau,1e-300)
            D["vtau"]=np.abs(op["vtau"].mean(1)-v)/sc
    return R,S,A,D,ou,op
if __name__=="__main__":
    for a in sys.argv[1:]:
        fid=int(a); R,S,A,D,ou,op=pointwise(fid)
        bad=np.zeros(len(R),bool)
        for k,d in D.items(): bad|=np.nan_to_num(d,nan=0)>1e-6
        degen=(S==0)&(A==0); tauW=(A==0)&~degen
        cats=collections.Counter(("degen" if dg else "alpha0" if tw else "interior") for dg,tw in zip(degen[bad],tauW[bad]))
        print(fid, lxc.name(fid), dict(cats), "of", len(R))
        ib=np.where(bad&~degen&(A>0))[0]
        for i in ib[:4]:
            print("   rho=%.2e s=%g a=%g"%(R[i],S[i],A[i]), {k:"%.1e"%D[k][i] for k in D})
