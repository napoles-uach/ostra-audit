# Explorer 1, 2026-11-23. Candidate fix, tested by monkeypatching PySCF 2.15.1: refuse NLC placeholders in the
# semilocal path of the XCFun backend with a clear error. Results (H2/cc-pVDZ, XCFun backend):
#   wb97m_v, scan_vv10, b3lyp+vv10 -> NotImplementedError (instead of SIGSEGV)
#   b3lyp + nlc=vv10 -1.16416231, scan + nlc=vv10 -1.15674014, b3lyp -1.17328112 (unchanged, equal to libxc)
from pyscf import gto, dft
from pyscf.dft import xcfun
def _guard(xc_code):
    if isinstance(xc_code, str):
        _, fn = xcfun.parse_xc(xc_code)
        if any(i >= 5000 for i, _ in fn):
            raise NotImplementedError(
                f"xc='{xc_code}' contains a VV10 functional; the XCFun backend has no implementation of its "
                "semilocal part. Use the libxc backend, or a semilocal xc with mf.nlc = 'vv10'.")
_xt, _ev = xcfun.xc_type, xcfun.eval_xc
def xc_type(xc_code): _guard(xc_code); return _xt(xc_code)
def eval_xc(xc_code, *a, **k): _guard(xc_code); return _ev(xc_code, *a, **k)
xcfun.xc_type, xcfun.eval_xc = xc_type, eval_xc
mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="cc-pvdz", verbose=0)
for xc, nlc in (("wb97m_v", "-"), ("scan_vv10", "-"), ("b3lyp+vv10", "-"), ("b3lyp", "vv10"), ("scan", "vv10"), ("b3lyp", "-")):
    mf = dft.RKS(mol); mf.xc = xc
    if nlc != "-": mf.nlc = nlc
    mf._numint.libxc = xcfun
    try: print(f"{xc:11s} nlc={nlc:5s} E={mf.kernel():.8f}")
    except Exception as e: print(f"{xc:11s} nlc={nlc:5s} {type(e).__name__}")
