# Explorer 1, 2026-11-23. What a filter-only fix would do: skip the NLC placeholder ids (>= 5000)
# in the XCFun backend and compare total energies with the libxc backend (H2, cc-pVDZ).
from pyscf import gto, dft
from pyscf.dft import xcfun
_parse = xcfun.parse_xc
def parse_sl(desc):
    hyb, fn = _parse(desc)
    return hyb, tuple((i, c) for i, c in fn if i < 5000)
xcfun.parse_xc = parse_sl
mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="cc-pvdz", verbose=0)
for xc in ("scan_vv10", "lc_vv10", "b97m_v", "wb97m_v"):
    e = []
    for be in ("libxc", "xcfun"):
        mf = dft.RKS(mol); mf.xc = xc
        if be == "xcfun": mf._numint.libxc = xcfun
        e.append(mf.kernel())
    print(f"{xc:10s} libxc {e[0]:.8f}  xcfun(filtered) {e[1]:.8f}  semilocal ids left: {parse_sl(xc)[1]}")
