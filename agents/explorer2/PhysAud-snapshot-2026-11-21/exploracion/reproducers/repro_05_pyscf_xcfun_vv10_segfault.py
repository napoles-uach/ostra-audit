#!/usr/bin/env python3
"""Minimal reproducer 05: PySCF 2.14.0 pyscf/dft/xcfun.py SIGSEGV (-11) on all 19 VV10/rVV10 functionals.
pyscf/dft/xcfun.py:368 inserts temporary dummy IDs 5000..5018 into XC_CODES for VV10_XC,
but eval_xc1 (xcfun.py:1037-1061) passes them unfiltered to C _itrf.XCFUN_eval_xc,
causing an immediate SIGSEGV crash (exit -11) with no Python exception.
"""
import subprocess, sys
from pyscf.dft import xcfun

vv10_names = sorted(k for k, v in xcfun.XC_CODES.items() if isinstance(v, int) and v >= 5000)
print(f"Total VV10/rVV10 dummy IDs >= 5000 in pyscf.dft.xcfun.XC_CODES: {len(vv10_names)}")
print("Sample:", vv10_names[:6])

print("\nTesting RKS(mol) with xcfun backend in isolated subprocesses:")
for xc in ["wb97m_v", "b97m_v", "scan_rvv10", "pbe+vv10"]:
    code = (
        "from pyscf import gto, dft\n"
        "mol = gto.M(atom='He 0 0 0', basis='sto-3g', verbose=0)\n"
        "mf = dft.RKS(mol); mf._numint.libxc = dft.xcfun; mf.xc = " + repr(xc) + "\n"
        "mf.kernel()\n"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    print(f"  xc={xc:12s} -> returncode={r.returncode} ({'SIGSEGV' if r.returncode == -11 else 'other'})")
