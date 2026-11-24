# OSTRA-013 isolated: E_c[REVTPSSC] with z^2 (released XCFun) vs z^3 (fixed), on the same self-consistent
# revTPSS density (libxc, which agrees with fixed XCFun to 1e-11 Eh). Also r4SCAN exchange as a validation of the
# first-order (non-self-consistent) estimate against the self-consistent difference.
import sys, json, subprocess, os, numpy
sp, fam = sys.argv[1], sys.argv[2]
X = "/home/claude/src_xcfun"
dmf = f"dm/{sp}__{fam}__libxc.npy"
if not os.path.exists(dmf):
    from pyscf import gto, dft
    from mols import MOLS, ATOMS
    atom, spin = MOLS[sp] if sp in MOLS else (f"{sp} 0 0 0", ATOMS[sp])
    mol = gto.M(atom=atom, spin=spin, basis="def2-tzvp", verbose=0)
    mf = (dft.UKS if spin else dft.RKS)(mol); mf.grids.level = 5; mf.conv_tol = 1e-10
    mf.xc = {"revtpss": "mgga_x_revtpss,mgga_c_revtpss", "r4scan": "mgga_x_r4scan,mgga_c_r2scan"}[fam]
    mf.kernel(); numpy.save(dmf, numpy.asarray(mf.make_rdm1()))
xc = {"revtpss": ",revtpssc", "r4scan": "r4scanx,"}[fam]
res = {}
for tag, pre in (("bundled", ""), ("patched", f"{X}/patched/b/lib/libxcfun.so.2")):
    o = subprocess.run(["python3", "nsc.py", dmf, xc, tag, "def2-tzvp"], capture_output=True, text=True,
                       env={**os.environ, "LD_PRELOAD": pre, "OMP_NUM_THREADS": "1"})
    res[tag] = json.loads(o.stdout)["Exc"] if o.stdout else None
print(json.dumps({"sp": sp, "fam": fam, **res}))
