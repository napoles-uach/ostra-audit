TITLE: [BUG REPORT] Bundled XCFun (a89b783): wrong results for REVTPSS, R4SCANX, TFK, TW and SPBEC with the xcfun backend; fix can go into libxcfun.patch

**Describe the bug**
PySCF builds XCFun from dftlibs/xcfun at `a89b783` and applies `pyscf/lib/libxcfun.patch`. Several functionals exposed in `pyscf.dft.xcfun.XC_CODES` give wrong results at that commit; details, references to the original papers and a tested patch are in dftlibs/xcfun#170. That repository has had no commits since November 2022, so a fix may be easier on the PySCF side.

Affected names (all in `XC_CODES`):

| Name | Problem (see dftlibs/xcfun#170) | Effect with `mf._numint.libxc = dft.xcfun` |
|---|---|---|
| `REVTPSSX` / `REVTPSS` | missing `XC_KINETIC` flag and two formula errors | NaN; SCF does not converge |
| `REVTPSSC` | z² instead of z³ | silent, up to 1.35 kcal/mol in atomization energies |
| `R4SCANX` / `R4SCAN` | pre-publication damping parameters (0.232 instead of 0.802 / 0.178) | silent, up to 2.0 kcal/mol in atomization energies |
| `TFK` | ignores spin polarization | -8 % at ζ = 0.4 (open-shell only) |
| `TW` | squares already-squared gradients | H atom T_W = 0.0075 Eh instead of 0.5 |
| `SPBEC` | β and γ swapped | correlation energy 64 % too negative |

The libxc backend (default) is not affected.

**To Reproduce**
```python
from pyscf import gto, dft
mol = gto.M(atom="O 0 0 0.117; H 0 0.757 -0.469; H 0 -0.757 -0.469", basis="def2-svp", verbose=0)
for xc in ("r4scanx,r4scanc", "tpssx,revtpssc"):
    e = []
    for backend in (dft.libxc, dft.xcfun):
        mf = dft.RKS(mol); mf.xc = xc; mf._numint.libxc = backend
        e.append(mf.kernel())
    print(xc, "libxc - xcfun = %.2f mEh" % ((e[0] - e[1]) * 1000))
```

**Expected behavior**
Both backends should agree (to the small differences from constant conventions, < 1e-3 mEh after the fix).

**Version and hardware**
- PySCF version: 2.15.1 (same XCFun commit and libxcfun.patch on master 6b288fd)
- PySCF install type: PyPI/pip
- Operating system: Linux x86_64 and macOS arm64

**Additional context**
The patch attached to dftlibs/xcfun#170 applies cleanly on top of `libxcfun.patch` (a89b783 + libxcfun.patch + fix). Built that way, XCFun's `testall` passes (the patch also updates two self-test references that had been generated with the incorrect code), unaffected functionals are bit-identical, and the fixed functionals agree with libxc 7.1.2. Two options:
1. Append the fix to `pyscf/lib/libxcfun.patch` (or apply it as a second patch).
2. If you prefer not to patch upstream code, remove or warn for the affected names in `XC_CODES` until XCFun releases a fix.

Happy to open a PR for either option.

---
Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of electronic-structure libraries. Analysis assisted by AI models (Claude, Anthropic; Gemini, Google); all results were checked and the report reviewed by a human before submission.
