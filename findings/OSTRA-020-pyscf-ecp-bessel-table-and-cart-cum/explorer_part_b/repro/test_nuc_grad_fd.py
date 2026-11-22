import ctypes, os, subprocess, sys
import numpy as np

patched_lib_path = sys.argv[1] if len(sys.argv) > 1 else None

def run_case(label, atom_fmt, r0, basis, ecp, do_fd=True, hs=(1e-3, 1e-4)):
    from pyscf import gto, scf
    import pyscf.gto.moleintor as mi
    import pyscf.gto.ecp as gecp
    if patched_lib_path:
        lib = ctypes.CDLL(patched_lib_path)
        mi.libcgto = lib
        gecp.libcgto = lib

    mol = gto.M(atom=atom_fmt.format(r=r0), basis=basis, ecp=ecp, unit="Bohr", verbose=0)
    max_lb = max(mol.bas_angular(i) for i in range(mol.nbas))
    ecp_ch = sorted(set(mol._ecpbas[:, 1].tolist()))
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.kernel()
    assert mf.converged, f"SCF did not converge for {label}"
    g = mf.nuc_grad_method().kernel()
    gz = g[1, 2]
    print(f"[{label}] max_lb={max_lb}, ecp_ch={ecp_ch}, nao={mol.nao}")
    print(f"  RHF E = {mf.e_tot:.12f} Eh")
    print(f"  Analytical grad_z(atom 1) = {gz:+.12e} Eh/Bohr")
    print(f"  Translational sum |g0+g1| = {np.linalg.norm(g[0]+g[1]):.3e}")
    if do_fd:
        for h in hs:
            es = []
            for dr in (-2*h, -h, +h, +2*h):
                m = gto.M(atom=atom_fmt.format(r=r0+dr), basis=basis, ecp=ecp, unit="Bohr", verbose=0)
                s = scf.RHF(m)
                s.conv_tol = 1e-12
                s.conv_tol_grad = 1e-10
                es.append(s.kernel())
            fd4 = (es[0] - 8*es[1] + 8*es[2] - es[3]) / (12*h)
            diff = abs(gz - fd4)
            rel = diff / abs(fd4)
            print(f"  FD4 (h={h:.0e} Bohr): {fd4:+.12e}  |anal-FD4|={diff:.3e}  rel={rel:.3e}")

if __name__ == "__main__":
    # Case A: compact diatomic with l=6 (i-function) on ECP atom + local & nonlocal ECP (fast FD4!)
    run_case(
        "NaH_l6_custom",
        "Na 0 0 0; H 0 0 {r:.12f}",
        3.5,
        basis={"Na": [[0, (0.5, 1.0)], [1, (0.3, 1.0)], [6, (1.2, 1.0)]], "H": "sto-3g"},
        ecp={"Na": [10, [[-1, [[], [[1.0, 2.0]]]], [0, [[], [], [[1.5, 3.0]]]], [1, [[], [], [[1.5, 3.0]]]]]]},
        do_fd=True, hs=(1e-2, 1e-3, 1e-4),
    )
    # Case B: REAL built-in PySCF molecule: AgH with built-in cc-pv5z-pp basis (l=6!) + built-in def2-svp ECP (local channel -1!)
    run_case(
        "AgH_builtin_ccpv5zpp_def2ecp",
        "Ag 0 0 0; H 0 0 {r:.12f}",
        3.06,
        basis={"Ag": "cc-pv5z-pp", "H": "sto-3g"},
        ecp={"Ag": "def2-svp"},
        do_fd=True, hs=(1e-3,),
    )
