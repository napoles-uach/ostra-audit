#!/usr/bin/env python3
"""Step 5 observable impact measurements in PySCF (def2-TZVP) for candidate findings.

Measures real-world molecular/atomic impact (energies, atomization energies, Kohn-Sham
matrix consistency, PES derivative noise, grid occupancy, and FHC API propagation) for:
  1. xc_E1_scaled_jet (src/expint_e1.c:140) in HSE06 & LC-wPBE (devel vs 7.0.0 / exact)
  2. s_scaling_2 (maple/scaling.mpl:16-19) at s=1 & s=15 in HSE06 & LC-wPBE (libxc vs smooth CP2K-form / unclipped)
  3. pc07_fab (maple/mgga_exc/mgga_k_pc07.mpl:27) at z=a/40 in r2SCAN-L & SCAN-L
  4. xc_func_set_fhc_enforcement (src/functionals.c:511-518) across allXC_FLAGS_ENFORCE_FHC meta-GGAs
  5. lda_c_pk09 (maple/lda_exc/lda_c_pk09.mpl:25,44) pole/NaN shell at rs ~ 22.74 on open-shell Li & N
"""
import ctypes
import math
import os
import subprocess
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
if "LIBXC_PATH" not in os.environ:
    os.environ["LIBXC_PATH"] = str(ROOT / "build-libxc" / "libxc.dylib")

import lxc  # noqa: E402
from pyscf import gto, dft  # noqa: E402

lxc.L.xc_func_set_fhc_enforcement.argtypes = [ctypes.c_void_p, ctypes.c_int]
lxc.L.xc_func_set_sigma_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
lxc.L.xc_func_set_tau_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]

PI = math.pi
X2S = 1.0 / (2.0 * (6.0 * PI**2) ** (1.0 / 3.0))
XC_FLAGS_NEEDS_TAU = 1 << 16
XC_FLAGS_ENFORCE_FHC = 1 << 17


def check_fhc_contract():
    import json
    sg = json.loads((ROOT / "scripts" / "source_graph.json").read_text())
    mix_fids = {int(fid) for fid, info in sg["funcs"].items() if info.get("aux_ids")}
    print("=== 1. API CONTRACT: xc_func_set_fhc_enforcement across all ENFORCE_FHC meta-GGAs ===")
    pt = {
        "rho": np.array([[0.12]]),
        "sigma": np.array([[0.10]]),  # tau_W = 0.10 / 0.96 = 0.10417 > tau=0.02 => alpha < 0
        "lapl": np.array([[0.0]]),
        "tau": np.array([[0.02]]),
    }
    single_ok, single_inv, mix_ok, mix_broken = [], [], [], []
    for fid in sorted(set(lxc.all_ids())):
        f1 = lxc.Func(fid, 1)
        if f1.family == "MGGA" and (f1.flags & XC_FLAGS_ENFORCE_FHC) and (f1.flags & XC_FLAGS_NEEDS_TAU) and (f1.flags & lxc.FLAG_EXC):
            o1 = f1.exc_vxc(pt)
            f2 = lxc.Func(fid, 1)
            lxc.L.xc_func_set_fhc_enforcement(f2.p, 0)
            o2 = f2.exc_vxc(pt)
            f2.close()
            delta = abs(float(o1["zk"][0]) - float(o2["zk"][0])) + abs(float(o1["vtau"][0, 0]) - float(o2["vtau"][0, 0]))
            if fid in mix_fids:
                (mix_ok if delta > 0.0 else mix_broken).append((fid, lxc.name(fid)))
            else:
                (single_ok if delta > 0.0 else single_inv).append((fid, lxc.name(fid)))
        f1.close()
    print(f"Single-kernel MGGAs: {len(single_ok)} respond to xc_func_set_fhc_enforcement(p, 0); {len(single_inv)} FHC-invariant ({[n for _, n in single_inv]})")
    print(f"Mixed/hybrid MGGAs (n_func_aux > 0): {len(mix_ok)} respond, {len(mix_broken)} BROKEN (100% silent no-op!)")
    print("  Sample broken mixed MGGAs:", [nm for _, nm in mix_broken[:20]])


def check_hse06_lcwpbe_molecular():
    print("\n=== 2. HSE06 & LC-wPBE MOLECULAR IMPACT (N2, H2O, CO; def2-TZVP, grid level 4) ===")
    mols = {
        "H2O": gto.M(atom="O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587", basis="def2-tzvp", verbose=0),
        "N2":  gto.M(atom="N 0 0 0; N 0 0 1.0977", basis="def2-tzvp", verbose=0),
        "CO":  gto.M(atom="C 0 0 0; O 0 0 1.1283", basis="def2-tzvp", verbose=0),
    }
    by_name = {lxc.name(i): i for i in lxc.all_ids()}
    Ha1, Ha2, Ha3, Ha4, Ha5 = 0.00979681, 0.0410834, 0.187440, 0.00120824, 0.0347188

    # First get converged PBE0 densities as clean reference densities on each molecule
    for mname, mol in mols.items():
        mf = dft.RKS(mol)
        mf.xc = "pbe0"
        mf.grids.level = 4
        mf.conv_tol = 1e-11
        mf.kernel()
        dm = mf.make_rdm1()
        grids = mf.grids
        ao = dft.numint.eval_ao(mol, grids.coords, deriv=1)
        rho_gga = dft.numint.eval_rho(mol, ao, dm, xctype="GGA")  # (4, N)
        rho = np.ascontiguousarray(rho_gga[0])
        sigma = np.ascontiguousarray(rho_gga[1]**2 + rho_gga[2]**2 + rho_gga[3]**2)
        wt = grids.weights

        mask = rho > 1e-12
        r_m, sig_m, w_m = rho[mask], sigma[mask], wt[mask]
        s_raw = X2S * np.sqrt(sig_m / 4.0) / ((r_m / 2.0) ** (4.0 / 3.0))

        # Check s_scaling_2 occupancy and energy shift from replacing s_bar on s in [1, 15] by smooth s (or CP2K s_bar)!
        # Note: we can mimic removing the s in [1, 15] shift in s_scaling_2 by adjusting sigma -> sigma_unscaled
        # such that s_scaling_2(s_adj) == s_true for s < 8.0!
        frac_s_1_15 = float(np.sum(w_m * r_m * ((s_raw >= 1.0) & (s_raw <= 15.0))) / np.sum(w_m * r_m))
        n_s_gt1 = int(np.sum((s_raw >= 1.0) & (s_raw <= 15.0)))
        n_s_gt15 = int(np.sum(s_raw > 15.0))

        for xcname, fid, omega in [("HSE06", by_name["hyb_gga_xc_hse06"], 0.11),
                                   ("LC-wPBE", by_name["hyb_gga_xc_lc_wpbe"], 0.40)]:
            # Compute aux5 distribution
            s_bar = np.where(
                s_raw < 1.0, s_raw,
                np.where(s_raw > 15.0, 8.572844, s_raw - np.log1p(np.exp(np.minimum(s_raw - 8.572844, 30.0))))
            )
            kF = (3.0 * PI**2 * r_m) ** (1.0 / 3.0)
            w_nu = omega / kF
            eb1 = np.where(w_nu < 14.0, 1.455915450052607, 2.0)
            Hs = (Ha1 * s_bar**2 + Ha2 * s_bar**4) / (1.0 + Ha3 * s_bar**4 + Ha4 * s_bar**5 + Ha5 * s_bar**6)
            aux5 = 9.0 * (s_bar**2 * Hs + eb1 * w_nu**2) / (4.0 * 1.0161144)
            n_jet_zone = int(np.sum((aux5 > 18.0) & (aux5 < 40.0) & (w_nu <= 14.0)))
            n_w14_zone = int(np.sum((w_nu > 13.0) & (w_nu < 15.0)))

            # Exact s_scaling_2 shift vs smooth CP2K s_scaling (s for s<=8.3, smax - sconst/s^2 for s>8.3)
            s_cp2k = np.where(s_raw <= 8.3, s_raw, 8.572844 - 18.79622316 / (s_raw ** 2))
            s_cp2k_clamped = np.minimum(s_cp2k, 8.571228)
            s_inv = s_cp2k_clamped.copy()
            hi_mask = s_cp2k_clamped >= 1.0
            st = s_inv[hi_mask] + np.log1p(np.exp(s_inv[hi_mask] - 8.572844))
            for _ in range(8):
                ex = np.exp(st - 8.572844)
                fv = st - np.log1p(ex) - s_cp2k_clamped[hi_mask]
                dfv = 1.0 / (1.0 + ex)
                st = st - fv / dfv
            s_inv[hi_mask] = st
            sig_cp2k = 4.0 * ((s_inv / X2S) * (r_m / 2.0) ** (4.0 / 3.0)) ** 2

            lxc.L.xc_func_set_ext_params.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]
            if xcname == "HSE06":
                # Exact HSE06 exchange in libxc: +1.0 * gga_x_wpbeh(omega=0) - 0.25 * gga_x_wpbeh(omega=0.11)
                f0 = lxc.Func(by_name["gga_x_wpbeh"], 1)
                lxc.L.xc_func_set_ext_params(ctypes.c_void_p(f0.p), (ctypes.c_double * 1)(0.0))
                f1 = lxc.Func(by_name["gga_x_wpbeh"], 1)
                lxc.L.xc_func_set_ext_params(ctypes.c_void_p(f1.p), (ctypes.c_double * 1)(omega))
                e_lib = f0.exc_vxc({"rho": r_m[:, None], "sigma": sig_m[:, None]})["zk"] - 0.25 * f1.exc_vxc({"rho": r_m[:, None], "sigma": sig_m[:, None]})["zk"]
                e_cp2k = f0.exc_vxc({"rho": r_m[:, None], "sigma": sig_cp2k[:, None]})["zk"] - 0.25 * f1.exc_vxc({"rho": r_m[:, None], "sigma": sig_cp2k[:, None]})["zk"]
                f0.close(); f1.close()
            else:
                # LC-wPBE exchange: +1.0 * gga_x_wpbeh(omega=0.40)
                f1 = lxc.Func(by_name["gga_x_wpbeh"], 1)
                lxc.L.xc_func_set_ext_params(ctypes.c_void_p(f1.p), (ctypes.c_double * 1)(omega))
                e_lib = f1.exc_vxc({"rho": r_m[:, None], "sigma": sig_m[:, None]})["zk"]
                e_cp2k = f1.exc_vxc({"rho": r_m[:, None], "sigma": sig_cp2k[:, None]})["zk"]
                f1.close()
            dE_s_scaling = float(np.sum(w_m * r_m * (e_lib - e_cp2k)))

            print(f"  {mname:3s} | {xcname:7s}: elec_frac(s in [1,15])={frac_s_1_15:.3f} (pts: s[1,15]={n_s_gt1}, s>15={n_s_gt15}), "
                  f"pts(18<aux5<40)={n_jet_zone:4d}, pts(13<w<15)={n_w14_zone:4d}, dE(s_scaling_2 - CP2K)={dE_s_scaling:+.3e} Eh ({dE_s_scaling*627.509:+.4f} kcal/mol)")

    # Compare Vxc matrix (devel vs 7.0.0) on H2O and N2 for LC-wPBE and HSE06!
    code_700 = """
import os, sys, numpy as np
sys.path.insert(0, "/Users/napolesd/Documents/Antigravity/PhysAud")
import lxc
dat = np.load(sys.argv[1])
by_name = {lxc.name(i): i for i in lxc.all_ids()}
f = lxc.Func(by_name[sys.argv[2]], 1)
out = f.exc_vxc({"rho": dat["rho"][:, None], "sigma": dat["sigma"][:, None]})
f.close()
np.savez(sys.argv[3], zk=out["zk"], vrho=out["vrho"][:,0], vsigma=out["vsigma"][:,0])
"""
    lib700 = str(ROOT / ".venv/lib/python3.13/site-packages/pyscf/lib/deps/lib/libxc.15.dylib")
    for mname in ("H2O", "N2"):
        mol = mols[mname]
        mf = dft.RKS(mol); mf.xc = "pbe0"; mf.grids.level = 4; mf.conv_tol = 1e-11; mf.kernel()
        dm = mf.make_rdm1()
        ao = dft.numint.eval_ao(mol, mf.grids.coords, deriv=1)
        rho_gga = dft.numint.eval_rho(mol, ao, dm, xctype="GGA")
        r_m = np.ascontiguousarray(rho_gga[0])
        sig_m = np.ascontiguousarray(rho_gga[1]**2 + rho_gga[2]**2 + rho_gga[3]**2)
        wt = mf.grids.weights
        in_npz = "/tmp/physaud_in.npz"
        out_npz = "/tmp/physaud_out700.npz"
        np.savez(in_npz, rho=r_m, sigma=sig_m)
        for xcname in ("hyb_gga_xc_hse06", "hyb_gga_xc_lc_wpbe"):
            f_dev = lxc.Func(by_name[xcname], 1)
            o_dev = f_dev.exc_vxc({"rho": r_m[:, None], "sigma": sig_m[:, None]})
            f_dev.close()
            subprocess.run([sys.executable, "-c", code_700, in_npz, xcname, out_npz],
                           env=dict(os.environ, LIBXC_PATH=lib700), check=True)
            o_700 = np.load(out_npz)
            dE = float(np.sum(wt * r_m * (o_dev["zk"] - o_700["zk"])))
            dvrho = o_dev["vrho"][:, 0] - o_700["vrho"]
            dvsig = o_dev["vsigma"][:, 0] - o_700["vsigma"]
            dvrho_max = float(np.max(np.abs(dvrho)))
            dVmat = dft.numint.eval_mat(mol, ao, wt, rho_gga, (dvrho, dvsig, None, None), xctype="GGA")
            dV_max = float(np.max(np.abs(dVmat)))
            dE_tr = float(np.einsum("ij,ji->", dVmat, dm))
            print(f"  [devel vs 7.0.0] {mname:3s} | {xcname:20s}: |dExc|={abs(dE):.3e} Eh, max|dvrho|={dvrho_max:.3e}, max|dV_KS|={dV_max:.3e} Eh, |Tr(dV*P)|={abs(dE_tr):.3e} Eh")


def check_r2scanl_molecular():
    print("\n=== 3. r2SCAN-L (pc07_fab zlo = a/40 = 0.044618) MOLECULAR IMPACT (H2O, N2, Ne; def2-TZVP) ===")
    by_name = {lxc.name(i): i for i in lxc.all_ids()}
    a_opt, b_opt = 1.784720, 0.258304
    zlo = a_opt / 40.0

    for mname, atom_str in [
        ("Ne", "Ne 0 0 0"),
        ("H2O", "O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587"),
        ("N2", "N 0 0 0; N 0 0 1.0977"),
    ]:
        mol = gto.M(atom=atom_str, basis="def2-tzvp", verbose=0)
        mf = dft.RKS(mol); mf.xc = "r2scan"; mf.grids.level = 4; mf.conv_tol = 1e-11; mf.kernel()
        dm = mf.make_rdm1()
        ao = dft.numint.eval_ao(mol, mf.grids.coords, deriv=2)
        rho_mgga = dft.numint.eval_rho(mol, ao, dm, xctype="MGGA", with_lapl=True)  # (6, N): rho, gx, gy, gz, lapl, tau
        r_m = np.ascontiguousarray(rho_mgga[0])
        sig_m = np.ascontiguousarray(rho_mgga[1]**2 + rho_mgga[2]**2 + rho_mgga[3]**2)
        lap_m = np.ascontiguousarray(rho_mgga[4])
        tau_m = np.ascontiguousarray(rho_mgga[5])
        wt = mf.grids.weights

        mask = r_m > 1e-12
        r, sig, lap, tau, w = r_m[mask], sig_m[mask], lap_m[mask], tau_m[mask], wt[mask]
        # Compute exact pc07_z on every grid point!
        # Unpolarized: xs = sqrt(sig/4)/(r/2)^(4/3), p = X2S^2 * xs^2, u = (lap/2)/(r/2)^(5/3), q = X2S^2 * u
        rs_s = r / 2.0
        p = (X2S ** 2) * (sig / 4.0) / (rs_s ** (8.0 / 3.0))
        q = (X2S ** 2) * (lap / 2.0) / (rs_s ** (5.0 / 3.0))
        fW = 5.0 * p / 3.0
        Delta = 8.0 * q**2 / 81.0 - p * q / 9.0 + 8.0 * p**2 / 243.0
        GE4 = 1.0 + 5.0 * p / 27.0 + 20.0 * q / 9.0 + Delta
        B = 1.0 - 40.0 * p / 27.0 + 20.0 * q / 9.0
        S = np.sqrt(1.0 + (Delta / (1.0 + fW))**2)
        z_pc = (Delta**2 * (1.0 + 2.0*fW)/(1.0 + fW)**2 + 2.0*Delta*B + B**2 + 2.0*fW*(Delta + B)) / (S * (GE4 + fW * S))

        in_band = (z_pc > 0.005) & (z_pc <= zlo)
        near_cut = np.abs(z_pc - zlo) < 0.005
        elec_below = float(np.sum(w * r * (z_pc <= zlo)) / np.sum(w * r))

        # Exact unclipped alpha_pc07 vs libxc clipped alpha_pc07 (0 for z <= zlo)!
        # Recall r2SCAN-L simply evaluates r2SCAN with alpha_reg replaced by pc07_alpha * K / (K + eta*fW*K) = pc07_alpha / (1 + eta*5*p/3)!
        # Let's evaluate mgga_x_r2scan + mgga_c_r2scan directly at tau corresponding to clipped vs unclipped pc07_alpha!
        pos = z_pc > 1e-4
        fab_smooth = np.zeros_like(z_pc)
        zp = z_pc[pos]
        # Avoid overflow for zp < zlo: (1 + exp(-a/(a-z)))^b / (1 + exp(a/z - a/(a-z)))^b
        arg_pos = a_opt / zp - a_opt / np.maximum(a_opt - zp, 1e-12)
        fab_smooth[pos] = np.where(
            zp >= a_opt, 1.0,
            ((1.0 + np.exp(-a_opt / np.maximum(a_opt - zp, 1e-12))) ** b_opt) * np.exp(-b_opt * np.logaddexp(0.0, arg_pos))
        )
        fab_libxc = np.where(z_pc <= zlo, 0.0, np.where(z_pc >= 39.0*zlo, 1.0, fab_smooth))
        alpha_smooth = np.maximum(z_pc, 0.0) * fab_smooth
        alpha_libxc = np.maximum(z_pc, 0.0) * fab_libxc

        # Feed tau_equiv = tau_W + K_FACTOR_C * (r/2)^(5/3) * 2 * alpha into mgga_x_r2scan (497) + mgga_c_r2scan (498)?
        # Wait: in mgga_x_r2scanl.mpl: alpha_r2scanl = pc07_alpha(p, u) / (1 + eta * 5*p/3), whereas mgga_x_r2scan uses (t - x^2/8)/(K + eta*x^2/8) = (t - K*fW)/(K*(1 + eta*fW)), which is IDENTICAL when t = K*(fW + pc07_alpha)!!
        K_C = 0.3 * (6.0 * PI**2) ** (2.0 / 3.0)
        tau_smooth = 2.0 * K_C * (fW + alpha_smooth) * (rs_s ** (5.0 / 3.0))
        tau_libxc = 2.0 * K_C * (fW + alpha_libxc) * (rs_s ** (5.0 / 3.0))

        fx = lxc.Func(by_name["mgga_x_r2scan"], 1)
        fc = lxc.Func(by_name["mgga_c_r2scan"], 1)
        fx_l = lxc.Func(by_name["mgga_x_r2scanl"], 1)
        fc_l = lxc.Func(by_name["mgga_c_r2scanl"], 1)
        zk_xl = fx_l.exc_vxc({"rho": r[:, None], "sigma": sig[:, None], "lapl": lap[:, None], "tau": tau[:, None]})["zk"]
        zk_cl = fc_l.exc_vxc({"rho": r[:, None], "sigma": sig[:, None], "lapl": lap[:, None], "tau": tau[:, None]})["zk"]
        zk_x_chk = fx.exc_vxc({"rho": r[:, None], "sigma": sig[:, None], "lapl": lap[:, None], "tau": tau_libxc[:, None]})["zk"]
        zk_c_chk = fc.exc_vxc({"rho": r[:, None], "sigma": sig[:, None], "lapl": lap[:, None], "tau": tau_libxc[:, None]})["zk"]
        zk_x_sm = fx.exc_vxc({"rho": r[:, None], "sigma": sig[:, None], "lapl": lap[:, None], "tau": tau_smooth[:, None]})["zk"]
        zk_c_sm = fc.exc_vxc({"rho": r[:, None], "sigma": sig[:, None], "lapl": lap[:, None], "tau": tau_smooth[:, None]})["zk"]
        for fn_obj in (fx, fc, fx_l, fc_l):
            fn_obj.close()

        # Verify identity check: zk_xl + zk_cl == zk_x_chk + zk_c_chk to ~1e-14!
        chk_err = float(np.max(np.abs((zk_xl + zk_cl) - (zk_x_chk + zk_c_chk))))
        dE_cut = float(np.sum(w * r * ((zk_x_chk + zk_c_chk) - (zk_x_sm + zk_c_sm))))
        print(f"  {mname:3s}: elec_frac(z<=zlo)={elec_below:.3f} (pts in (0.005,zlo]={int(np.sum(in_band)):4d}, |z-zlo|<0.005={int(np.sum(near_cut)):4d}), "
              f"map_verify={chk_err:.2e}, dE_xc(zlo=a/40 clip)={dE_cut:+.3e} Eh ({dE_cut*627.509:+.5f} kcal/mol)")


def check_pk09_atomic():
    print("\n=== 4. PK09 (lda_c_pk09) POLE / NaN AT rs ~ 22.74 ON OPEN-SHELL ATOMS (Li, N, O2; def2-TZVP) ===")
    by_name = {lxc.name(i): i for i in lxc.all_ids()}
    for name_sys, atom_str, spin in [("Li", "Li 0 0 0", 1), ("N", "N 0 0 0", 3), ("O2", "O 0 0 0; O 0 0 1.2075", 2)]:
        mol = gto.M(atom=atom_str, basis="def2-tzvp", spin=spin, verbose=0)
        mf = dft.UKS(mol); mf.xc = "lda_x,lda_c_pw"; mf.grids.level = 4; mf.kernel()
        dm = mf.make_rdm1()
        ao = dft.numint.eval_ao(mol, mf.grids.coords, deriv=0)
        ra = dft.numint.eval_rho(mol, ao, dm[0], xctype="LDA")
        rb = dft.numint.eval_rho(mol, ao, dm[1], xctype="LDA")
        rtot = ra + rb
        mask = rtot > 1e-14
        rs_g = ((3.0 / (4.0 * PI)) / rtot[mask]) ** (1.0 / 3.0)
        near = int(np.sum((rs_g > 22.60) & (rs_g < 22.90)))
        f_pk = lxc.Func(by_name["lda_c_pk09"], 2)
        out_pk = f_pk.exc_vxc({"rho": np.ascontiguousarray(np.vstack([ra[mask], rb[mask]]).T)})
        f_pk.close()
        n_nan_zk = int(np.sum(~np.isfinite(out_pk["zk"])))
        n_nan_v = int(np.sum(~np.isfinite(out_pk["vrho"])))
        max_v = float(np.nanmax(np.abs(out_pk["vrho"])))
        # Also try 1 UKS SCF step with lda_x,lda_c_pk09 in PySCF!
        mf_pk = dft.UKS(mol); mf_pk.xc = "lda_x,lda_c_pk09"; mf_pk.grids.level = 4; mf_pk.max_cycle = 5; mf_pk.verbose = 0
        e_pk = mf_pk.kernel()
        print(f"  {name_sys:3s} (spin={spin}): pts in rs=[22.60,22.90]={near:4d}, NaN(zk)={n_nan_zk}, NaN(vrho)={n_nan_v}, max|vrho|={max_v:.3e} Eh, UKS_5cyc_E={e_pk}")


if __name__ == "__main__":
    check_fhc_contract()
    check_hse06_lcwpbe_molecular()
    check_r2scanl_molecular()
    check_pk09_atomic()
