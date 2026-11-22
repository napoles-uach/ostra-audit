import subprocess, sys, math, numpy as np
import mpmath as mp
mp.mp.dps = 50

def exact_Kl(l, z):
    z = mp.mpf(z)
    return float(mp.exp(-z) * mp.sqrt(mp.pi / (2*z)) * mp.besseli(l + mp.mpf("0.5"), z))

print("=== 50-digit mpmath exact K_l(z) reference ===")
for z in [16.0, 16.01, 17.0, 20.0]:
    vals = {l: exact_Kl(l, z) for l in [0, 12, 16, 24, 28, 32]}
    print(f"EXACT z={z:.2f}: " + " ".join(f"l{l}={v:+.10e}" for l, v in vals.items()))

print("\n=== Testing omega_nuc[CART_CUM] stack overflow in PySCF ECPscalar_ipnuc at l=6 (lmax=13 > 12) ===")
code = """
from pyscf import gto
ecp_local = [60, [[-1, [[], [], [[0.5, 5.0]]]]]]
mol = gto.M(atom="Au 0 0 0; H 0 0 2.5", basis={"Au": [[0,[1.0,1.0]]], "H": [[6,[1.0,1.0]]]}, ecp={"Au": ecp_local}, spin=0, unit="Bohr")
ip = mol.intor("ECPscalar_ipnuc")
print("Unexpected success:", ip.shape)
"""
res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
print(f"Subprocess exit code: {res.returncode} (negative = signal, e.g. -6 = SIGABRT stack smash)")
if res.stdout: print("stdout:", res.stdout.strip())
if res.stderr: print("stderr:", res.stderr.strip())
