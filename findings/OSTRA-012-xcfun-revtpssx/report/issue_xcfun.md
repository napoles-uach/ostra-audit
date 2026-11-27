TITLE: Errors in seven functionals: REVTPSSX, REVTPSSC, R4SCANX, TFK, TW, SPBEC (formulas/parameters) and SCANX/SCANC edge cases

## Summary
Comparing every XCFun functional exposed through PySCF with libxc, and checking each discrepancy against the original paper and exact constraints, we found errors in seven functionals on master (a89b783). Six have a one- or two-line fix, tested below; the SCAN edge cases (item 7) are reported without a patch. Two built-in tests (R4SCANX, TFK) have reference values generated with the incorrect code, so they pass today and need the new values included in the patch.

| # | Functional | Problem | Effect |
|---|---|---|---|
| 1 | REVTPSSX | registered without `XC_KINETIC`; denominator `(1 + sqrt(e) p)` instead of `(1 + sqrt(e) p)^2`; `p` instead of `p0 = p/|∇n|²` inside the square root | NaN as registered; SCF does not converge |
| 2 | REVTPSSC | `z^2` instead of `z^3` in ε_c = ε_revPKZB (1 + d ε_revPKZB z³) (TPSSC uses `pow3`) | silent, up to 1.35 kcal/mol in atomization energies |
| 3 | R4SCANX | `DX_DAMP4_P = DX_DAMP4_A = 0.232` instead of the published 0.802 and 0.178 | silent, up to 2.0 kcal/mol in atomization energies |
| 4 | TFK | `C_F n^(5/3)` ignores spin polarization | -8 % at ζ = 0.4, -37 % at ζ = 1 |
| 5 | TW | `pow(gaa + gbb, 2.0) / n`: squares quantities that are already squared gradients and drops `gab` | H atom T_W = 0.0075 Eh instead of 0.5 |
| 6 | SPBEC | β and γ swapped (`gamm = 0.066725`, `beta = 0.031091`) and β/γ used in place of γ inside G | E_c 64 % too negative; wrong t → ∞ limit |
| 7 | SCANX, SCANC | ±inf at α = 1 exactly; SCANX v_σ halved at s = 0 | rare points; no patch proposed |

## Details

**1. REVTPSSX** (src/functionals/revtpssx.cpp, revtpssx_eps.hpp). revTPSS keeps the TPSS exchange form with new constants (Perdew et al., PRL 103, 026403 (2009); erratum PRL 106, 179902 changes only tables). Compared with `tpssx_eps.hpp`, three differences: the missing `XC_KINETIC` flag (so hosts treat it as a GGA), the unsquared denominator and `p` used where TPSSX uses `p0`. With the three changes the energy density matches libxc MGGA_X_REVTPSS to 2e-8; the remainder comes from the number of digits of c and e (libxc: 2.35203946, 2.16769874; with those, 1e-15).

**2. REVTPSSC** (revtpssc_eps.hpp:107-109). The paper keeps "the form of TPSS correlation unchanged", ε_c = ε_c^revPKZB [1 + d ε_c^revPKZB z³], d = 2.8. The code uses `pow2`. The error is invisible at z = 0 (uniform gas) and z = 1 (one-electron densities). After the fix: libxc to 1e-6 relative; the rest is β = 0.066725 vs the exact PBE β used by libxc.

**3. R4SCANX** (SCAN_like_eps.hpp:156-157). Furness et al., J. Chem. Phys. 156, 034109 (2022) (arXiv:2110.00647), give d_p2 = 0.361, d_p4 = 0.802, d_α4 = 0.178. #140 (2020) added r4SCAN as "as yet unpublished", so 0.232 is presumably a pre-publication value; the self-test reference values were generated with it. With the published values: libxc MGGA_X_R4SCAN to 2e-15.

**4. TFK** (tfk.cpp). The other kinetic functionals in XCFun (BTK, PW91K, VWK) are spin-resolved; TFK should satisfy T[n_a, n_b] = (T[2n_a] + T[2n_b])/2, i.e. C_F 2^(2/3) (n_a^(5/3) + n_b^(5/3)). The test reference at (1, 0.8), 7.64755771625168 = C_F 1.8^(5/3), is the spin-unaware value; the patch updates it.

**5. TW** (tw.cpp:21). `gaa` and `gbb` are |∇n_a|², |∇n_b|². The fix `0.125 * (gaa/a + gbb/b)` is the spin-scaled von Weizsäcker functional and equals libxc GGA_K_VW and XCFun's own VWK (vonw.cpp). Making TW an alias of VWK would be an alternative. (The form |∇n|²/(8n) of the total density is not spin-scaled: open-shell Li differs from GGA_K_VW by 46 mEh.)

| T_W (Eh) | H atom | Li atom (UKS) | H2O |
|---|---|---|---|
| libxc GGA_K_VW = XCFun VWK | 0.49998 | 7.178193 | 57.28 |
| XCFun TW, master | 0.0075 | 414.1 | 1.3e7 |
| XCFun TW, patched | 0.49998 | 7.178193 | 57.28 |

**6. SPBEC** (spbec.cpp:22-27). Swart, Solà, Bickelhaupt, J. Chem. Phys. 131, 094103 (2009), Eq. (9): H = γφ³ ln[1 + (β/γ) t²/(1 + G t²)], G = (β/γ)/(exp(−ε_c/(γφ³)) − 1), with β ≈ 0.066725 and γ ≈ 0.031091 "kept at their original values". H2O (def2-SVP, PBE density):

| E_c (Eh) | |
|---|---|
| Eq. (9) with VWN5 ε_c | -0.368141 |
| Eq. (9) with PW92 ε_c | -0.365746 |
| libxc GGA_C_SPBE (PW92) | -0.365745 |
| XCFun SPBEC, master | -0.597502 |
| XCFun SPBEC, patched (keeps VWN5) | -0.368140 |

The current code also violates the rapidly varying limit ε_c → 0 for t → ∞ (n = 0.1, |∇n| = 1e4: -4.6e-2 instead of -4e-11); the patch restores it. The patch keeps XCFun's VWN5 LDA part; switching to PW92 would make it identical to libxc. Your call.

**7. SCANX / SCANC.** At α = 1 exactly (uniform-gas-like points), SCANX returns +inf and SCANC −inf (v_σ NaN), where libxc is finite; the switching function involves 1/(1−α). At s = 0, α = 1.2, SCANX v_σ is −0.0508 while libxc and central finite differences of the energy give −0.0993. Such points are rare on molecular grids; we have not prepared a fix.

## How it was tested
- Patch below applied to master a89b783, built with `XCFUN_MAX_ORDER=3`, loaded into PySCF 2.15.1.
- `testall`: all tests ok (with the updated TFK and R4SCANX references; the R4SCANX values were computed with libxc at the test point and match the patched code to 1e-11).
- Against libxc 7.1.2 on H2O, OH and stretched N2 densities (integrated E_xc): REVTPSS 1e-3 mEh, R4SCANX 4e-8 mEh, TFK and TW 1e-8 mEh, SPBEC 2.4 mEh (VWN5 vs PW92, see 6).
- Thirteen unaffected functionals (PBE, B3LYP, TPSS, SCAN, M06-L, PW91, LDA, BP86, B3LYP5, M06, PBE0, CAM-B3LYP, TPSSh) bit-identical, values and first and second derivatives at 2000 random points.
- Impact (PySCF, 26 small molecules, def2-TZVP, released vs patched): atomization energies change by up to 2.0 kcal/mol for r4SCAN and 1.35 kcal/mol for revTPSS correlation, without any warning.

<details><summary>Patch</summary>

```diff
diff --git a/src/functionals/SCAN_like_eps.hpp b/src/functionals/SCAN_like_eps.hpp
index 27613b3..038a049 100644
--- a/src/functionals/SCAN_like_eps.hpp
+++ b/src/functionals/SCAN_like_eps.hpp
@@ -153,8 +153,8 @@ static num SCAN_X_Fx(const num p,
   const parameter CFDX1 = 1.24;
 
   const parameter D_DAMP2 = 0.361;
-  const parameter DX_DAMP4_P = 0.232;
-  const parameter DX_DAMP4_A = 0.232;
+  const parameter DX_DAMP4_P = 0.802;
+  const parameter DX_DAMP4_A = 0.178;
   const parameter B1 = 0.156632;
   const parameter B2 = 0.12083;
   const parameter B3 = 0.5;
diff --git a/src/functionals/r4SCANx.cpp b/src/functionals/r4SCANx.cpp
index ee6feff..3782c31 100644
--- a/src/functionals/r4SCANx.cpp
+++ b/src/functionals/r4SCANx.cpp
@@ -38,11 +38,11 @@ FUNCTIONAL(XC_R4SCANX) = {"r4SCAN exchange functional",
                           1,
                           1e-11,
                           {0.217, 0.0632, 0.191, 0.0535, 0.015, 0.267, 0.0328},
-                          {-1.62076202198e-01,
-                           -8.60008367983e-01,
+                          {-1.62076025565e-01,
+                           -8.60019219236e-01,
                            -5.57023311096e-01,
-                           -3.79265181913e-02,
+                           -3.79374373972e-02,
                            0.00000000000e+00,
                            -1.04606093722e-01,
-                           5.89149583569e-02,
+                           5.89332768915e-02,
                            5.50319523444e-02}};
diff --git a/src/functionals/revtpssc_eps.hpp b/src/functionals/revtpssc_eps.hpp
index c3f2bb0..12f7e2b 100644
--- a/src/functionals/revtpssc_eps.hpp
+++ b/src/functionals/revtpssc_eps.hpp
@@ -104,8 +104,8 @@ static num epsc_revpkzb(const densvars<num> & d) /* eps_c^{rev PKZB} it is the e
 
 template <typename num> static num revtpssc_eps(const densvars<num> & d) {
   num eps_pkzb = epsc_revpkzb(d);
-  num tauwtau2 = pow2(d.gnn / (8.0 * d.n * d.tau));
+  num tauwtau3 = pow3(d.gnn / (8.0 * d.n * d.tau));
   parameter dd = 2.8;
-  return eps_pkzb * (1 + dd * eps_pkzb * tauwtau2);
+  return eps_pkzb * (1 + dd * eps_pkzb * tauwtau3);
 }
 } // namespace revtpssc_eps
diff --git a/src/functionals/revtpssx.cpp b/src/functionals/revtpssx.cpp
index 2dc8a7e..ca6a703 100644
--- a/src/functionals/revtpssx.cpp
+++ b/src/functionals/revtpssx.cpp
@@ -28,5 +28,5 @@ FUNCTIONAL(XC_REVTPSSX) = {
     "for Condensed Matter Physics and Quantum Chemistry\n"
     "Phys. Rev. Lett. 103 (2009) 026403\n"
     "Implemented by Andrea Debnarova\n",
-    XC_DENSITY | XC_GRADIENT,
+    XC_DENSITY | XC_GRADIENT | XC_KINETIC,
     ENERGY_FUNCTION(revtpssx)};
diff --git a/src/functionals/revtpssx_eps.hpp b/src/functionals/revtpssx_eps.hpp
index 5e88898..e12625d 100644
--- a/src/functionals/revtpssx_eps.hpp
+++ b/src/functionals/revtpssx_eps.hpp
@@ -30,7 +30,8 @@ static num x(const num & d_n, const num & d_gnn, const num & d_tau) {
   parameter b = 0.40;
   parameter e = 2.1677;
   parameter c = 2.35204;
-  num p = d_gnn / (4 * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 8.0 / 3.0));
+  num p0 = 1.0 / (4 * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 8.0 / 3.0));
+  num p = d_gnn * p0;
   num tauw = d_gnn / (8.0 * d_n);
   num z = tauw / d_tau;
   num z2 = pow2(z);
@@ -41,11 +42,11 @@ static num x(const num & d_n, const num & d_gnn, const num & d_tau) {
   x_a += 146.0 * pow2(q_b) / 2025.0;
   // Take the gradient out of the square root
   x_a -= 73.0 / 405.0 * q_b * d_gnn *
-         sqrt((0.5 * 0.6 * 0.6) * pow(8 * d_n * d_tau, -2) + 0.5 * p * p);
+         sqrt((0.5 * 0.6 * 0.6) * pow(8 * d_n * d_tau, -2) + 0.5 * p0 * p0);
   // instead of x_a -= 73.0*q_b*sqrt(0.5*0.6*0.6*z2 + 0.5*p*p)/405.0;
   x_a += pow2(10.0 * p / 81.0) / kapa;
   x_a += sqrt(e) * 0.6 * 0.6 * z2 * 20.0 / 81.0 + e * mu * pow3(p);
-  return x_a / (1 + sqrt(e) * p);
+  return x_a / pow2(1 + sqrt(e) * p);
 }
 
 template <typename num>
diff --git a/src/functionals/spbec.cpp b/src/functionals/spbec.cpp
index fd3bf75..4f2d814 100644
--- a/src/functionals/spbec.cpp
+++ b/src/functionals/spbec.cpp
@@ -19,12 +19,12 @@
 /* "Simplified" PBEc, spbec. Implemented by ulfek */
 
 /* These are the constants in the paper, not very high precision.. */
-static const parameter gamm = 0.066725;
-static const parameter beta = 0.031091;
+static const parameter gamm = 0.031091;
+static const parameter beta = 0.066725;
 static const parameter beta_gamma = beta / gamm;
 
 template <typename num> static num G(const num & eps, const num & phi3) {
-  return beta_gamma / expm1(-eps / (beta_gamma * phi3));
+  return beta_gamma / expm1(-eps / (gamm * phi3));
 }
 
 template <typename num>
diff --git a/src/functionals/tfk.cpp b/src/functionals/tfk.cpp
index 9ad7473..1f216ad 100644
--- a/src/functionals/tfk.cpp
+++ b/src/functionals/tfk.cpp
@@ -20,7 +20,7 @@
 template <typename num> static num tfk(const densvars<num> & d) {
   using xcfun_constants::CF;
 
-  return CF * pow(d.n, 5.0 / 3.0);
+  return CF * pow(2.0, 2.0 / 3.0) * (pow(d.a, 5.0 / 3.0) + pow(d.b, 5.0 / 3.0));
 }
 
 FUNCTIONAL(XC_TFK) = {"Thomas-Fermi Kinetic Energy Functional",
@@ -32,9 +32,9 @@ FUNCTIONAL(XC_TFK) = {"Thomas-Fermi Kinetic Energy Functional",
                       1,
                       1e-5,
                       {1., .8},
-                      {7.64755771625168,
-                       7.08107195949229,
-                       7.08107195949229,
-                       2.62261924425641,
-                       2.62261924425641,
-                       2.62261924425641}};
+                      {7.70003416201095,
+                       7.59633312057599,
+                       6.54632143680283,
+                       5.06422208038400,
+                       0.00000000000000,
+                       5.45526786400236}};
diff --git a/src/functionals/tw.cpp b/src/functionals/tw.cpp
index e7f3e57..864d82a 100644
--- a/src/functionals/tw.cpp
+++ b/src/functionals/tw.cpp
@@ -18,7 +18,7 @@
 //  von Weizsacker kinetic energy functional
 
 template <typename num> static num tw(const densvars<num> & d) {
-  return 1. / 8. * pow(d.gaa + d.gbb, 2.0) / d.n;
+  return 0.125 * (d.gaa / d.a + d.gbb / d.b);
 }
 
 FUNCTIONAL(XC_TW) = {
```
</details>

PySCF bundles this XCFun as an alternative backend (`mf._numint.libxc = dft.xcfun`); we will let them know once there is a fix here.

---
Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of electronic-structure libraries. Analysis assisted by AI models (Claude, Anthropic; Gemini, Google); all results were checked and the report reviewed by a human before submission.
