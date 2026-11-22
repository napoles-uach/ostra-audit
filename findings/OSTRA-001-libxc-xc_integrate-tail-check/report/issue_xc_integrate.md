**Title:** `xc_integrate`: tail-check cut-off truncates integrals whose integrand is negligible at mid-interval (breaks `LDA_X_1D_SOFT` for β·k_F ≳ 38.6)

### Summary

The early-termination test in the tanh-sinh integrator (`src/integrate.c`) stops the abscissa loop as soon as a single term `w·(|f(x₊)|+|f(x₋)|)` falls below tolerance. Because the loop runs from the interval midpoint outwards, an integrand that is negligible at the midpoint but concentrated near an endpoint is truncated before the endpoint region is ever sampled. The result is silently wrong by orders of magnitude.

This currently affects `LDA_X_1D_SOFT`, whose integrands `2K₀(t)` and `2t·K₀(t)` on `[0, R]` (R = πβρ) decay as e⁻ᵗ. For R/2 ≳ 38.6, i.e. β·k_F ≳ 38.6, the exchange energy collapses and even turns positive (figure 1).

- **Affected versions:** current `devel` (`506ec7ba`; tested at `caa72bd`, and `src/integrate.c` is unchanged between the two) and release 7.1.2 (the PySCF-bundled build shows identical numbers).
- **Build:** source distribution, CMake, `-DCMAKE_BUILD_TYPE=Release -DMAXORDER=1 -DBUILD_SHARED_LIBS=ON`, GCC on Linux (aarch64). I also reproduced it with the integrator compiled standalone (`gcc -O2`).
- **Responsible lines:** the tail checks in `src/integrate.c` (level-0 loop, around line 114; refinement loop, around line 141).

### Reproducer

```c
#include <stdio.h>
#include <xc.h>
#include <xc_funcs.h>
int main(void) {
  double rho[] = {10.0, 24.5, 24.7, 30.0, 100.0, 1000.0}, zk[6], vrho[6];
  xc_func_type f;
  xc_func_init(&f, XC_LDA_X_1D_SOFT, XC_UNPOLARIZED);
  xc_lda_exc_vxc(&f, 6, rho, zk, vrho);
  for (int i = 0; i < 6; i++) printf("rho = %8.2f  eps_x = % .10f  v_x = % .10f\n", rho[i], zk[i], vrho[i]);
  xc_func_end(&f);
}
```

Current `devel`:
```
rho =    24.50  eps_x = -0.4958644415  v_x = -0.5000000000
rho =    24.70  eps_x =  0.0031255097  v_x = -0.0009765625   <- positive exchange energy
rho =    30.00  eps_x = -0.0004815794  v_x = -0.0004881758
rho =   100.00  eps_x = -0.0000193063  v_x = -0.0000193212
rho =  1000.00  eps_x = -0.0000000000  v_x = -0.0000000000
```
The expected high-density limit is ε_x → −1/(2β) = −0.5 (β = 1). The threshold is the same in β·k_F for β = 0.5, 1 and 2.

Calling `xc_integrate` directly on ∫₀ᴿ 2K₀(t) dt (exact value → π) gives:

| R | `xc_integrate` | exact |
|---|---|---|
| 77.0 | 3.141592653589 | 3.141592653590 |
| 77.5 | 0.006135923152 (= π/512) | 3.141592653590 |
| 100 | 0.003067019959 (= π/1024) | 3.141592653590 |

I got the same numbers by compiling `integrate.c` standalone with an independent K₀, so the problem lies in the quadrature and not in `xc_bessel_K0`.

### Root cause

In each refinement level the first new abscissa (k = 1, t = h) lies next to the midpoint x = R/2, where 2K₀(R/2) ~ e^(−R/2). Once that falls below `tol·(1+|T|)`, the `break` fires at k = 1 and no term near x = 0 is added. `new_terms` is then ≈ 0 and `T = T/2 + h·new_terms` only halves the previous estimate. This produces the exact π·2⁻ᵐ plateaus, and the convergence test accepts the result (figure 2).

The comment above the test already describes the intended condition ("once the *weight* is well below relative tolerance"), but the code tests the product weight × integrand.

### Proposed fix

Require the weight itself to be negligible before stopping:

```diff
--- a/src/integrate.c
+++ b/src/integrate.c
@@ -111,7 +111,7 @@
     T += w * (pts[0] + pts[1]);
     /* tail check: once the weight is well below relative tolerance,
        further terms cannot move the integral inside the target. */
-    if(w * (fabs(pts[0]) + fabs(pts[1])) <= tol * (1.0 + fabs(T)))
+    if(w <= tol && w * (fabs(pts[0]) + fabs(pts[1])) <= tol * (1.0 + fabs(T)))
       break;
   }
   T *= h;
@@ -138,7 +138,7 @@
       pts[1] = mid - xt;
       func(pts, 2, ex);
       new_terms += w * (pts[0] + pts[1]);
-      if(w * (fabs(pts[0]) + fabs(pts[1])) <= tol * (1.0 + fabs(T)))
+      if(w <= tol && w * (fabs(pts[0]) + fabs(pts[1])) <= tol * (1.0 + fabs(T)))
         break;
     }
     T = 0.5 * T + h * new_terms;
```

Since w ≤ 2ε is only reached at t ≈ 3.0 (just below `XC_TANH_SINH_T_MAX` ≈ 3.13), the loop at most continues to the existing t_max guard. No abscissa outside the current range is evaluated.

### Verification

- **Standalone integrator.** ∫₀ᴿ 2K₀ = π and ∫₀ᴿ 2tK₀ = 2(1 − R·K₁(R)) are reproduced to < 10⁻¹¹ relative error for R up to 10⁴. Without the fix the result collapses at R ≈ 77.2 (figure 3).
- **Non-regression on all `xc_integrate` callers.** I compared ε and the potential before and after the patch at 2000 random points (ρ ∈ [10⁻⁶, 10³], s ∈ [10⁻³, 30], both spin channels, ζ ∈ [−0.9, 0.9]) for `gga_c_ft97`, `gga_x_ft97_a`, `gga_x_ft97_b`, `gga_x_fd_lb94`, `gga_x_fd_revlb94` and `lda_x_1d_exponential`. The results are **bit-identical** (largest difference 1.3·10⁻¹⁶). `lda_x_1d_soft` changes only above the threshold (ρ ≥ 24.6, or ρ_total ≥ 14.1 in the polarized case), where it now follows the exact result, e.g. ε_x(ρ = 1000) = −0.4998986788.
- **Potential vs. energy consistency.** For `lda_x_1d_soft` with the fix, v_x matches a finite-difference derivative of ρ·ε_x for both spin channels (previously inconsistent above the threshold).
- **Smooth integrands.** The number of evaluations is unchanged for e^(−x²) on [0, 3] and −ln x on [0, 1].

### Broader test of the integrator

I ran a battery of integrands with known integrals through the standalone `xc_integrate`, with and without the patch:

| integrand | current | patched |
|---|---|---|
| e⁻ˣ on [0, R], R = 80 … 1000 | rel. error ≈ 1 | ≤ 2·10⁻¹³ |
| e^(−(R−x)) on [0, R] (mass at the right end) | rel. error ≈ 1 for R ≥ 80 | ≤ 2·10⁻¹³ |
| x⁵e⁻ˣ on [0, R], R ≥ 300 | rel. error ≈ 1 | ≤ 2·10⁻¹⁵ |
| narrow Gaussian at 2 % of [0, R], **already R = 10** | rel. error 4 %, up to 100 % for larger R | ≤ 2·10⁻¹⁵ |
| −e⁻ˣ ln x on [0, 100] (= γ) | rel. error 0.58 | 6·10⁻¹³ |
| e^(−x²) on [−L, L], −ln x, 1/(1+x²), √(1−x²), cos 10x, 1/√x | unchanged | unchanged |

So the problem is not specific to long intervals. Any integrand whose mass sits away from the midpoint can be truncated. The patch never made a result worse in these tests.

**Cost.** For the current callers the runtime is unchanged (FT97, FD-LB94, 1D-exponential) or about +50 % (1D-soft, which now actually integrates the high-density region). The worst case I found is a narrow integrand centred on a wide interval (e^(−x²) on [−70, 70]): 77 → 3211 evaluations, because the loop now always runs out to w ≤ tol. If that matters, an alternative is to stop only after the weight *and* several consecutive terms are negligible. I chose the minimal change.

### Practical impact

For β = 1 the failure needs ρ > 24.6 bohr⁻¹, which typical 1D model systems do not reach, so the effect on published results is probably small. The issue is in the shared integrator, though, so any future caller with an endpoint-concentrated integrand on a long interval would be affected silently.

---
*Disclosure: the issue was found and analysed with the help of an AI assistant (Claude). It ran a systematic check of potential vs. energy consistency across all functionals, then isolated and tested the fix. All the numbers above come from the reproducer and from scripts I can share.*
