TITLE: gCP: parameter tables indexed out of range for Fr/Ra/Ac (E = Inf), Z >= 95 (ERROR STOP) and Kr/Xe/Rn with MINIS; gradient silently skipped when the virial is not requested

## Summary
Two independent problems in the geometric counterpoise (gCP) code, both silent for the caller:

1. `effective_atomic_number` (src/dftd3/gcp/param.f90) maps Z = 72-89 to `Z - 50`, so Fr/Ra/Ac get 37-39, beyond the 36-entry parameter tables. Z > 94 is not mapped at all. `emiss_hf_minis` has 35 entries, so Kr/Xe/Rn (zeff = 36) with HF/MINIS also read past its end.
2. `get_geometric_counterpoise` computes derivatives only if `gradient` **and** `sigma` are present (`grad = present(gradient) .and. present(sigma)`, src/dftd3/gcp.f90:104). With `dftd3_get_counterpoise(..., gradient, NULL)` the gradient stays zero without any message. The dispersion routine handles the two outputs independently.

Checked on main 84cb7cd (v1.7.0).

## Reproducer and results
With the released Python package (dftd3 1.7.0), X-O at 2.2 Angstrom:

```python
import numpy as np
from dftd3.interface import GeometricCounterpoise
B = 1 / 0.52917721
for z, m, b in ((87, "r2scan3c", None), (88, "r2scan3c", None), (89, "r2scan3c", None), (95, "r2scan3c", None)):
    kw = dict(method=m) if b is None else dict(method=m, basis=b)
    r = GeometricCounterpoise(np.array([z, 8]), np.array([[0, 0, 0], [0, 0, 2.2 * B]]), **kw).get_counterpoise(grad=True)
    print(z, r["energy"], abs(r["gradient"]).max())
```

| Element | Method | Result |
|---|---|---|
| Fr, Ra, Ac | r2SCAN-3c | E = inf, gradient NaN |
| Am, ..., Lr (Z >= 95) | r2SCAN-3c | `ERROR STOP invalid shell combination in ssovl` (the host process exits) |
| Kr, Xe | HF/MINIS | E = 0 (read past the 35-entry table; flagged by a bounds-checked build) |

Gradient without the virial (water dimer, HF-3c):

```python
from dftd3 import library
g = GeometricCounterpoise(numbers, positions, method="hf3c")
e = np.array(0.0); grad = np.zeros((6, 3))
library.get_counterpoise(g._mol, g._gcp, e, grad, None)   # sigma = NULL
```

max |grad| = 0.0 instead of 0.0215 Eh/bohr. The Python class, the PySCF plugin and the `s-dftd3` CLI always pass both arrays, so they are not affected; direct C API (and Fortran) callers are.

## Suggested fix
<details><summary>Patch (src/dftd3/gcp/param.f90, src/dftd3/gcp.f90)</summary>

```diff
diff --git a/src/dftd3/gcp.f90 b/src/dftd3/gcp.f90
index 5e8d3a2..bf289d1 100644
--- a/src/dftd3/gcp.f90
+++ b/src/dftd3/gcp.f90
@@ -98,17 +98,21 @@ subroutine get_geometric_counterpoise_atomic(mol, param, cutoff, energies, gradi
    !> Work partition of the atom pairs, absent selects the complete work
    type(work_partition), intent(in), optional :: partition
 
-   real(wp), allocatable :: lattr(:, :)
+   real(wp), allocatable :: lattr(:, :), gtmp(:, :), stmp(:, :)
    logical :: grad
 
-   grad = present(gradient) .and. present(sigma)
+   grad = present(gradient) .or. present(sigma)
+   if (grad) then
+      allocate(gtmp(3, mol%nat), source=0.0_wp)
+      allocate(stmp(3, 3), source=0.0_wp)
+   end if
 
    if (allocated(param%emiss) .and. allocated(param%slater) .and. allocated(param%xv)) then
       call get_lattice_points(mol%periodic, mol%lattice, cutoff%gcp, lattr)
       if (grad) then
          call gcp_deriv(mol, lattr, cutoff%gcp, param%zeff, param%emiss, param%slater, &
             & param%xv, param%rvdw, param%sigma, param%alpha, param%beta, param%damp, &
-            & param%dmp_scal, param%dmp_exp, energies, gradient, sigma, partition)
+            & param%dmp_scal, param%dmp_exp, energies, gtmp, stmp, partition)
       else
          call gcp_energy(mol, lattr, cutoff%gcp, param%zeff, param%emiss, param%slater, &
             & param%xv, param%rvdw, param%sigma, param%alpha, param%beta, param%damp, &
@@ -123,7 +127,7 @@ subroutine get_geometric_counterpoise_atomic(mol, param, cutoff, energies, gradi
    if (param%srb) then
       if (grad) then
          call srb_deriv(mol, lattr, cutoff%srb, mol%num, param%rvdw_srb, param%rscal, param%qscal, &
-            & rexp_srb, zexp_srb, energies, gradient, sigma, partition)
+            & rexp_srb, zexp_srb, energies, gtmp, stmp, partition)
       else
          call srb_energy(mol, lattr, cutoff%srb, mol%num, param%rvdw_srb, param%rscal, param%qscal, &
             & rexp_srb, zexp_srb, energies, partition)
@@ -133,12 +137,18 @@ subroutine get_geometric_counterpoise_atomic(mol, param, cutoff, energies, gradi
    if (param%base) then
       if (grad) then
          call srb_deriv(mol, lattr, cutoff%srb, param%zeff, param%rvdw, param%rscal, param%qscal, &
-            & rexp_base, zexp_base, energies, gradient, sigma, partition)
+            & rexp_base, zexp_base, energies, gtmp, stmp, partition)
       else
          call srb_energy(mol, lattr, cutoff%srb, param%zeff, param%rvdw, param%rscal, param%qscal, &
             & rexp_base, zexp_base, energies, partition)
       end if
    end if
+
+   if (grad) then
+      if (present(gradient)) gradient(:, :) = gradient + gtmp
+      if (present(sigma)) sigma(:, :) = sigma + stmp
+   end if
+
 end subroutine get_geometric_counterpoise_atomic
 
 
diff --git a/src/dftd3/gcp/param.f90 b/src/dftd3/gcp/param.f90
index d091571..861229a 100644
--- a/src/dftd3/gcp/param.f90
+++ b/src/dftd3/gcp/param.f90
@@ -124,7 +124,7 @@ module dftd3_gcp_param
       & 1.232920_wp, 1.343390_wp, 1.448280_wp, 1.613360_wp, 1.768140_wp, 1.992010_wp, 2.233110_wp, 2.493230_wp, &
       & 3.029640_wp, 3.389980_wp, &  ! H-Ca
       & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
-      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
+      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
 
    real(wp), parameter :: emiss_hf_631gd(*) = [ &! H-Ca + Br (no 3d)
       & 0.010083_wp, 0.008147_wp, &
@@ -1059,10 +1059,12 @@ pure function effective_atomic_number(number) result(zeff)
       zeff = number - 18
    case(55:57)
       zeff = number - 18*2
-   case(58:71, 90:94)
+   case(58:71, 89:94)
       zeff = 21
-   case(72:89)
+   case(72:86)
       zeff = number - 2*18 - 14
+   case(87:88)
+      zeff = number - 2*18 - 32
    end select
 end function effective_atomic_number
 
```
</details>

- Fr/Ra are mapped like Cs/Ba (onto K/Ca) and Ac joins the actinides mapped to 21, by analogy with the existing cases; `emiss_hf_minis` gets a 36th (zero) entry. Please check that this mapping is what you intend.
- Z >= 95 is left as is: it may be better to reject unsupported elements in `get_gcp_param` with an error instead of reaching `error stop`.
- `get_geometric_counterpoise` accumulates derivatives in local arrays when either output is requested, and adds them to whichever is present.

Tested by building main with and without the patch (meson, Python extension):

| Test | Unpatched | Patched |
|---|---|---|
| Fr-O, Ra-O, Ac-O (r2SCAN-3c) | inf / NaN | finite; Fr and Ra equal to K and Ca, Ac equal to Th |
| `get_counterpoise(..., sigma=NULL)`, water dimer HF-3c | gradient 0 | identical to the call with sigma |
| Energies, gradients and virials of H2CO and the water dimer with HF-3c, PBEh-3c, r2SCAN-3c, B97-3c, HF/def2-SVP, B3LYP/6-31G* | — | bit-identical |
| `meson test` / Python tests | — | 78 OK, 0 fail / 36 passed |

## Related, minor
`Bint` / `bint_all` choose between the 12-term series and the closed form with `abs(za - zb) < 0.1`, while the argument is `x = (zb - za) R / 2`. Against an mpmath reference, the series has relative errors of 3.5e-4 to 1.8e-3 at x = 5 and the closed form 4e-3 to 0.9 at x = 0.05-0.1 for k = 7-8, so each branch can be used where it is least accurate. In our Cl-Ge test (r2SCAN-3c) the analytic gradient still matches finite differences of the energy to 3e-5 or better, so the effect on energies and gradients seems small; switching on `abs(x)` (and the comment "12 terms seem accurate enough") may be worth a look.

---
Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of electronic-structure libraries. Analysis assisted by AI models (Claude, Anthropic; Gemini, Google); all results were checked and the report reviewed by a human before submission.
