TITLE: get_charges leaves dqdr unset when requested without dqdL; no error for elements without EEQ parameters (Z > 103)

## Summary
1. `get_charges` computes charge derivatives only if `dqdr` **and** `dqdL` are present (`grad = (present(dqdr) .and. present(dqdL)) .or. present(cache)`, src/multicharge/charge.f90:67). A caller asking only for the Cartesian derivatives gets an `intent(out)` array that is never written. The same pattern is used for `gradient`/`sigma` and `dqdr`/`dqdL` in `solve` (model/type.F90:181-183).
2. The EEQ2019 and EEQBC2025 tables stop at Z = 103; the getters return -1 beyond that, and `new_eeq2019_model` / `new_eeqbc2025_model` do not check. Superheavy elements then get unphysical charges without an error. dftd4 rejects only Z = 104-111 (src/dftd4/model/d4.f90:106), so Cn-Og reach this path through dftd4.

Checked on main aa2daad (after #81). The same two-optional-outputs convention in mctc-lib is reported in grimme-lab/mctc-lib#166.

## Reproducer and results
Fortran program against the library (attached below):

| Test | Result on main |
|---|---|
| H2O, `get_charges(..., dqdr=dqdr)` vs `get_charges(..., dqdr=dqdr, dqdL=dqdL)` | dqdr keeps whatever it held (999 in the test) |
| Cn-O, Og-O (3.8 bohr), EEQ2019 | q(Cn) = -0.89, q(O) = +0.89; q(Og) = -0.88; no error |
| Same with dftd4 4.3.1 (Python) | q(Cn) = -0.87, q(O) = +0.87; dispersion energy returned |

## Suggested fix
<details><summary>Patch (charge.f90, model/type.F90, param.f90)</summary>

```diff
diff --git a/src/multicharge/charge.f90 b/src/multicharge/charge.f90
index a0d9022..cc9de18 100644
--- a/src/multicharge/charge.f90
+++ b/src/multicharge/charge.f90
@@ -64,7 +64,7 @@ subroutine get_charges(mchrg_model, mol, error, qvec, dqdr, dqdL, cache)
    real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
    real(wp), allocatable :: trans(:, :)
 
-   grad = (present(dqdr) .and. present(dqdL)) .or. present(cache)
+   grad = present(dqdr) .or. present(dqdL) .or. present(cache)
 
    allocate(cn(mol%nat), qloc(mol%nat))
    if (grad) then
diff --git a/src/multicharge/model/type.F90 b/src/multicharge/model/type.F90
index c020488..6ae7a25 100644
--- a/src/multicharge/model/type.F90
+++ b/src/multicharge/model/type.F90
@@ -179,8 +179,8 @@ subroutine solve(self, mol, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
 
    ! Calculate gradient if the respective arrays are present
    dcn = present(dcndr) .and. present(dcndL)
-   grad = present(gradient) .and. present(sigma) .and. dcn
-   cpq = ((present(dqdr) .and. present(dqdL)) .or. present(response)) .and. dcn
+   grad = (present(gradient) .or. present(sigma)) .and. dcn
+   cpq = (present(dqdr) .or. present(dqdL) .or. present(response)) .and. dcn
 
    cache => local_cache
    if (present(response)) then
@@ -268,11 +268,15 @@ subroutine solve(self, mol, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
    end if
 
    if (grad) then
-      gradient = 0.0_wp
-      call gemv(dadr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=0.5_wp)
-      call gemv(dxdr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=-1.0_wp)
-      call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
-      call gemv(dxdL, vrhs, sigma, beta=1.0_wp, alpha=-1.0_wp)
+      if (present(gradient)) then
+         gradient = 0.0_wp
+         call gemv(dadr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=0.5_wp)
+         call gemv(dxdr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=-1.0_wp)
+      end if
+      if (present(sigma)) then
+         call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
+         call gemv(dxdL, vrhs, sigma, beta=1.0_wp, alpha=-1.0_wp)
+      end if
    end if
 
    if (cpq) then
@@ -288,7 +292,7 @@ subroutine solve(self, mol, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
          if (present(dqdr)) dqdr(:, :, :) = cache%dqdr
          call move_alloc(vrhs, cache%qvec)
          call move_alloc(ainv, cache%ainv)
-      else
+      else if (present(dqdr)) then
          call gemm(dadr, ainv(:, :mol%nat), dqdr, alpha=-1.0_wp)
       end if
    end if
diff --git a/src/multicharge/param.f90 b/src/multicharge/param.f90
index 605e3cd..eccd583 100644
--- a/src/multicharge/param.f90
+++ b/src/multicharge/param.f90
@@ -15,7 +15,7 @@
 
 module multicharge_param
    use mctc_data, only: get_covalent_rad, get_pauling_en, get_vdw_rad
-   use mctc_env, only: error_type, wp
+   use mctc_env, only: error_type, fatal_error, wp
    use mctc_io, only: structure_type
    use mctc_io_convert, only: autoaa
    use multicharge_model, only: mchrg_model_type, &
@@ -64,6 +64,12 @@ subroutine new_eeq2019_model(mol, model, error)
    rad = get_eeq_rad(mol%num)
    rcov = get_covalent_rad(mol%num)
 
+   if (any(rad <= 0.0_wp)) then
+      call fatal_error(error, "EEQ2019 parameters are not available for element '"// &
+         & trim(mol%sym(minloc(rad, 1)))//"'")
+      return
+   end if
+
    allocate(eeq)
    call new_eeq_model(eeq, mol=mol, error=error, chi=chi, &
       & rad=rad, eta=eta, kcnchi=kcnchi, cutoff=cutoff, &
@@ -102,6 +108,11 @@ subroutine new_eeqbc2025_model(mol, model, error)
    cap = get_eeqbc_cap(mol%num)
    rcov = get_eeqbc_cov_radii(mol%num)
    avg_cn = get_eeqbc_avg_cn(mol%num)
+   if (any(rad <= 0.0_wp)) then
+      call fatal_error(error, "EEQBC2025 parameters are not available for element '"// &
+         & trim(mol%sym(minloc(rad, 1)))//"'")
+      return
+   end if
    ! Electronegativities normalized to Fluorine
    ! with actinides (Th-Lr) set to average of 1.30
    en = get_pauling_en(mol%num)
```
</details>

Each optional output is handled on its own, and the two constructors return an error naming the element when a parameter is missing. Tested by building main with and without the patch:

| Test | Unpatched | Patched |
|---|---|---|
| dqdr requested alone | not written | equal to dqdr from the call with dqdL |
| Cn-O, Og-O | charges with the wrong sign, no error | `EEQ2019 parameters are not available for element 'Cn'` |
| Hg-O, Lr-O charges | — | unchanged |
| `meson test` | 50 OK | 50 OK, 0 fail |

dftd4 may also want to extend its element check to Z > 111 or pass the multicharge error through.

<details><summary>test_multicharge.f90</summary>

```fortran
program t
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, new
   use multicharge, only : mchrg_model_type, new_eeq2019_model, get_charges
   implicit none
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(error_type), allocatable :: error
   real(wp), allocatable :: q(:), dqdr(:,:,:), dqdL(:,:,:), dqdr2(:,:,:)
   real(wp) :: xyz(3,3)
   ! OSTRA-026: water, dqdr requested alone vs together with dqdL
   xyz = reshape([0.0_wp,0.0_wp,0.221_wp, 0.0_wp,1.430_wp,-0.886_wp, 0.0_wp,-1.430_wp,-0.886_wp],[3,3])
   call new(mol, [8,1,1], xyz)
   call new_eeq2019_model(mol, model, error)
   allocate(q(3), dqdr(3,3,3), dqdL(3,3,3), dqdr2(3,3,3))
   call get_charges(model, mol, error, q, dqdr=dqdr, dqdL=dqdL)
   dqdr2 = 999.0_wp
   call get_charges(model, mol, error, q, dqdr=dqdr2)
   print '(a,es10.3,a,es10.3)', "026 max|dqdr| (with dqdL) = ", maxval(abs(dqdr)), &
      & "   max|dqdr(alone) - dqdr(with dqdL)| = ", maxval(abs(dqdr2 - dqdr))
   ! OSTRA-021: element beyond the EEQ2019 parameter range
   call run(80); call run(103); call run(112); call run(118)
contains
   subroutine run(zz)
      integer, intent(in) :: zz
      type(structure_type) :: m2
      real(wp) :: x2(3,2), q2(2)
      class(mchrg_model_type), allocatable :: md
      type(error_type), allocatable :: err
      x2 = reshape([0.0_wp,0.0_wp,0.0_wp, 0.0_wp,0.0_wp,3.8_wp],[3,2])
      call new(m2, [zz, 8], x2)
      call new_eeq2019_model(m2, md, err)
      if (allocated(err)) then
         print '(a,i4,2a)', "021 Z=", zz, "  error: ", err%message; return
      end if
      call get_charges(md, m2, err, q2)
      if (allocated(err)) then
         print '(a,i4,2a)', "021 Z=", zz, "  error: ", err%message; return
      end if
      print '(a,i4,a,2f9.4)', "021 Z=", zz, "  q(Z), q(O) = ", q2
   end subroutine run
end program t
```
</details>

---
Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of electronic-structure libraries. Analysis assisted by AI models (Claude, Anthropic; Gemini, Google); all results were checked and the report reviewed by a human before submission.
