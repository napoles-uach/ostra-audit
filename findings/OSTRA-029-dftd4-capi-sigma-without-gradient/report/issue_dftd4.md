TITLE: get_dispersion crashes when only one of gradient / sigma is requested (C API: gradient = NULL, sigma != NULL); d4_model couples CN and charge derivatives

## Summary
1. `get_dispersion` and `get_dispersion2` (src/dftd4/disp.f90:89, 207) set `grad = present(gradient) .or. present(sigma)`
   but then write to both arrays unconditionally (`gradient(:, :) = 0.0_wp`, `sigma(:, :) = 0.0_wp`, and the
   `get_dispersion2/3` / `add_coordination_number_derivs` calls). Requesting only one of them dereferences the absent
   argument and crashes. The C wrapper already works around one direction: it allocates `sigma` when only the gradient
   is requested (src/dftd4/api.f90:731-740, with the comment "Would need to refactor dispersion subroutines to make sigma
   truly optional"). The opposite case, `dftd4_get_dispersion(error, mol, disp, param, &energy, NULL, sigma)`, still
   crashes. The Python API and the CLI always pass both and are not affected.
2. `d4_model` / `d4s_model`: `weight_references`, `get_atomic_c6` and `get_polarizabilities` compute derivatives only if
   the CN **and** the charge derivatives are all present (e.g. `if (present(gwdcn) .and. present(gwdq))`,
   src/dftd4/model/d4.f90:270). Requesting only one leaves that `intent(out)` array unwritten. These are public Fortran
   routines; dftd4 itself always requests both.

Checked on main 4a5bbf0 (released 4.3.1 behaves the same).

## Reproducer and results
NH3, PBE0-D4, outputs pre-filled with 999 (programs attached below):

| Call | main |
|---|---|
| C: `dftd4_get_dispersion(..., gradient, sigma)` | correct |
| C: `dftd4_get_dispersion(..., gradient, NULL)` | correct (wrapper allocates sigma) |
| C: `dftd4_get_dispersion(..., NULL, sigma)` | segmentation fault |
| Fortran: `get_dispersion(..., energy, gradient=g)` | segmentation fault |
| Fortran: `d4%weight_references(mol, cn, q, gwvec, gwdcn=a)` | `a` not written (999) |
| Fortran: `d4%weight_references(mol, cn, q, gwvec, gwdq=b)` | `b` not written (999) |

## Suggested fix
- `get_dispersion` / `get_dispersion2`: accumulate into local `gradient_local(3, nat)` / `sigma_local(3, 3)` buffers
  and copy out only the arrays that are present. With that, the workaround in the C wrapper is no longer needed.
- `d4_model` / `d4s_model`: test and write each derivative on its own (`present(gwdcn) .or. present(gwdq)` for the
  branch, `if (present(...))` for each store).

<details><summary>Patch (api.f90, disp.f90, model/d4.f90, model/d4s.f90)</summary>

```diff
diff --git a/src/dftd4/api.f90 b/src/dftd4/api.f90
index 27d06b7..9517bf6 100644
--- a/src/dftd4/api.f90
+++ b/src/dftd4/api.f90
@@ -731,12 +731,6 @@ subroutine get_dispersion_api(verror, vmol, vdisp, vparam, &
    has_sigma = present(c_sigma)
    if (has_sigma) then
       sigma = c_sigma(:3, :3)
-   ! Still needs to be passed into dispersion subroutines,
-   ! just won't be returned through the API.
-   ! Would need to refactor dispersion
-   ! subroutines to make sigma truly optional.
-   else if (has_grad) then
-      allocate(sigma(3,3))
    end if
 
    ! Evaluate energy, gradient (optional), and sigma (optional) analytically
diff --git a/src/dftd4/disp.f90 b/src/dftd4/disp.f90
index 61bc205..c808514 100644
--- a/src/dftd4/disp.f90
+++ b/src/dftd4/disp.f90
@@ -79,6 +79,7 @@ subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma, par
    real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
    real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
    real(wp), allocatable :: lattr(:, :)
+   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)
    type(error_type), allocatable :: local_error
 
    energy = ieee_value(0.0_wp, ieee_quiet_nan)
@@ -123,17 +124,17 @@ subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma, par
       allocate(dEdcn(mol%nat), dEdq(mol%nat))
       dEdcn(:) = 0.0_wp
       dEdq(:) = 0.0_wp
-      gradient(:, :) = 0.0_wp
-      sigma(:, :) = 0.0_wp
+      allocate(gradient_local(3, mol%nat), source=0.0_wp)
+      allocate(sigma_local(3, 3), source=0.0_wp)
    end if
 
    call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
    call param%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, &
-      & disp%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, &
-      & sigma, partition)
+      & disp%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient_local, &
+      & sigma_local, partition)
    if (grad) then
-      call d4_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
-      call d4_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
+      call d4_gemv(dqdr, dEdq, gradient_local, beta=1.0_wp)
+      call d4_gemv(dqdL, dEdq, sigma_local, beta=1.0_wp)
    end if
 
    q(:) = 0.0_wp
@@ -142,16 +143,18 @@ subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma, par
 
    call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
    call param%get_dispersion3(mol, lattr, cutoff%disp3, cutoff%width3, &
-      & disp%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, &
-      & sigma, partition)
+      & disp%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient_local, &
+      & sigma_local, partition)
    if (grad) then
       call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
       call add_coordination_number_derivs(mol, lattr, cutoff%cn, &
-         & disp%rcov, disp%en, dEdcn, gradient, sigma, error=local_error)
+         & disp%rcov, disp%en, dEdcn, gradient_local, sigma_local, error=local_error)
       if (allocated(local_error)) then
          if (present(error)) call move_alloc(local_error, error)
          return
       end if
+      if (present(gradient)) gradient(:, :) = gradient_local(:, :)
+      if (present(sigma)) sigma(:, :) = sigma_local(:, :)
    end if
 
    energy = sum(energies)
@@ -197,6 +200,7 @@ subroutine get_dispersion2(mol, disp, param, cutoff, energy, gradient, sigma, pa
    real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
    real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
    real(wp), allocatable :: lattr(:, :)
+   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)
    type(error_type), allocatable :: local_error
 
    energy = ieee_value(0.0_wp, ieee_quiet_nan)
@@ -241,24 +245,26 @@ subroutine get_dispersion2(mol, disp, param, cutoff, energy, gradient, sigma, pa
       allocate(dEdcn(mol%nat), dEdq(mol%nat))
       dEdcn(:) = 0.0_wp
       dEdq(:) = 0.0_wp
-      gradient(:, :) = 0.0_wp
-      sigma(:, :) = 0.0_wp
+      allocate(gradient_local(3, mol%nat), source=0.0_wp)
+      allocate(sigma_local(3, 3), source=0.0_wp)
    end if
 
    call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
    call param%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, &
-      & disp%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, &
-      & sigma, partition)
+      & disp%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient_local, &
+      & sigma_local, partition)
    if (grad) then
-      call d4_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
-      call d4_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
+      call d4_gemv(dqdr, dEdq, gradient_local, beta=1.0_wp)
+      call d4_gemv(dqdL, dEdq, sigma_local, beta=1.0_wp)
       call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
       call add_coordination_number_derivs(mol, lattr, cutoff%cn, &
-         & disp%rcov, disp%en, dEdcn, gradient, sigma, error=local_error)
+         & disp%rcov, disp%en, dEdcn, gradient_local, sigma_local, error=local_error)
       if (allocated(local_error)) then
          if (present(error)) call move_alloc(local_error, error)
          return
       end if
+      if (present(gradient)) gradient(:, :) = gradient_local(:, :)
+      if (present(sigma)) sigma(:, :) = sigma_local(:, :)
    end if
 
    energy = sum(energies)
diff --git a/src/dftd4/model/d4.f90 b/src/dftd4/model/d4.f90
index 45aaa10..5e962eb 100644
--- a/src/dftd4/model/d4.f90
+++ b/src/dftd4/model/d4.f90
@@ -267,10 +267,10 @@ subroutine weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)
    real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk, wf, zi, gi, maxcn
    real(wp), parameter :: eps_norm = tiny(1.0_wp)**0.5_wp
 
-   if (present(gwdcn) .and. present(gwdq)) then
+   if (present(gwdcn) .or. present(gwdq)) then
       gwvec(:, :, :) = 0.0_wp
-      gwdcn(:, :, :) = 0.0_wp
-      gwdq(:, :, :) = 0.0_wp
+      if (present(gwdcn)) gwdcn(:, :, :) = 0.0_wp
+      if (present(gwdq)) gwdq(:, :, :) = 0.0_wp
 
       !$omp parallel do default(none) schedule(runtime) &
       !$omp shared(gwvec, gwdcn, gwdq, mol, self, cn, q) private(iat, izp, iref, &
@@ -317,18 +317,22 @@ subroutine weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)
             end if
 
             gwvec(iref, iat, 1) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
-            gwdq(iref, iat, 1) = gwk * dzeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
+            if (present(gwdq)) then
+               gwdq(iref, iat, 1) = gwk * dzeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
+            end if
 
             ! This expression behaves differently for -O0 and -O3 optimization
             ! levels for tiny norm values. It yields incorrect values for -O0
             ! due to fp underflow. These tiny values can occur if the CN is very
             ! far from all reference CNs. To ensure consistent behavior, we set
             ! the norm to zero above.
-            dgwk = norm * (expd - expw * dnorm * norm)
-            if (is_exceptional(dgwk) .or. norm == 0.0_wp) then
-               dgwk = 0.0_wp
+            if (present(gwdcn)) then
+               dgwk = norm * (expd - expw * dnorm * norm)
+               if (is_exceptional(dgwk) .or. norm == 0.0_wp) then
+                  dgwk = 0.0_wp
+               end if
+               gwdcn(iref, iat, 1) = dgwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
             end if
-            gwdcn(iref, iat, 1) = dgwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
          end do
       end do
 
@@ -509,11 +513,11 @@ subroutine get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)
    integer :: iat, jat, izp, jzp, iref, jref
    real(wp) :: refc6, dc6, dc6dcni, dc6dcnj, dc6dqi, dc6dqj
 
-   if (present(gwdcn).and.present(dc6dcn) &
-      & .and.present(gwdq).and.present(dc6dq)) then
+   if ((present(gwdcn) .and. present(dc6dcn)) &
+      & .or. (present(gwdq) .and. present(dc6dq))) then
       c6(:, :) = 0.0_wp
-      dc6dcn(:, :) = 0.0_wp
-      dc6dq(:, :) = 0.0_wp
+      if (present(gwdcn) .and. present(dc6dcn)) dc6dcn(:, :) = 0.0_wp
+      if (present(gwdq) .and. present(dc6dq)) dc6dq(:, :) = 0.0_wp
 
       !$omp parallel do default(none) schedule(runtime) &
       !$omp shared(c6, dc6dcn, dc6dq, mol, self, gwvec, gwdcn, gwdq) &
@@ -532,18 +536,26 @@ subroutine get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)
                do jref = 1, self%ref(jzp)
                   refc6 = self%c6(iref, jref, izp, jzp)
                   dc6 = dc6 + gwvec(iref, iat, 1) * gwvec(jref, jat, 1) * refc6
-                  dc6dcni = dc6dcni + gwdcn(iref, iat, 1) * gwvec(jref, jat, 1) * refc6
-                  dc6dcnj = dc6dcnj + gwvec(iref, iat, 1) * gwdcn(jref, jat, 1) * refc6
-                  dc6dqi = dc6dqi + gwdq(iref, iat, 1) * gwvec(jref, jat, 1) * refc6
-                  dc6dqj = dc6dqj + gwvec(iref, iat, 1) * gwdq(jref, jat, 1) * refc6
+                  if (present(gwdcn) .and. present(dc6dcn)) then
+                     dc6dcni = dc6dcni + gwdcn(iref, iat, 1) * gwvec(jref, jat, 1) * refc6
+                     dc6dcnj = dc6dcnj + gwvec(iref, iat, 1) * gwdcn(jref, jat, 1) * refc6
+                  end if
+                  if (present(gwdq) .and. present(dc6dq)) then
+                     dc6dqi = dc6dqi + gwdq(iref, iat, 1) * gwvec(jref, jat, 1) * refc6
+                     dc6dqj = dc6dqj + gwvec(iref, iat, 1) * gwdq(jref, jat, 1) * refc6
+                  end if
                end do
             end do
             c6(iat, jat) = dc6
             c6(jat, iat) = dc6
-            dc6dcn(iat, jat) = dc6dcni
-            dc6dcn(jat, iat) = dc6dcnj
-            dc6dq(iat, jat) = dc6dqi
-            dc6dq(jat, iat) = dc6dqj
+            if (present(gwdcn) .and. present(dc6dcn)) then
+               dc6dcn(iat, jat) = dc6dcni
+               dc6dcn(jat, iat) = dc6dcnj
+            end if
+            if (present(gwdq) .and. present(dc6dq)) then
+               dc6dq(iat, jat) = dc6dqi
+               dc6dq(jat, iat) = dc6dqj
+            end if
          end do
       end do
 
@@ -668,11 +680,11 @@ subroutine get_polarizabilities(self, mol, gwvec, gwdcn, gwdq, alpha, dadcn, dad
    integer :: iat, izp, iref
    real(wp) :: refa, da, dadcni, dadqi
 
-   if (present(gwdcn).and.present(dadcn) &
-      & .and.present(gwdq).and.present(dadq)) then
+   if ((present(gwdcn) .and. present(dadcn)) &
+      & .or. (present(gwdq) .and. present(dadq))) then
       alpha(:) = 0.0_wp
-      dadcn(:) = 0.0_wp
-      dadq(:) = 0.0_wp
+      if (present(gwdcn) .and. present(dadcn)) dadcn(:) = 0.0_wp
+      if (present(gwdq) .and. present(dadq)) dadq(:) = 0.0_wp
 
       !$omp parallel do default(none) schedule(runtime) &
       !$omp shared(alpha, dadcn, dadq, mol, self, gwvec, gwdcn, gwdq) &
@@ -685,12 +697,12 @@ subroutine get_polarizabilities(self, mol, gwvec, gwdcn, gwdq, alpha, dadcn, dad
          do iref = 1, self%ref(izp)
             refa = self%aiw(1, iref, izp)
             da = da + gwvec(iref, iat, 1) * refa
-            dadcni = dadcni + gwdcn(iref, iat, 1) * refa
-            dadqi = dadqi + gwdq(iref, iat, 1) * refa
+            if (present(gwdcn) .and. present(dadcn)) dadcni = dadcni + gwdcn(iref, iat, 1) * refa
+            if (present(gwdq) .and. present(dadq)) dadqi = dadqi + gwdq(iref, iat, 1) * refa
          end do
          alpha(iat) = da
-         dadcn(iat) = dadcni
-         dadq(iat) = dadqi
+         if (present(gwdcn) .and. present(dadcn)) dadcn(iat) = dadcni
+         if (present(gwdq) .and. present(dadq)) dadq(iat) = dadqi
       end do
 
    else
diff --git a/src/dftd4/model/d4s.f90 b/src/dftd4/model/d4s.f90
index f35570a..bf5e3b2 100644
--- a/src/dftd4/model/d4s.f90
+++ b/src/dftd4/model/d4s.f90
@@ -264,10 +264,10 @@ subroutine weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)
    real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk, wf, zi, gi, maxcn
    real(wp), parameter :: eps_norm = tiny(1.0_wp)**0.5_wp
 
-   if (present(gwdcn) .and. present(gwdq)) then
+   if (present(gwdcn) .or. present(gwdq)) then
       gwvec(:, :, :) = 0.0_wp
-      gwdcn(:, :, :) = 0.0_wp
-      gwdq(:, :, :) = 0.0_wp
+      if (present(gwdcn)) gwdcn(:, :, :) = 0.0_wp
+      if (present(gwdq)) gwdq(:, :, :) = 0.0_wp
 
       !$omp parallel do default(none) schedule(runtime) &
       !$omp shared(gwvec, gwdcn, gwdq, mol, self, cn, q) &
@@ -319,13 +319,17 @@ subroutine weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)
                end if
 
                gwvec(iref, iat, jat) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
-               gwdq(iref, iat, jat) = gwk * dzeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
+               if (present(gwdq)) then
+                  gwdq(iref, iat, jat) = gwk * dzeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
+               end if
 
-               dgwk = norm * (expd - expw * dnorm * norm)
-               if (is_exceptional(dgwk) .or. norm == 0.0_wp) then
-                  dgwk = 0.0_wp
+               if (present(gwdcn)) then
+                  dgwk = norm * (expd - expw * dnorm * norm)
+                  if (is_exceptional(dgwk) .or. norm == 0.0_wp) then
+                     dgwk = 0.0_wp
+                  end if
+                  gwdcn(iref, iat, jat) = dgwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
                end if
-               gwdcn(iref, iat, jat) = dgwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
             end do
 
          end do
@@ -517,11 +521,11 @@ subroutine get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)
    integer :: iat, jat, izp, jzp, iref, jref
    real(wp) :: refc6, dc6, dc6dcni, dc6dcnj, dc6dqi, dc6dqj
 
-   if (present(gwdcn).and.present(dc6dcn) &
-      & .and.present(gwdq).and.present(dc6dq)) then
+   if ((present(gwdcn) .and. present(dc6dcn)) &
+      & .or. (present(gwdq) .and. present(dc6dq))) then
       c6(:, :) = 0.0_wp
-      dc6dcn(:, :) = 0.0_wp
-      dc6dq(:, :) = 0.0_wp
+      if (present(gwdcn) .and. present(dc6dcn)) dc6dcn(:, :) = 0.0_wp
+      if (present(gwdq) .and. present(dc6dq)) dc6dq(:, :) = 0.0_wp
 
       !$omp parallel do default(none) schedule(runtime) &
       !$omp shared(c6, dc6dcn, dc6dq, mol, self, gwvec, gwdcn, gwdq) &
@@ -540,18 +544,26 @@ subroutine get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)
                do jref = 1, self%ref(jzp)
                   refc6 = self%c6(iref, jref, izp, jzp)
                   dc6 = dc6 + gwvec(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
-                  dc6dcni = dc6dcni + gwdcn(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
-                  dc6dcnj = dc6dcnj + gwvec(iref, iat, jat) * gwdcn(jref, jat, iat) * refc6
-                  dc6dqi = dc6dqi + gwdq(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
-                  dc6dqj = dc6dqj + gwvec(iref, iat, jat) * gwdq(jref, jat, iat) * refc6
+                  if (present(gwdcn) .and. present(dc6dcn)) then
+                     dc6dcni = dc6dcni + gwdcn(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
+                     dc6dcnj = dc6dcnj + gwvec(iref, iat, jat) * gwdcn(jref, jat, iat) * refc6
+                  end if
+                  if (present(gwdq) .and. present(dc6dq)) then
+                     dc6dqi = dc6dqi + gwdq(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
+                     dc6dqj = dc6dqj + gwvec(iref, iat, jat) * gwdq(jref, jat, iat) * refc6
+                  end if
                end do
             end do
             c6(iat, jat) = dc6
             c6(jat, iat) = dc6
-            dc6dcn(iat, jat) = dc6dcni
-            dc6dcn(jat, iat) = dc6dcnj
-            dc6dq(iat, jat) = dc6dqi
-            dc6dq(jat, iat) = dc6dqj
+            if (present(gwdcn) .and. present(dc6dcn)) then
+               dc6dcn(iat, jat) = dc6dcni
+               dc6dcn(jat, iat) = dc6dcnj
+            end if
+            if (present(gwdq) .and. present(dc6dq)) then
+               dc6dq(iat, jat) = dc6dqi
+               dc6dq(jat, iat) = dc6dqj
+            end if
          end do
       end do
 
@@ -676,11 +688,11 @@ subroutine get_polarizabilities(self, mol, gwvec, gwdcn, gwdq, alpha, dadcn, dad
    integer :: iat, izp, iref
    real(wp) :: refa, da, dadcni, dadqi
 
-   if (present(gwdcn).and.present(dadcn) &
-      & .and.present(gwdq).and.present(dadq)) then
+   if ((present(gwdcn) .and. present(dadcn)) &
+      & .or. (present(gwdq) .and. present(dadq))) then
       alpha(:) = 0.0_wp
-      dadcn(:) = 0.0_wp
-      dadq(:) = 0.0_wp
+      if (present(gwdcn) .and. present(dadcn)) dadcn(:) = 0.0_wp
+      if (present(gwdq) .and. present(dadq)) dadq(:) = 0.0_wp
 
       !$omp parallel do default(none) schedule(runtime) &
       !$omp shared(alpha, dadcn, dadq, mol, self, gwvec, gwdcn, gwdq) &
@@ -693,12 +705,12 @@ subroutine get_polarizabilities(self, mol, gwvec, gwdcn, gwdq, alpha, dadcn, dad
          do iref = 1, self%ref(izp)
             refa = self%aiw(1, iref, izp)
             da = da + gwvec(iref, iat, iat) * refa
-            dadcni = dadcni + gwdcn(iref, iat, iat) * refa
-            dadqi = dadqi + gwdq(iref, iat, iat) * refa
+            if (present(gwdcn) .and. present(dadcn)) dadcni = dadcni + gwdcn(iref, iat, iat) * refa
+            if (present(gwdq) .and. present(dadq)) dadqi = dadqi + gwdq(iref, iat, iat) * refa
          end do
          alpha(iat) = da
-         dadcn(iat) = dadcni
-         dadq(iat) = dadqi
+         if (present(gwdcn) .and. present(dadcn)) dadcn(iat) = dadcni
+         if (present(gwdq) .and. present(dadq)) dadq(iat) = dadqi
       end do
 
    else
```
</details>

| Test | Unpatched | Patched |
|---|---|---|
| C API, sigma only | segmentation fault | identical to the call with both (difference 0.0; 2.7e-20 with OpenMP, reduction order) |
| Fortran `get_dispersion`, gradient only | segmentation fault | identical (0.0) |
| `weight_references`, gwdcn only / gwdq only | not written | identical (0.0) |
| Calls with all outputs or none | — | identical results |
| `meson test`, debug, no OpenMP | 76 OK, 1 expected fail | 76 OK, 1 expected fail |
| `meson test`, release, OpenMP (4 threads) | 75 OK, 1 expected fail | 75 OK, 1 expected fail |

The underlying convention (two optional outputs required together) also appears in mctc-lib and is reported in
grimme-lab/mctc-lib#166. dftd4's own calls pass both derivatives or none, but the public wrapper
`dftd4_ncoord::get_coordination_number` (src/dftd4/ncoord.f90:80) forwards the optional arguments and inherits that
behaviour until mctc-lib is fixed.

<details><summary>test_capi_sigma.c</summary>

```c
#include <stdio.h>
#include <math.h>
#include "dftd4.h"
static double md(double*a,double*b,int n){double m=0;for(int i=0;i<n;i++){double d=fabs(a[i]-b[i]);if(d>m)m=d;}return m;}
int main(void){
  int num[4]={7,1,1,1};
  double xyz[12]={0,0,0.22, 0,1.77,-0.51, 1.53,-0.885,-0.51, -1.53,-0.885,-0.51};
  dftd4_error err=dftd4_new_error();
  dftd4_structure mol=dftd4_new_structure(err,4,num,xyz,NULL,NULL,NULL);
  dftd4_model disp=dftd4_new_d4_model(err,mol);
  dftd4_param par=dftd4_load_rational_damping(err,"pbe0",true);
  double e0,e2,g0[12],s0[9],s2[9];
  dftd4_get_dispersion(err,mol,disp,par,&e0,g0,s0);
  for(int i=0;i<9;i++) s2[i]=999.0;
  dftd4_get_dispersion(err,mol,disp,par,&e2,NULL,s2);
  printf("NH3: max|sigma|(both)=%.3e  max|sigma(sigma only)-sigma(both)|=%.3e  dE=%.1e\n",md(s0,(double[9]){0},9),md(s2,s0,9),e2-e0);
  return 0;}
```
</details>

<details><summary>test_fortran_gradient.f90</summary>

```fortran
program p29f
  use mctc_env, only: wp, error_type
  use mctc_io, only: structure_type, new
  use dftd4, only: d4_model, new_d4_model, rational_damping_param, get_rational_damping, &
     & get_dispersion, realspace_cutoff, damping_param
  implicit none
  type(structure_type) :: mol
  type(d4_model) :: d4
  type(error_type), allocatable :: err
  class(damping_param), allocatable :: par
  real(wp) :: xyz(3,4), e0, e1, g0(3,4), s0(3,3), g1(3,4)
  xyz = reshape([0.0_wp,0.0_wp,0.22_wp, 0.0_wp,1.77_wp,-0.51_wp, 1.53_wp,-0.885_wp,-0.51_wp, -1.53_wp,-0.885_wp,-0.51_wp],[3,4])
  call new(mol, [7,1,1,1], xyz)
  call new_d4_model(err, d4, mol)
  call get_rational_damping("pbe0", par, s9=1.0_wp)
  call get_dispersion(mol, d4, par, realspace_cutoff(), e0, g0, s0)
  g1 = 999.0_wp
  call get_dispersion(mol, d4, par, realspace_cutoff(), e1, gradient=g1)
  print '(a,es10.3)', 'Fortran gradient only: max|g - g(both)| = ', maxval(abs(g1-g0))
end program
```
</details>

<details><summary>test_model_weights.f90</summary>

```fortran
program p33
  use mctc_env, only: wp, error_type
  use mctc_io, only: structure_type, new
  use dftd4_model, only: d4_model, new_d4_model
  implicit none
  type(structure_type) :: mol
  type(d4_model) :: d4
  type(error_type), allocatable :: err
  real(wp) :: xyz(3,4), cn(4), q(4)
  real(wp), allocatable :: gw(:,:,:), gdcn(:,:,:), gdq(:,:,:), a(:,:,:), b(:,:,:)
  integer :: mref
  xyz = reshape([0.0_wp,0.0_wp,0.22_wp, 0.0_wp,1.77_wp,-0.51_wp, 1.53_wp,-0.885_wp,-0.51_wp, -1.53_wp,-0.885_wp,-0.51_wp],[3,4])
  call new(mol, [7,1,1,1], xyz)
  call new_d4_model(err, d4, mol)
  cn = [2.9_wp, 0.95_wp, 0.95_wp, 0.95_wp]; q = [-0.6_wp, 0.2_wp, 0.2_wp, 0.2_wp]
  mref = maxval(d4%ref)
  allocate(gw(mref,4,d4%ncoup)); allocate(gdcn, gdq, a, b, mold=gw)
  call d4%weight_references(mol, cn, q, gw, gdcn, gdq)          ! reference: both
  a = 999.0_wp; call d4%weight_references(mol, cn, q, gw, gwdcn=a)   ! CN derivative only
  b = 999.0_wp; call d4%weight_references(mol, cn, q, gw, gwdq=b)    ! q derivative only
  print '(a,es10.3,a,es10.3,a,es10.3,a,es10.3)', 'max|gwdcn|=',maxval(abs(gdcn)),'  gwdcn-only err=',maxval(abs(a-gdcn)), &
     & '  max|gwdq|=',maxval(abs(gdq)),'  gwdq-only err=',maxval(abs(b-gdq))
end program
```
</details>

---
Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of electronic-structure libraries. Analysis assisted by AI models (Claude, Anthropic; Gemini, Google); all results were checked and the report reviewed by a human before submission.
