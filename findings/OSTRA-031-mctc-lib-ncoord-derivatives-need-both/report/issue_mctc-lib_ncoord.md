TITLE: get_coordination_number leaves dcndr / dcndL / dcndrlist unset unless both derivatives are requested

## Summary
`ncoord_type%get_coordination_number` (and therefore `get_cn` / `get_cn_legacy`) computes derivatives only when
both optional outputs are present:

```fortran
if (present(dcndr) .and. present(dcndL)) then        ! src/mctc/ncoord/type.f90:222
if (present(dcndrlist) .and. present(dcndL)) then    ! src/mctc/ncoord/type.f90:216
```

A caller asking only for `dcndr` (or only for `dcndL`, or only for `dcndrlist`) gets an `intent(out)` array that is
never written, with no error. When `cut > 0`, `cut_coordination_number` is then called on that unset array and
rescales whatever it contains.

Checked on main 71a21af.

## Reproducer and results
H2O, erf counting function, outputs pre-filled with 999 (program attached below):

| Call | cut = 0 | cut = 1.5 |
|---|---|---|
| `get_cn(..., dcndr=r)` vs reference with both | r stays 999 everywhere | r = 443-653 (rescaled garbage) |
| `get_cn(..., dcndL=l)` vs reference with both | max error 1.0e3 | max error 6.5e2 |

The same happens for `dcndrlist` with a neighbour list. dftd4 (`dftd4_ncoord`) and simple-dftd3 (`dftd3_ncoord`)
wrap this routine and show the same behaviour.

## Suggested fix
Compute into a local buffer when only one of the two derivative arrays is requested, and pass only the
arrays that exist to `cut_coordination_number`:

<details><summary>Patch (src/mctc/ncoord/type.f90)</summary>

```diff
diff --git a/src/mctc/ncoord/type.f90 b/src/mctc/ncoord/type.f90
index db77eac..cc8582f 100644
--- a/src/mctc/ncoord/type.f90
+++ b/src/mctc/ncoord/type.f90
@@ -212,22 +212,40 @@ contains
       !> Derivative of the CN with respect to the Cartesian coordinates
       real(wp), intent(out), optional :: dcndrlist(:, :)
 
+      real(wp), allocatable :: dcndr_loc(:, :, :), dcndL_loc(:, :, :), dcndrlist_loc(:, :)
+
       if (present(list)) then
          if (present(dcndrlist) .and. present(dcndL)) then
             call ncoord_d_list(self, mol, trans, cn, dcndrlist, dcndL, list)
+         else if (present(dcndrlist)) then
+            allocate(dcndL_loc(3, 3, mol%nat))
+            call ncoord_d_list(self, mol, trans, cn, dcndrlist, dcndL_loc, list)
+         else if (present(dcndL)) then
+            allocate(dcndrlist_loc(3, size(list%nlat)))
+            call ncoord_d_list(self, mol, trans, cn, dcndrlist_loc, dcndL, list)
          else
             call ncoord_list(self, mol, trans, cn, list)
          end if
       else
          if (present(dcndr) .and. present(dcndL)) then
             call ncoord_d(self, mol, trans, cn, dcndr, dcndL)
+         else if (present(dcndr)) then
+            allocate(dcndL_loc(3, 3, mol%nat))
+            call ncoord_d(self, mol, trans, cn, dcndr, dcndL_loc)
+         else if (present(dcndL)) then
+            allocate(dcndr_loc(3, mol%nat, mol%nat))
+            call ncoord_d(self, mol, trans, cn, dcndr_loc, dcndL)
          else
             call ncoord(self, mol, trans, cn)
          end if
       end if
 
       if (self%cut > 0.0_wp) then
-         call cut_coordination_number(self%cut, cn, dcndr, dcndL, dcndrlist, list)
+         if (present(list)) then
+            call cut_coordination_number(self%cut, cn, dcndL=dcndL, dcndrlist=dcndrlist, list=list)
+         else
+            call cut_coordination_number(self%cut, cn, dcndr=dcndr, dcndL=dcndL)
+         end if
       end if
 
    end subroutine get_coordination_number
```
</details>

| Test | Unpatched | Patched |
|---|---|---|
| dcndr alone, cut = 0 / 1.5 | not written / rescaled garbage | identical to the call with both (difference 0.0) |
| dcndL alone, cut = 0 / 1.5 | not written / rescaled garbage | identical (0.0) |
| dcndrlist alone (neighbour list) | not written | identical (0.0) |
| Derivatives vs 4-point finite differences (exp, dftd4, erf) | — | agree to ~1e-11 |
| Calls with both or with none | — | bit-identical to unpatched |
| `meson test` | 38 OK | 38 OK, 0 fail |

## Related: the same convention elsewhere in the ecosystem
Requiring two optional outputs together (`present(a) .and. present(b)`) is used widely in mctc-lib, multicharge,
dftd4, simple-dftd3 and tblite. In internal kernels that are always called with both arguments it is harmless. In
public entry points it leaves the requested output unwritten. We found five such entry points and will report the
ones outside mctc-lib in their own repositories:

- multicharge `get_charges` (dqdr without dqdL)
- dftd4 `get_dispersion` / C API `dftd4_get_dispersion` (sigma without gradient: segmentation fault)
- dftd4 `d4_model` / `d4s_model` (CN derivatives without charge derivatives, and the reverse)
- simple-dftd3 `d3_model%weight_references` / `get_atomic_c6` (second derivatives)
- tblite `xtb_singlepoint` (gradient without sigma)

A test that calls each public routine with every combination of optional outputs, starting from sentinel-filled
arrays, catches all of them. We are happy to contribute such tests if that is useful.

<details><summary>test_ncoord.f90</summary>

```fortran
program p31
  use mctc_env, only: wp, error_type
  use mctc_io, only: structure_type, new
  use mctc_ncoord, only: ncoord_type, new_ncoord, cn_count
  implicit none
  type(structure_type) :: mol
  class(ncoord_type), allocatable :: nc
  type(error_type), allocatable :: err
  real(wp) :: xyz(3,3), cn(3), r1(3,3,3), l1(3,3,3), r2(3,3,3), l2(3,3,3), cuts(2)
  integer :: k
  xyz = reshape([0.0_wp,0.0_wp,-0.7_wp, 0.0_wp,1.43_wp,0.4_wp, 0.0_wp,-1.43_wp,0.4_wp],[3,3])
  call new(mol, [8,1,1], xyz)
  cuts = [0.0_wp, 1.5_wp]
  do k = 1, 2
    if (cuts(k) > 0) then
      call new_ncoord(nc, mol, cn_count%erf, err, cut=cuts(k))
    else
      call new_ncoord(nc, mol, cn_count%erf, err)
    end if
    call nc%get_cn(mol, cn, err, dcndr=r1, dcndL=l1)          ! reference: both
    r2 = 999.0_wp; call nc%get_cn(mol, cn, err, dcndr=r2)       ! dcndr only
    l2 = 999.0_wp; call nc%get_cn(mol, cn, err, dcndL=l2)       ! dcndL only
    print '(a,f4.1,a,es10.3,a,es10.3,a,2es11.3)', 'cut=',cuts(k),'  max|dcndr_only-ref|=',maxval(abs(r2-r1)), &
      & '  max|dcndL_only-ref|=',maxval(abs(l2-l1)), '  range dcndr_only:',minval(r2),maxval(r2)
  end do
end program
```
</details>

---
Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of electronic-structure libraries. Analysis assisted by AI models (Claude, Anthropic; Gemini, Google); all results were checked and the report reviewed by a human before submission.
