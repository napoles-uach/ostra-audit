program bitid_tblite
   use, intrinsic :: iso_fortran_env, only : wp => real64, int64
   use mctc_env, only : error_type
   use mctc_io, only : structure_type, new
   use tblite_context_type, only : context_type
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_ceh_ceh, only : get_effective_qat
   implicit none
   type(structure_type) :: mol
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(error_type), allocatable :: err
   integer, parameter :: nat = 3
   real(wp) :: xyz(3, nat), e0, e1, g(3, nat), s(3, 3)
   real(wp) :: cn_en(nat), dcn_endr(3, nat, nat), dcn_endL(3, 3, nat)
   real(wp) :: qat(nat, 1), dqatdr(3, nat, nat, 1), dqatdL(3, 3, nat, 1)
   integer :: i, j, k

   xyz(:, 1) = [0.00_wp, 0.00_wp, 0.12_wp]
   xyz(:, 2) = [0.00_wp, 1.43_wp, -0.95_wp]
   xyz(:, 3) = [0.15_wp, -1.40_wp, -0.98_wp]
   call new(mol, [8, 1, 1], xyz)
   call new_gfn2_calculator(calc, mol, err)
   ctx%verbosity = 0
   call new_wavefunction(wfn, nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   call xtb_singlepoint(ctx, mol, calc, wfn, 1.0e-4_wp, e0)
   call new_wavefunction(wfn, nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   call xtb_singlepoint(ctx, mol, calc, wfn, 1.0e-4_wp, e1, g, s)
   print '(a,2(1x,z16.16))', "E ", transfer(e0, 0_int64), transfer(e1, 0_int64)
   do j = 1, nat; do i = 1, 3
      print '(a,2(i0,1x),z16.16)', "G ", i, j, transfer(g(i, j), 0_int64)
   end do; end do
   do j = 1, 3; do i = 1, 3
      print '(a,2(i0,1x),z16.16)', "S ", i, j, transfer(s(i, j), 0_int64)
   end do; end do

   cn_en = [1.5_wp, 0.8_wp, 0.8_wp]
   dcn_endr = 0.25_wp; dcn_endL = 0.35_wp
   call get_effective_qat(mol, cn_en, qat, dcn_endr, dcn_endL, dqatdr, dqatdL)
   do k = 1, nat; do j = 1, nat; do i = 1, 3
      print '(a,3(i0,1x),z16.16)', "QDR ", i, j, k, transfer(dqatdr(i, j, k, 1), 0_int64)
   end do; end do; end do
end program bitid_tblite
