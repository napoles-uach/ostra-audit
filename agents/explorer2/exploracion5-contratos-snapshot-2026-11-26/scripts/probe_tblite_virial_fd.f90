program probe_tblite_virial_fd
   use, intrinsic :: iso_fortran_env, only : wp => real64
   use mctc_env, only : error_type
   use mctc_io, only : structure_type, new
   use tblite_context_type, only : context_type
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none

   type(structure_type) :: mol, mol_p
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(error_type), allocatable :: err
   integer, parameter :: nat = 3
   integer :: num(nat)
   real(wp) :: xyz(3, nat), e_ref, g_ref(3, nat), s_ref(3, 3)
   real(wp) :: s_only(3, 3), e_only
   real(wp) :: ep1, em1, ep2, em2, eps(3, 3), h, fd_val, max_fd_s
   integer :: a, b, ih
   real(wp), parameter :: steps(3) = [1.0e-3_wp, 1.0e-4_wp, 1.0e-5_wp]

   num = [8, 1, 1]
   xyz(:, 1) = [0.00_wp, 0.00_wp, 0.12_wp]
   xyz(:, 2) = [0.00_wp, 1.43_wp, -0.95_wp]
   xyz(:, 3) = [0.15_wp, -1.40_wp, -0.98_wp]
   call new(mol, num, xyz)
   call new_gfn2_calculator(calc, mol, err)
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   ctx%verbosity = 0
   call xtb_singlepoint(ctx, mol, calc, wfn, 1.0e-4_wp, e_ref, g_ref, s_ref)

   s_only = 999.0_wp
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   call xtb_singlepoint(ctx, mol, calc, wfn, 1.0e-4_wp, e_only, sigma=s_only)

   print '(a,es15.8,a,es12.5,a,es12.5)', "GFN2 virial oracle: E=", e_ref, &
      & " max|s_ref|=", maxval(abs(s_ref)), " max|s_only-s_ref|=", maxval(abs(s_only - s_ref))

   do ih = 1, size(steps)
      h = steps(ih)
      max_fd_s = 0.0_wp
      do a = 1, 3
         do b = 1, 3
            eps = 0.0_wp; eps(a, b) = h
            call new(mol_p, num, xyz + matmul(eps, xyz))
            call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
            call xtb_singlepoint(ctx, mol_p, calc, wfn, 1.0e-4_wp, ep1)
            eps = 0.0_wp; eps(a, b) = -h
            call new(mol_p, num, xyz + matmul(eps, xyz))
            call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
            call xtb_singlepoint(ctx, mol_p, calc, wfn, 1.0e-4_wp, em1)
            eps = 0.0_wp; eps(a, b) = 2.0_wp * h
            call new(mol_p, num, xyz + matmul(eps, xyz))
            call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
            call xtb_singlepoint(ctx, mol_p, calc, wfn, 1.0e-4_wp, ep2)
            eps = 0.0_wp; eps(a, b) = -2.0_wp * h
            call new(mol_p, num, xyz + matmul(eps, xyz))
            call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
            call xtb_singlepoint(ctx, mol_p, calc, wfn, 1.0e-4_wp, em2)
            fd_val = (-ep2 + 8.0_wp*ep1 - 8.0_wp*em1 + em2) / (12.0_wp * h)
            max_fd_s = max(max_fd_s, abs(fd_val - s_ref(a, b)))
         end do
      end do
      print '(a,es8.1,a,es11.3)', "  4pt-FD virial h=", h, " max|s_ref-FD|=", max_fd_s
   end do
end program probe_tblite_virial_fd
