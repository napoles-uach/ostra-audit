program repro_bint_baux
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use dftd3_cutoff, only : realspace_cutoff
   use dftd3_gcp_param, only : gcp_param, get_gcp_param
   use dftd3_gcp, only : get_geometric_counterpoise, get_geometric_counterpoise_hessian
   implicit none
   type(structure_type) :: mol
   type(gcp_param) :: param
   type(realspace_cutoff) :: cutoff
   real(wp) :: e, grad(3,2), gp1(3,2), gm1(3,2), gp2(3,2), gm2(3,2), sig(3,3), hess(6,6), fd4, h, r0
   integer :: ih, ir
   real(wp), parameter :: hs(4) = [1.0e-2_wp, 1.0e-3_wp, 1.0e-4_wp, 1.0e-5_wp]
   real(wp), parameter :: rs(5) = [0.8_wp, 1.0_wp, 1.5_wp, 2.0_wp, 3.0_wp]

   mol%nat = 2; mol%nid = 2; mol%id = [1, 2]; mol%num = [17, 32]; mol%periodic = .false.
   allocate(mol%xyz(3, 2))
   mol%xyz(:, 1) = 0.0_wp
   cutoff = realspace_cutoff()
   call get_gcp_param(param, mol, method="r2scan3c")
   write(*,'(A,2F10.6,A,F10.6)') 'Cl-Ge r2SCAN-3c slater:', param%slater(1), param%slater(2), &
      ' |dz|=', abs(param%slater(1)-param%slater(2))

   write(*,'(A)') '--- Energy, Gradient, Hessian sweep (Cl-Ge, |dz|=0.1218 >= 0.1) ---'
   do ir = 1, 5
      r0 = rs(ir)
      mol%xyz(:, 2) = [0.0_wp, 0.0_wp, r0]
      e = 0.0_wp; grad = 0.0_wp; sig = 0.0_wp; hess = 0.0_wp
      call get_geometric_counterpoise(mol, param, cutoff, e, grad, sig)
      call get_geometric_counterpoise_hessian(mol, param, cutoff, hess)
      write(*,'(A,F4.2,A,F7.4,A,ES22.14,A,ES22.14,A,ES22.14)') &
         'R=', r0, ' |bx|=', abs(param%slater(1)-param%slater(2))*r0*0.5_wp, &
         ' E=', e, ' Gz=', grad(3,2), ' Hzz=', hess(6,6)
   end do

   write(*,'(A)') '--- FD4 step-size scan of d(grad_z)/dR vs analytical Hzz at R=2.0 Bohr ---'
   r0 = 2.0_wp; mol%xyz(:, 2) = [0.0_wp, 0.0_wp, r0]
   hess = 0.0_wp
   call get_geometric_counterpoise_hessian(mol, param, cutoff, hess)
   write(*,'(A,ES22.14)') 'Analytical Hzz:', hess(6,6)
   do ih = 1, 4
      h = hs(ih)
      e = 0.0_wp; gp1 = 0.0_wp; gm1 = 0.0_wp; gp2 = 0.0_wp; gm2 = 0.0_wp; sig = 0.0_wp
      mol%xyz(3,2) = r0 + h;   call get_geometric_counterpoise(mol, param, cutoff, e, gp1, sig)
      mol%xyz(3,2) = r0 - h;   call get_geometric_counterpoise(mol, param, cutoff, e, gm1, sig)
      mol%xyz(3,2) = r0 + 2*h; call get_geometric_counterpoise(mol, param, cutoff, e, gp2, sig)
      mol%xyz(3,2) = r0 - 2*h; call get_geometric_counterpoise(mol, param, cutoff, e, gm2, sig)
      fd4 = (-gp2(3,2) + 8.0_wp*gp1(3,2) - 8.0_wp*gm1(3,2) + gm2(3,2)) / (12.0_wp * h)
      write(*,'(A,ES8.1,A,ES22.14,A,ES11.3,A,ES11.3)') &
         '  h=', h, ' FD4=', fd4, ' abs_diff=', abs(hess(6,6)-fd4), ' rel=', abs(hess(6,6)-fd4)/abs(fd4)
   end do
end program repro_bint_baux
