program repro_gcp_grad_contract
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use dftd3_cutoff, only : realspace_cutoff
   use dftd3_gcp_param, only : gcp_param, get_gcp_param
   use dftd3_gcp, only : get_geometric_counterpoise
   implicit none
   type(structure_type) :: mol
   type(gcp_param) :: param
   type(realspace_cutoff) :: cutoff
   real(wp) :: e1, e2, ep, em, h, fd, grad_with_sig(3, 2), grad_no_sig(3, 2), sigma(3, 3)

   mol%nat = 2; mol%nid = 2; mol%id = [1, 2]; mol%num = [6, 7]; mol%periodic = .false.
   allocate(mol%xyz(3, 2))
   mol%xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   mol%xyz(:, 2) = [0.0_wp, 0.0_wp, 2.5_wp]
   cutoff = realspace_cutoff()
   call get_gcp_param(param, mol, method="r2scan3c")

   e1 = 0.0_wp; grad_with_sig = 0.0_wp; sigma = 0.0_wp
   call get_geometric_counterpoise(mol, param, cutoff, e1, gradient=grad_with_sig, sigma=sigma)

   ! Typical non-periodic call: pass gradient, omit sigma (identical to C API sigma=NULL)
   e2 = 0.0_wp; grad_no_sig = -999.0_wp
   call get_geometric_counterpoise(mol, param, cutoff, e2, gradient=grad_no_sig)

   ! Also test with zero-initialized gradient (what C API dftd3_get_counterpoise sees!)
   block
      real(wp) :: e3, grad_zero(3, 2)
      e3 = 0.0_wp; grad_zero = 0.0_wp
      call get_geometric_counterpoise(mol, param, cutoff, e3, gradient=grad_zero)
      h = 1.0e-5_wp
      ep = 0.0_wp; mol%xyz(3,2) = 2.5_wp + h; call get_geometric_counterpoise(mol, param, cutoff, ep)
      em = 0.0_wp; mol%xyz(3,2) = 2.5_wp - h; call get_geometric_counterpoise(mol, param, cutoff, em)
      fd = (ep - em) / (2.0_wp * h)
      write(*,'(A,ES18.10)') 'Energy (with sigma):     ', e1
      write(*,'(A,ES18.10)') 'Energy (without sigma):  ', e2
      write(*,'(A,ES18.10)') 'Grad_z (FD reference):   ', fd
      write(*,'(A,ES18.10)') 'Grad_z (with sigma):     ', grad_with_sig(3,2)
      write(*,'(A,ES18.10)') 'Grad_z (no sigma, -999): ', grad_no_sig(3,2)
      write(*,'(A,ES18.10)') 'Grad_z (no sigma, 0.0):  ', grad_zero(3,2)
      write(*,'(A,ES12.4)')  'Abs diff (no sigma vs FD):', abs(grad_zero(3,2) - fd)
   end block
end program repro_gcp_grad_contract
