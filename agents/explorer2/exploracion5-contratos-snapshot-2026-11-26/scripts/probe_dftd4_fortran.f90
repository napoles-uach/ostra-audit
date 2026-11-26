program probe_dftd4_fortran
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, new
   use dftd4, only : d4_model, new_d4_model, d4s_model, new_d4s_model, &
      & damping_param, get_rational_damping, realspace_cutoff, &
      & get_dispersion, get_lattice_points, get_coordination_number
   use dftd4_disp, only : get_dispersion2
   implicit none

   real(wp), parameter :: SENTINEL = 999.0_wp
   character(len=64) :: mode
   type(structure_type) :: mol, mol_p
   type(d4_model) :: d4
   type(d4s_model) :: d4s
   class(damping_param), allocatable :: param
   type(realspace_cutoff) :: cutoff
   type(error_type), allocatable :: err
   integer :: num(3), nat, mref, iat, idir, a, b, ih
   real(wp) :: xyz(3, 3), lat(3, 3), xyz_tmp(3, 3), lat_tmp(3, 3), eps_mat(3, 3)
   logical :: pbc(3)
   real(wp) :: e_ref, g_ref(3, 3), s_ref(3, 3), e, g(3, 3), s(3, 3)
   real(wp) :: ep1, em1, ep2, em2, fd_val, max_fd_g, max_fd_s, h
   real(wp) :: hs(4) = [1.0e-2_wp, 1.0e-3_wp, 1.0e-4_wp, 1.0e-5_wp]
   real(wp), allocatable :: lattr(:, :), cn(:), q(:), dcndr_ref(:, :, :), dcndL_ref(:, :, :)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: gw_ref(:, :, :), gwdcn_ref(:, :, :), gwdq_ref(:, :, :)
   real(wp), allocatable :: gw(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), allocatable :: c6_ref(:, :), dc6dcn_ref(:, :), dc6dq_ref(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: a_ref(:), dadcn_ref(:), dadq_ref(:)
   real(wp), allocatable :: alpha(:), dadcn(:), dadq(:)
   real(wp) :: en_ref(3), dEdcn_ref(3), dEdq_ref(3), gd2_ref(3, 3), sd2_ref(3, 3)
   real(wp) :: en(3), dEdcn(3), dEdq(3), gd2(3, 3), sd2(3, 3)

   mode = "all_safe"
   if (command_argument_count() >= 1) call get_command_argument(1, mode)

   nat = 3
   num = [8, 1, 1]
   xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   xyz(:, 2) = [1.80_wp, 0.35_wp, 0.15_wp]
   xyz(:, 3) = [-0.45_wp, 1.75_wp, -0.25_wp]
   lat(:, 1) = [8.0_wp, 0.2_wp, 0.1_wp]
   lat(:, 2) = [0.3_wp, 8.2_wp, 0.15_wp]
   lat(:, 3) = [0.1_wp, 0.25_wp, 8.5_wp]
   pbc = [.true., .true., .true.]
   call new(mol, num, xyz, lattice=lat, periodic=pbc)
   call new_d4_model(err, d4, mol)
   if (allocated(err)) error stop err%message
   call new_d4s_model(err, d4s, mol)
   if (allocated(err)) error stop err%message
   call get_rational_damping("pbe", param, s9=1.0_wp)
   cutoff = realspace_cutoff()

   ! Crash-prone cases isolated by mode
   if (trim(mode) == "disp_grad_only") then
      e = SENTINEL; g = SENTINEL
      print '(a)', "Calling dftd4 get_dispersion(gradient only)..."
      call get_dispersion(mol, d4, param, cutoff, e, gradient=g, error=err)
      print '(a,es12.5,a,es12.5)', "RETURNED: E=", e, " max|g|=", maxval(abs(g))
      stop
   else if (trim(mode) == "disp_sigma_only") then
      e = SENTINEL; s = SENTINEL
      print '(a)', "Calling dftd4 get_dispersion(sigma only)..."
      call get_dispersion(mol, d4, param, cutoff, e, sigma=s, error=err)
      print '(a,es12.5,a,es12.5)', "RETURNED: E=", e, " max|s|=", maxval(abs(s))
      stop
   else if (trim(mode) == "disp2_grad_only") then
      e = SENTINEL; g = SENTINEL
      print '(a)', "Calling dftd4 get_dispersion2(gradient only)..."
      call get_dispersion2(mol, d4, param, cutoff, e, gradient=g, error=err)
      print '(a,es12.5,a,es12.5)', "RETURNED: E=", e, " max|g|=", maxval(abs(g))
      stop
   else if (trim(mode) == "disp2_sigma_only") then
      e = SENTINEL; s = SENTINEL
      print '(a)', "Calling dftd4 get_dispersion2(sigma only)..."
      call get_dispersion2(mol, d4, param, cutoff, e, sigma=s, error=err)
      print '(a,es12.5,a,es12.5)', "RETURNED: E=", e, " max|s|=", maxval(abs(s))
      stop
   end if

   print '(a)', "=== 1. dftd4 get_dispersion reference + 4pt FD physical oracle ==="
   e_ref = SENTINEL; g_ref = SENTINEL; s_ref = SENTINEL
   call get_dispersion(mol, d4, param, cutoff, e_ref, g_ref, s_ref, error=err)
   if (allocated(err)) error stop err%message
   print '(a,es16.8,a,es12.5,a,es12.5)', "REF E=", e_ref, " max|g|=", maxval(abs(g_ref)), &
      & " max|s|=", maxval(abs(s_ref))
   e = SENTINEL
   call get_dispersion(mol, d4, param, cutoff, e, error=err)
   print '(a,es12.5)', "  (F,F) |E-ref|=", abs(e - e_ref)

   do ih = 1, 4
      h = hs(ih)
      max_fd_g = 0.0_wp
      do iat = 1, nat
         do idir = 1, 3
            xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) + h
            call new(mol_p, num, xyz_tmp, lattice=lat, periodic=pbc)
            call get_dispersion(mol_p, d4, param, cutoff, ep1, error=err)
            xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) - h
            call new(mol_p, num, xyz_tmp, lattice=lat, periodic=pbc)
            call get_dispersion(mol_p, d4, param, cutoff, em1, error=err)
            xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) + 2.0_wp*h
            call new(mol_p, num, xyz_tmp, lattice=lat, periodic=pbc)
            call get_dispersion(mol_p, d4, param, cutoff, ep2, error=err)
            xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) - 2.0_wp*h
            call new(mol_p, num, xyz_tmp, lattice=lat, periodic=pbc)
            call get_dispersion(mol_p, d4, param, cutoff, em2, error=err)
            fd_val = (-ep2 + 8.0_wp*ep1 - 8.0_wp*em1 + em2) / (12.0_wp*h)
            max_fd_g = max(max_fd_g, abs(fd_val - g_ref(idir, iat)))
         end do
      end do
      max_fd_s = 0.0_wp
      do a = 1, 3
         do b = 1, 3
            eps_mat = 0.0_wp; eps_mat(a, b) = h
            call new(mol_p, num, xyz + matmul(eps_mat, xyz), lattice=lat + matmul(eps_mat, lat), periodic=pbc)
            call get_dispersion(mol_p, d4, param, cutoff, ep1, error=err)
            eps_mat = 0.0_wp; eps_mat(a, b) = -h
            call new(mol_p, num, xyz + matmul(eps_mat, xyz), lattice=lat + matmul(eps_mat, lat), periodic=pbc)
            call get_dispersion(mol_p, d4, param, cutoff, em1, error=err)
            eps_mat = 0.0_wp; eps_mat(a, b) = 2.0_wp*h
            call new(mol_p, num, xyz + matmul(eps_mat, xyz), lattice=lat + matmul(eps_mat, lat), periodic=pbc)
            call get_dispersion(mol_p, d4, param, cutoff, ep2, error=err)
            eps_mat = 0.0_wp; eps_mat(a, b) = -2.0_wp*h
            call new(mol_p, num, xyz + matmul(eps_mat, xyz), lattice=lat + matmul(eps_mat, lat), periodic=pbc)
            call get_dispersion(mol_p, d4, param, cutoff, em2, error=err)
            fd_val = (-ep2 + 8.0_wp*ep1 - 8.0_wp*em1 + em2) / (12.0_wp*h)
            max_fd_s = max(max_fd_s, abs(fd_val - s_ref(a, b)))
         end do
      end do
      print '(a,es8.1,a,es11.3,a,es11.3)', "  4pt-FD h=", h, " max|g-FD|=", max_fd_g, " max|s-FD|=", max_fd_s
   end do

   print '(a)', ""
   print '(a)', "=== 2. dftd4_ncoord::get_coordination_number (2^2 combinations) ==="
   allocate(cn(nat), dcndr_ref(3, nat, nat), dcndL_ref(3, 3, nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, d4%rcov, d4%en, cn, dcndr_ref, dcndL_ref, err)
   dcndr = SENTINEL
   call get_coordination_number(mol, lattr, cutoff%cn, d4%rcov, d4%en, cn, dcndr=dcndr, error=err)
   print '(a,es12.5,a,es12.5)', "  (T,F) dcndr only: max|err|=", maxval(abs(dcndr - dcndr_ref)), " min=", minval(dcndr)
   dcndL = SENTINEL
   call get_coordination_number(mol, lattr, cutoff%cn, d4%rcov, d4%en, cn, dcndL=dcndL, error=err)
   print '(a,es12.5,a,es12.5)', "  (F,T) dcndL only: max|err|=", maxval(abs(dcndL - dcndL_ref)), " min=", minval(dcndL)

   print '(a)', ""
   print '(a)', "=== 3. dftd4 d4_model: weight_references, get_atomic_c6, get_polarizabilities ==="
   allocate(q(nat)); q = [-0.4_wp, 0.2_wp, 0.2_wp]
   mref = maxval(d4%ref)
   allocate(gw_ref(mref, nat, d4%ncoup), gwdcn_ref(mref, nat, d4%ncoup), gwdq_ref(mref, nat, d4%ncoup))
   allocate(gw(mref, nat, d4%ncoup), gwdcn(mref, nat, d4%ncoup), gwdq(mref, nat, d4%ncoup))
   call d4%weight_references(mol, cn, q, gw_ref, gwdcn_ref, gwdq_ref)
   gw = SENTINEL; gwdcn = SENTINEL
   call d4%weight_references(mol, cn, q, gw, gwdcn=gwdcn)
   print '(a,es12.5,a,es12.5)', "  d4%weight_references (gwdcn only): max|err|=", &
      & maxval(abs(gwdcn - gwdcn_ref)), " min=", minval(gwdcn)
   gw = SENTINEL; gwdq = SENTINEL
   call d4%weight_references(mol, cn, q, gw, gwdq=gwdq)
   print '(a,es12.5,a,es12.5)', "  d4%weight_references (gwdq only):  max|err|=", &
      & maxval(abs(gwdq - gwdq_ref)), " min=", minval(gwdq)

   allocate(c6_ref(nat, nat), dc6dcn_ref(nat, nat), dc6dq_ref(nat, nat))
   allocate(c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat))
   call d4%get_atomic_c6(mol, gw_ref, gwdcn_ref, gwdq_ref, c6_ref, dc6dcn_ref, dc6dq_ref)
   c6 = SENTINEL; dc6dcn = SENTINEL
   call d4%get_atomic_c6(mol, gw_ref, gwdcn=gwdcn_ref, c6=c6, dc6dcn=dc6dcn)
   print '(a,es12.5,a,es12.5)', "  d4%get_atomic_c6 (dc6dcn only):    max|err|=", &
      & maxval(abs(dc6dcn - dc6dcn_ref)), " min=", minval(dc6dcn)
   c6 = SENTINEL; dc6dq = SENTINEL
   call d4%get_atomic_c6(mol, gw_ref, gwdq=gwdq_ref, c6=c6, dc6dq=dc6dq)
   print '(a,es12.5,a,es12.5)', "  d4%get_atomic_c6 (dc6dq only):     max|err|=", &
      & maxval(abs(dc6dq - dc6dq_ref)), " min=", minval(dc6dq)

   allocate(a_ref(nat), dadcn_ref(nat), dadq_ref(nat), alpha(nat), dadcn(nat), dadq(nat))
   call d4%get_polarizabilities(mol, gw_ref, gwdcn_ref, gwdq_ref, a_ref, dadcn_ref, dadq_ref)
   alpha = SENTINEL; dadcn = SENTINEL
   call d4%get_polarizabilities(mol, gw_ref, gwdcn=gwdcn_ref, alpha=alpha, dadcn=dadcn)
   print '(a,es12.5,a,es12.5)', "  d4%get_polarizabilities (dadcn):   max|err|=", &
      & maxval(abs(dadcn - dadcn_ref)), " min=", minval(dadcn)
   alpha = SENTINEL; dadq = SENTINEL
   call d4%get_polarizabilities(mol, gw_ref, gwdq=gwdq_ref, alpha=alpha, dadq=dadq)
   print '(a,es12.5,a,es12.5)', "  d4%get_polarizabilities (dadq):    max|err|=", &
      & maxval(abs(dadq - dadq_ref)), " min=", minval(dadq)

   ! Also check d4s_model (with its own d4s references!)
   deallocate(gw_ref, gwdcn_ref, gwdq_ref, gw, gwdcn, gwdq)
   mref = maxval(d4s%ref)
   allocate(gw_ref(mref, nat, d4s%ncoup), gwdcn_ref(mref, nat, d4s%ncoup), gwdq_ref(mref, nat, d4s%ncoup))
   allocate(gw(mref, nat, d4s%ncoup), gwdcn(mref, nat, d4s%ncoup), gwdq(mref, nat, d4s%ncoup))
   call d4s%weight_references(mol, cn, q, gw_ref, gwdcn_ref, gwdq_ref)
   call d4s%get_atomic_c6(mol, gw_ref, gwdcn_ref, gwdq_ref, c6_ref, dc6dcn_ref, dc6dq_ref)
   call d4s%get_polarizabilities(mol, gw_ref, gwdcn_ref, gwdq_ref, a_ref, dadcn_ref, dadq_ref)
   gwdcn = SENTINEL
   call d4s%weight_references(mol, cn, q, gw, gwdcn=gwdcn)
   print '(a,es12.5,a,es12.5)', "  d4s%weight_references (gwdcn only):max|err|=", &
      & maxval(abs(gwdcn - gwdcn_ref)), " min=", minval(gwdcn)
   gwdq = SENTINEL
   call d4s%weight_references(mol, cn, q, gw, gwdq=gwdq)
   print '(a,es12.5,a,es12.5)', "  d4s%weight_references (gwdq only): max|err|=", &
      & maxval(abs(gwdq - gwdq_ref)), " min=", minval(gwdq)
   dc6dcn = SENTINEL
   call d4s%get_atomic_c6(mol, gw_ref, gwdcn=gwdcn_ref, c6=c6, dc6dcn=dc6dcn)
   print '(a,es12.5,a,es12.5)', "  d4s%get_atomic_c6 (dc6dcn only):   max|err|=", &
      & maxval(abs(dc6dcn - dc6dcn_ref)), " min=", minval(dc6dcn)
   dc6dq = SENTINEL
   call d4s%get_atomic_c6(mol, gw_ref, gwdq=gwdq_ref, c6=c6, dc6dq=dc6dq)
   print '(a,es12.5,a,es12.5)', "  d4s%get_atomic_c6 (dc6dq only):    max|err|=", &
      & maxval(abs(dc6dq - dc6dq_ref)), " min=", minval(dc6dq)
   dadcn = SENTINEL
   call d4s%get_polarizabilities(mol, gw_ref, gwdcn=gwdcn_ref, alpha=alpha, dadcn=dadcn)
   print '(a,es12.5,a,es12.5)', "  d4s%get_polarizabilities (dadcn):  max|err|=", &
      & maxval(abs(dadcn - dadcn_ref)), " min=", minval(dadcn)
   dadq = SENTINEL
   call d4s%get_polarizabilities(mol, gw_ref, gwdq=gwdq_ref, alpha=alpha, dadq=dadq)
   print '(a,es12.5,a,es12.5)', "  d4s%get_polarizabilities (dadq):   max|err|=", &
      & maxval(abs(dadq - dadq_ref)), " min=", minval(dadq)

   print '(a)', ""
   print '(a)', "=== 4. dftd4 rational_damping_param: get_dispersion2 & get_dispersion3 ==="
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   en_ref = 0.0_wp; dEdcn_ref = 0.0_wp; dEdq_ref = 0.0_wp; gd2_ref = 0.0_wp; sd2_ref = 0.0_wp
   call param%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, d4%r4r2, c6_ref, &
      & dc6dcn_ref, dc6dq_ref, en_ref, dEdcn_ref, dEdq_ref, gd2_ref, sd2_ref)
   ! Omit sigma only
   en = 0.0_wp; dEdcn = SENTINEL; dEdq = SENTINEL; gd2 = SENTINEL
   call param%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, d4%r4r2, c6_ref, &
      & dc6dcn=dc6dcn_ref, dc6dq=dc6dq_ref, energy=en, dEdcn=dEdcn, dEdq=dEdq, gradient=gd2)
   print '(a,es12.5,a,es12.5)', "  disp2 (no sigma): max|gd2-ref|=", maxval(abs(gd2 - gd2_ref)), " min(gd2)=", minval(gd2)
   ! Omit gradient only
   en = 0.0_wp; dEdcn = SENTINEL; dEdq = SENTINEL; sd2 = SENTINEL
   call param%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, d4%r4r2, c6_ref, &
      & dc6dcn=dc6dcn_ref, dc6dq=dc6dq_ref, energy=en, dEdcn=dEdcn, dEdq=dEdq, sigma=sd2)
   print '(a,es12.5,a,es12.5)', "  disp2 (no grad):  max|sd2-ref|=", maxval(abs(sd2 - sd2_ref)), " min(sd2)=", minval(sd2)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   en_ref = 0.0_wp; dEdcn_ref = 0.0_wp; dEdq_ref = 0.0_wp; gd2_ref = 0.0_wp; sd2_ref = 0.0_wp
   call param%get_dispersion3(mol, lattr, cutoff%disp3, cutoff%width3, d4%r4r2, c6_ref, &
      & dc6dcn_ref, dc6dq_ref, en_ref, dEdcn_ref, dEdq_ref, gd2_ref, sd2_ref)
   en = 0.0_wp; dEdcn = SENTINEL; dEdq = SENTINEL; gd2 = SENTINEL
   call param%get_dispersion3(mol, lattr, cutoff%disp3, cutoff%width3, d4%r4r2, c6_ref, &
      & dc6dcn=dc6dcn_ref, dc6dq=dc6dq_ref, energy=en, dEdcn=dEdcn, dEdq=dEdq, gradient=gd2)
   print '(a,es12.5,a,es12.5)', "  disp3 (no sigma): max|gd3-ref|=", maxval(abs(gd2 - gd2_ref)), " min(gd3)=", minval(gd2)
end program probe_dftd4_fortran
