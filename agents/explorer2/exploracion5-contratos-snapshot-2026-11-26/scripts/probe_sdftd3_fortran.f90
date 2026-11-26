program probe_sdftd3_fortran
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, new
   use dftd3, only : d3_model, new_d3_model, damping_param, d3_param, &
      & rational_damping_param, new_rational_damping, get_rational_damping, &
      & zero_damping_param, new_zero_damping, get_zero_damping, &
      & mzero_damping_param, new_mzero_damping, get_mzero_damping, &
      & optimizedpower_damping_param, new_optimizedpower_damping, get_optimizedpower_damping, &
      & cso_damping_param, new_cso_damping, get_cso_damping, &
      & z_damping_param, new_z_damping, get_z_damping, &
      & realspace_cutoff, get_dispersion, get_lattice_points, get_coordination_number
   implicit none

   real(wp), parameter :: SENTINEL = 999.0_wp
   character(len=64) :: mode
   type(structure_type) :: mol
   type(d3_model) :: d3
   type(d3_param) :: inp
   type(rational_damping_param) :: rat
   type(zero_damping_param) :: zer
   type(mzero_damping_param) :: mzer
   type(optimizedpower_damping_param) :: op
   type(cso_damping_param) :: cso
   type(z_damping_param) :: zd
   type(realspace_cutoff) :: cutoff
   type(error_type), allocatable :: err
   integer :: num(3), nat, mref, mask
   real(wp) :: xyz(3, 3), lat(3, 3)
   logical :: pbc(3)
   real(wp) :: e_ref, g_ref(3, 3), s_ref(3, 3), h_ref(9, 9)
   real(wp) :: e, g(3, 3), s(3, 3), h(9, 9)
   real(wp), allocatable :: lattr(:, :), cn(:), dcndr_ref(:, :, :), dcndL_ref(:, :, :)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: gw_ref(:, :), gwdcn_ref(:, :), gwd2cn_ref(:, :)
   real(wp), allocatable :: gw(:, :), gwdcn(:, :), gwd2cn(:, :)
   real(wp), allocatable :: c6_ref(:, :), dc6dcn_ref(:, :), d2c6_ref(:, :), d2mix_ref(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), d2c6(:, :), d2mix(:, :)
   real(wp) :: en_ref(3), dEdcn_ref(3), gd_ref(3, 3), sd_ref(3, 3)
   real(wp) :: en(3), dEdcn(3), gd(3, 3), sd(3, 3)

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
   call new_d3_model(d3, mol)
   call get_rational_damping(inp, "pbe", err, s9=1.0_wp); call new_rational_damping(rat, inp)
   call get_zero_damping(inp, "pbe", err, s9=1.0_wp); call new_zero_damping(zer, inp)
   call get_mzero_damping(inp, "pbe", err, s9=1.0_wp); call new_mzero_damping(mzer, inp)
   call get_optimizedpower_damping(inp, "pbe", err, s9=1.0_wp); call new_optimizedpower_damping(op, inp)
   call get_cso_damping(inp, "pbe", err, s9=1.0_wp); call new_cso_damping(cso, inp)
   call get_z_damping(inp, "pbe", err, s9=1.0_wp); call new_z_damping(zd, inp)
   cutoff = realspace_cutoff()

   allocate(cn(nat), dcndr_ref(3, nat, nat), dcndL_ref(3, 3, nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, d3%rcov, cn, dcndr_ref, dcndL_ref)

   mref = maxval(d3%ref)
   allocate(gw_ref(mref, nat), gwdcn_ref(mref, nat), gwd2cn_ref(mref, nat))
   allocate(gw(mref, nat), gwdcn(mref, nat), gwd2cn(mref, nat))
   call d3%weight_references(mol, cn, gw_ref, gwdcn_ref, gwd2cn_ref)

   allocate(c6_ref(nat, nat), dc6dcn_ref(nat, nat), d2c6_ref(nat, nat), d2mix_ref(nat, nat))
   allocate(c6(nat, nat), dc6dcn(nat, nat), d2c6(nat, nat), d2mix(nat, nat))
   call d3%get_atomic_c6(mol, gw_ref, gwdcn_ref, c6_ref, dc6dcn_ref, gwd2cn_ref, d2c6_ref, d2mix_ref)

   if (trim(mode) == "c6_second_no_first") then
      print '(a)', "Calling d3%get_atomic_c6(gwd2cn, d2c6dcn2, d2c6dcnij) WITHOUT gwdcn/dc6dcn..."
      c6 = SENTINEL; d2c6 = SENTINEL; d2mix = SENTINEL
      call d3%get_atomic_c6(mol, gw_ref, c6=c6, gwd2cn=gwd2cn_ref, d2c6dcn2=d2c6, d2c6dcnij=d2mix)
      print '(a,es12.5)', "RETURNED: max|d2c6-ref|=", maxval(abs(d2c6 - d2c6_ref))
      stop
   end if

   print '(a)', "=== 1. simple-dftd3 CONTROL: get_dispersion (all 2^3 grad/sigma/hessian masks) ==="
   e_ref = SENTINEL; g_ref = SENTINEL; s_ref = SENTINEL; h_ref = SENTINEL
   call get_dispersion(err, mol, d3, rat, cutoff, e_ref, g_ref, s_ref, h_ref)
   if (allocated(err)) error stop err%message
   do mask = 0, 7
      e = SENTINEL; g = SENTINEL; s = SENTINEL; h = SENTINEL
      select case(mask)
      case(0); call get_dispersion(err, mol, d3, rat, cutoff, e)
      case(1); call get_dispersion(err, mol, d3, rat, cutoff, e, gradient=g)
      case(2); call get_dispersion(err, mol, d3, rat, cutoff, e, sigma=s)
      case(3); call get_dispersion(err, mol, d3, rat, cutoff, e, gradient=g, sigma=s)
      case(4); call get_dispersion(err, mol, d3, rat, cutoff, e, hessian=h)
      case(5); call get_dispersion(err, mol, d3, rat, cutoff, e, gradient=g, hessian=h)
      case(6); call get_dispersion(err, mol, d3, rat, cutoff, e, sigma=s, hessian=h)
      case(7); call get_dispersion(err, mol, d3, rat, cutoff, e, gradient=g, sigma=s, hessian=h)
      end select
      if (allocated(err)) error stop err%message
      print '(a,i1,a,es10.3,a,es10.3,a,es10.3,a,es10.3)', "  mask=", mask, &
         & " |dE|=", abs(e - e_ref), &
         & " |dg|=", merge(maxval(abs(g - g_ref)), 0.0_wp, btest(mask, 0)), &
         & " |ds|=", merge(maxval(abs(s - s_ref)), 0.0_wp, btest(mask, 1)), &
         & " |dh|=", merge(maxval(abs(h - h_ref)), 0.0_wp, btest(mask, 2))
   end do

   print '(a)', ""
   print '(a)', "=== 2. simple-dftd3 d3%weight_references (2^2 combinations of gwdcn, gwd2cn) ==="
   print '(a,es12.5,a,es12.5)', "  ref max|gwdcn|=", maxval(abs(gwdcn_ref)), " max|gwd2cn|=", maxval(abs(gwd2cn_ref))
   gw = SENTINEL; gwdcn = SENTINEL
   call d3%weight_references(mol, cn, gw, gwdcn=gwdcn)
   print '(a,es12.5)', "  (T,F) gwdcn only:  max|err|=", maxval(abs(gwdcn - gwdcn_ref))
   gw = SENTINEL; gwd2cn = SENTINEL
   call d3%weight_references(mol, cn, gw, gwd2cn=gwd2cn)
   print '(a,es12.5,a,es12.5)', "  (F,T) gwd2cn only: max|err|=", maxval(abs(gwd2cn - gwd2cn_ref)), &
      & " max|gwd2cn|=", maxval(abs(gwd2cn))

   print '(a)', ""
   print '(a)', "=== 3. simple-dftd3 d3%get_atomic_c6 (partial optional combinations) ==="
   c6 = SENTINEL; dc6dcn = SENTINEL; d2c6 = SENTINEL
   call d3%get_atomic_c6(mol, gw_ref, gwdcn=gwdcn_ref, c6=c6, dc6dcn=dc6dcn, &
      & gwd2cn=gwd2cn_ref, d2c6dcn2=d2c6)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  (gwdcn,dc6dcn,gwd2cn,d2c6dcn2; no d2c6dcnij): |dc6|=", &
      & maxval(abs(dc6dcn - dc6dcn_ref)), " |d2c6|=", maxval(abs(d2c6 - d2c6_ref)), " min(d2c6)=", minval(d2c6)
   c6 = SENTINEL; dc6dcn = SENTINEL; d2mix = SENTINEL
   call d3%get_atomic_c6(mol, gw_ref, gwdcn=gwdcn_ref, c6=c6, dc6dcn=dc6dcn, &
      & gwd2cn=gwd2cn_ref, d2c6dcnij=d2mix)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  (gwdcn,dc6dcn,gwd2cn,d2c6dcnij; no d2c6dcn2): |dc6|=", &
      & maxval(abs(dc6dcn - dc6dcn_ref)), " |d2mix|=", maxval(abs(d2mix - d2mix_ref)), " min(d2mix)=", minval(d2mix)

   print '(a)', ""
   print '(a)', "=== 4. simple-dftd3 ncoord::get_coordination_number (inherits mctc-lib) ==="
   dcndr = SENTINEL
   call get_coordination_number(mol, lattr, cutoff%cn, d3%rcov, cn, dcndr=dcndr)
   print '(a,es12.5,a,es12.5)', "  (T,F) dcndr only: max|err|=", maxval(abs(dcndr - dcndr_ref)), " min=", minval(dcndr)
   dcndL = SENTINEL
   call get_coordination_number(mol, lattr, cutoff%cn, d3%rcov, cn, dcndL=dcndL)
   print '(a,es12.5,a,es12.5)', "  (F,T) dcndL only: max|err|=", maxval(abs(dcndL - dcndL_ref)), " min=", minval(dcndL)

   print '(a)', ""
   print '(a)', "=== 5. simple-dftd3 damping%get_dispersion2/3 (all 6 damping types + ATM) ==="
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call check_damp("rational      ", rat)
   call check_damp("zero          ", zer)
   call check_damp("mzero         ", mzer)
   call check_damp("optimizedpower", op)
   call check_damp("cso           ", cso)
   call check_damp("z             ", zd)

contains

   subroutine check_damp(label, p)
      character(len=*), intent(in) :: label
      class(damping_param), intent(in) :: p
      en_ref = 0.0_wp; dEdcn_ref = 0.0_wp; gd_ref = 0.0_wp; sd_ref = 0.0_wp
      call p%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, d3%rvdw, d3%r4r2, &
         & c6_ref, dc6dcn_ref, en_ref, dEdcn_ref, gd_ref, sd_ref)
      en = 0.0_wp; dEdcn = SENTINEL; gd = SENTINEL
      call p%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, d3%rvdw, d3%r4r2, &
         & c6_ref, dc6dcn_ref, en, dEdcn=dEdcn, gradient=gd)
      print '(a,a,a,es10.3,a,es10.3,a,es10.3)', "  ", trim(label), " disp2(no sigma): ref_g=", &
         & maxval(abs(gd_ref)), " |gd-ref|=", maxval(abs(gd - gd_ref)), " min(gd)=", minval(gd)
      if (trim(label) == "rational") then
         en_ref = 0.0_wp; dEdcn_ref = 0.0_wp; gd_ref = 0.0_wp; sd_ref = 0.0_wp
         call p%get_dispersion3(mol, lattr, cutoff%disp3, cutoff%width3, d3%rvdw, d3%r4r2, &
            & c6_ref, dc6dcn_ref, en_ref, dEdcn_ref, gd_ref, sd_ref)
         en = 0.0_wp; dEdcn = SENTINEL; gd = SENTINEL
         call p%get_dispersion3(mol, lattr, cutoff%disp3, cutoff%width3, d3%rvdw, d3%r4r2, &
            & c6_ref, dc6dcn_ref, en, dEdcn=dEdcn, gradient=gd)
         print '(a,es10.3,a,es10.3,a,es10.3)', "  atm          disp3(no sigma): ref_g=", &
            & maxval(abs(gd_ref)), " |gd-ref|=", maxval(abs(gd - gd_ref)), " min(gd)=", minval(gd)
      end if
   end subroutine check_damp
end program probe_sdftd3_fortran
