program bitid_dftd4
   use, intrinsic :: iso_fortran_env, only : wp => real64, int64
   use mctc_env, only : error_type
   use mctc_io, only : structure_type, new
   use dftd4
   implicit none
   type(structure_type) :: mol
   type(d4_model) :: d4
   type(d4s_model) :: d4s
   class(damping_param), allocatable :: param
   type(realspace_cutoff) :: cutoff
   type(error_type), allocatable :: err
   integer, parameter :: nat = 3
   real(wp) :: xyz(3, nat), lat(3, 3), e0, e1, g(3, nat), s(3, 3)
   real(wp) :: cn(nat), q(nat), c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat)
   real(wp) :: alpha(nat), dadcn(nat), dadq(nat)
   real(wp), allocatable :: gw(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   integer :: i, j, k, mref

   xyz(:, 1) = [0.00_wp, 0.00_wp, 0.12_wp]
   xyz(:, 2) = [0.00_wp, 1.43_wp, -0.95_wp]
   xyz(:, 3) = [0.15_wp, -1.40_wp, -0.98_wp]
   lat = reshape([12.0_wp, 0.0_wp, 0.0_wp, 0.4_wp, 12.0_wp, 0.0_wp, 0.0_wp, 0.2_wp, 12.0_wp], [3, 3])
   call new(mol, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .true.])
   call new_d4_model(err, d4, mol)
   call new_d4s_model(err, d4s, mol)
   call get_rational_damping("pbe", param, s9=1.0_wp)
   cutoff = realspace_cutoff(disp2=40.0_wp, disp3=30.0_wp, cn=30.0_wp)
   call get_dispersion(mol, d4, param, cutoff, e0, error=err)
   call get_dispersion(mol, d4, param, cutoff, e1, g, s, error=err)
   print '(a,2(1x,z16.16))', "E ", transfer(e0, 0_int64), transfer(e1, 0_int64)
   do j = 1, nat; do i = 1, 3
      print '(a,2(i0,1x),z16.16)', "G ", i, j, transfer(g(i, j), 0_int64)
   end do; end do
   do j = 1, 3; do i = 1, 3
      print '(a,2(i0,1x),z16.16)', "S ", i, j, transfer(s(i, j), 0_int64)
   end do; end do

   cn = [1.6_wp, 0.8_wp, 0.8_wp]; q = [-0.4_wp, 0.2_wp, 0.2_wp]
   mref = maxval(d4s%ref)
   allocate(gw(mref, nat, d4s%ncoup), gwdcn(mref, nat, d4s%ncoup), gwdq(mref, nat, d4s%ncoup))
   call d4s%weight_references(mol, cn, q, gw, gwdcn, gwdq)
   call d4s%get_atomic_c6(mol, gw, gwdcn, gwdq, c6, dc6dcn, dc6dq)
   call d4s%get_polarizabilities(mol, gw, gwdcn, gwdq, alpha, dadcn, dadq)
   do k = 1, d4s%ncoup; do j = 1, nat; do i = 1, mref
      print '(a,3(i0,1x),3(1x,z16.16))', "GW ", i, j, k, &
         & transfer(gw(i,j,k), 0_int64), transfer(gwdcn(i,j,k), 0_int64), transfer(gwdq(i,j,k), 0_int64)
   end do; end do; end do
   do j = 1, nat; do i = 1, nat
      print '(a,2(i0,1x),3(1x,z16.16))', "C6 ", i, j, &
         & transfer(c6(i,j), 0_int64), transfer(dc6dcn(i,j), 0_int64), transfer(dc6dq(i,j), 0_int64)
   end do; end do
   do i = 1, nat
      print '(a,i0,3(1x,z16.16))', "AL ", i, &
         & transfer(alpha(i), 0_int64), transfer(dadcn(i), 0_int64), transfer(dadq(i), 0_int64)
   end do
end program bitid_dftd4
