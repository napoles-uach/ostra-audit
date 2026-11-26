program bitid_mctc
   use, intrinsic :: iso_fortran_env, only : wp => real64, int64
   use mctc_env, only : error_type
   use mctc_io, only : structure_type, new
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   implicit none
   type(structure_type) :: mol
   class(ncoord_type), allocatable :: nc
   type(error_type), allocatable :: err
   integer, parameter :: nat = 3
   real(wp) :: xyz(3, nat), lat(3, 3), cn0(nat), cn1(nat), dcndr(3, nat, nat), dcndL(3, 3, nat)
   integer :: i, j, k

   xyz(:, 1) = [0.00_wp, 0.00_wp, 0.12_wp]
   xyz(:, 2) = [0.00_wp, 1.43_wp, -0.95_wp]
   xyz(:, 3) = [0.15_wp, -1.40_wp, -0.98_wp]
   lat = reshape([10.0_wp, 0.0_wp, 0.0_wp, 0.5_wp, 10.0_wp, 0.0_wp, 0.0_wp, 0.3_wp, 10.0_wp], [3, 3])
   call new(mol, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .true.])
   call new_ncoord(nc, mol, cn_count%erf, err, cut=6.0_wp)
   call nc%get_cn(mol, cn0)
   call nc%get_cn(mol, cn1, dcndr, dcndL)
   do i = 1, nat
      print '(a,i0,1x,z16.16,1x,z16.16)', "CN ", i, transfer(cn0(i), 0_int64), transfer(cn1(i), 0_int64)
   end do
   do k = 1, nat; do j = 1, nat; do i = 1, 3
      print '(a,3(i0,1x),z16.16)', "DR ", i, j, k, transfer(dcndr(i, j, k), 0_int64)
   end do; end do; end do
   do k = 1, nat; do j = 1, 3; do i = 1, 3
      print '(a,3(i0,1x),z16.16)', "DL ", i, j, k, transfer(dcndL(i, j, k), 0_int64)
   end do; end do; end do
end program bitid_mctc
