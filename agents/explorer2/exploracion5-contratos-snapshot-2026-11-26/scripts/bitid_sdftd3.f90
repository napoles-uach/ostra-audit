program bitid_sdftd3
   use, intrinsic :: iso_fortran_env, only : wp => real64, int64
   use mctc_env, only : error_type
   use mctc_io, only : structure_type, new
   use dftd3
   implicit none
   type(structure_type) :: mol
   type(d3_model) :: d3
   type(error_type), allocatable :: err
   integer, parameter :: nat = 3
   real(wp) :: xyz(3, nat), cn(nat), c6(nat, nat), dc6dcn(nat, nat)
   real(wp) :: d2c6dcn2(nat, nat), d2c6dcnij(nat, nat)
   real(wp), allocatable :: gw(:, :), gwdcn(:, :), gwd2cn(:, :)
   integer :: i, j, mref

   xyz(:, 1) = [0.00_wp, 0.00_wp, 0.12_wp]
   xyz(:, 2) = [0.00_wp, 1.43_wp, -0.95_wp]
   xyz(:, 3) = [0.15_wp, -1.40_wp, -0.98_wp]
   call new(mol, [8, 1, 1], xyz)
   call new_d3_model(d3, mol)
   cn = [1.8_wp, 0.9_wp, 0.9_wp]
   mref = maxval(d3%ref)
   allocate(gw(mref, nat), gwdcn(mref, nat), gwd2cn(mref, nat))
   call d3%weight_references(mol, cn, gw, gwdcn, gwd2cn)
   call d3%get_atomic_c6(mol, gw, gwdcn, c6, dc6dcn, gwd2cn, d2c6dcn2, d2c6dcnij)
   do j = 1, nat; do i = 1, mref
      print '(a,2(i0,1x),3(1x,z16.16))', "GW ", i, j, &
         & transfer(gw(i,j), 0_int64), transfer(gwdcn(i,j), 0_int64), transfer(gwd2cn(i,j), 0_int64)
   end do; end do
   do j = 1, nat; do i = 1, nat
      print '(a,2(i0,1x),4(1x,z16.16))', "C6 ", i, j, &
         & transfer(c6(i,j), 0_int64), transfer(dc6dcn(i,j), 0_int64), &
         & transfer(d2c6dcn2(i,j), 0_int64), transfer(d2c6dcnij(i,j), 0_int64)
   end do; end do
end program bitid_sdftd3
