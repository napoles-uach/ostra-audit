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
