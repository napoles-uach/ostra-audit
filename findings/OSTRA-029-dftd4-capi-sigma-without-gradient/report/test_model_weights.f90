program p33
  use mctc_env, only: wp, error_type
  use mctc_io, only: structure_type, new
  use dftd4_model, only: d4_model, new_d4_model
  implicit none
  type(structure_type) :: mol
  type(d4_model) :: d4
  type(error_type), allocatable :: err
  real(wp) :: xyz(3,4), cn(4), q(4)
  real(wp), allocatable :: gw(:,:,:), gdcn(:,:,:), gdq(:,:,:), a(:,:,:), b(:,:,:)
  integer :: mref
  xyz = reshape([0.0_wp,0.0_wp,0.22_wp, 0.0_wp,1.77_wp,-0.51_wp, 1.53_wp,-0.885_wp,-0.51_wp, -1.53_wp,-0.885_wp,-0.51_wp],[3,4])
  call new(mol, [7,1,1,1], xyz)
  call new_d4_model(err, d4, mol)
  cn = [2.9_wp, 0.95_wp, 0.95_wp, 0.95_wp]; q = [-0.6_wp, 0.2_wp, 0.2_wp, 0.2_wp]
  mref = maxval(d4%ref)
  allocate(gw(mref,4,d4%ncoup)); allocate(gdcn, gdq, a, b, mold=gw)
  call d4%weight_references(mol, cn, q, gw, gdcn, gdq)          ! reference: both
  a = 999.0_wp; call d4%weight_references(mol, cn, q, gw, gwdcn=a)   ! CN derivative only
  b = 999.0_wp; call d4%weight_references(mol, cn, q, gw, gwdq=b)    ! q derivative only
  print '(a,es10.3,a,es10.3,a,es10.3,a,es10.3)', 'max|gwdcn|=',maxval(abs(gdcn)),'  gwdcn-only err=',maxval(abs(a-gdcn)), &
     & '  max|gwdq|=',maxval(abs(gdq)),'  gwdq-only err=',maxval(abs(b-gdq))
end program
