program t
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, new
   use multicharge, only : mchrg_model_type, new_eeq2019_model, get_charges
   implicit none
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(error_type), allocatable :: error
   real(wp), allocatable :: q(:), dqdr(:,:,:), dqdL(:,:,:), dqdr2(:,:,:)
   real(wp) :: xyz(3,3)
   ! OSTRA-026: water, dqdr requested alone vs together with dqdL
   xyz = reshape([0.0_wp,0.0_wp,0.221_wp, 0.0_wp,1.430_wp,-0.886_wp, 0.0_wp,-1.430_wp,-0.886_wp],[3,3])
   call new(mol, [8,1,1], xyz)
   call new_eeq2019_model(mol, model, error)
   allocate(q(3), dqdr(3,3,3), dqdL(3,3,3), dqdr2(3,3,3))
   call get_charges(model, mol, error, q, dqdr=dqdr, dqdL=dqdL)
   dqdr2 = 999.0_wp
   call get_charges(model, mol, error, q, dqdr=dqdr2)
   print '(a,es10.3,a,es10.3)', "026 max|dqdr| (with dqdL) = ", maxval(abs(dqdr)), &
      & "   max|dqdr(alone) - dqdr(with dqdL)| = ", maxval(abs(dqdr2 - dqdr))
   ! OSTRA-021: element beyond the EEQ2019 parameter range
   call run(80); call run(103); call run(112); call run(118)
contains
   subroutine run(zz)
      integer, intent(in) :: zz
      type(structure_type) :: m2
      real(wp) :: x2(3,2), q2(2)
      class(mchrg_model_type), allocatable :: md
      type(error_type), allocatable :: err
      x2 = reshape([0.0_wp,0.0_wp,0.0_wp, 0.0_wp,0.0_wp,3.8_wp],[3,2])
      call new(m2, [zz, 8], x2)
      call new_eeq2019_model(m2, md, err)
      if (allocated(err)) then
         print '(a,i4,2a)', "021 Z=", zz, "  error: ", err%message; return
      end if
      call get_charges(md, m2, err, q2)
      if (allocated(err)) then
         print '(a,i4,2a)', "021 Z=", zz, "  error: ", err%message; return
      end if
      print '(a,i4,a,2f9.4)', "021 Z=", zz, "  q(Z), q(O) = ", q2
   end subroutine run
end program t
