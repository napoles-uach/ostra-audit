program repro_dqdr_contract
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use multicharge_model, only : mchrg_model_type
   use multicharge_param, only : new_eeq2019_model
   use multicharge_charge, only : get_charges
   implicit none
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(error_type), allocatable :: err
   real(wp) :: q(2), qp(2), qm(2), dqdr_full(3,2,2), dqdL_full(3,3,2), dqdr_noL(3,2,2), h, fd

   mol%nat = 2; mol%nid = 2; mol%id = [1, 2]; mol%num = [8, 1]; mol%periodic = .false.; mol%charge = 0.0_wp
   allocate(mol%xyz(3, 2))
   mol%xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   mol%xyz(:, 2) = [0.0_wp, 0.0_wp, 1.8_wp]
   call new_eeq2019_model(mol, model, err)
   if (allocated(err)) error stop err%message

   dqdr_full = 0.0_wp; dqdL_full = 0.0_wp
   call get_charges(model, mol, err, q, dqdr=dqdr_full, dqdL=dqdL_full)

   ! Non-periodic call: request dqdr (intent(out)), omit dqdL!
   dqdr_noL = -999.0_wp
   call get_charges(model, mol, err, q, dqdr=dqdr_noL)

   h = 1.0e-5_wp
   mol%xyz(3, 2) = 1.8_wp + h; call get_charges(model, mol, err, qp)
   mol%xyz(3, 2) = 1.8_wp - h; call get_charges(model, mol, err, qm)
   fd = (qp(1) - qm(1)) / (2.0_wp * h)

   write(*,'(A,2ES16.8)') 'Charges Q(O,H):         ', q(1), q(2)
   write(*,'(A,ES18.10)') 'dq_O/dz_H (FD ref):     ', fd
   write(*,'(A,ES18.10)') 'dq_O/dz_H (with dqdL):  ', dqdr_full(3, 2, 1)
   write(*,'(A,ES18.10)') 'dq_O/dz_H (no dqdL):    ', dqdr_noL(3, 2, 1)
   write(*,'(A,ES12.4)')  'Abs diff (no dqdL vs FD):', abs(dqdr_noL(3, 2, 1) - fd)
end program repro_dqdr_contract
