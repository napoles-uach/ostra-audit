program repro_maxelem
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use multicharge_model, only : mchrg_model_type
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2025_model
   use multicharge_charge, only : get_charges
   implicit none
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(error_type), allocatable :: err
   real(wp) :: q(2)
   integer :: ztest(3) = [103, 112, 118], iz
   character(len=2) :: sym(3) = ["Lr", "Cn", "Og"]

   mol%nat = 2; mol%nid = 2; mol%id = [1, 2]; mol%periodic = .false.; mol%charge = 0.0_wp
   allocate(mol%xyz(3, 2))
   mol%xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   mol%xyz(:, 2) = [0.0_wp, 0.0_wp, 4.0_wp]

   do iz = 1, 3
      mol%num = [ztest(iz), 8]
      if (allocated(model)) deallocate(model)
      if (allocated(err)) deallocate(err)
      call new_eeq2019_model(mol, model, err)
      if (allocated(err)) then
         write(*,'(A,A2,A,I3,A,A)') 'EEQ2019 ', sym(iz), ' (Z=', ztest(iz), '): ERROR -> ', err%message
      else
         call get_charges(model, mol, err, q)
         write(*,'(A,A2,A,I3,A,4F7.2,A,2F10.5)') 'EEQ2019 ', sym(iz), ' (Z=', ztest(iz), &
            '): NO ERROR! (chi,rad,eta,kcn)=', model%chi(1), model%rad(1), model%eta(1), model%kcnchi(1), &
            '  Q=', q(1), q(2)
      end if
   end do
end program repro_maxelem
