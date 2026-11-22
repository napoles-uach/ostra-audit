program repro_gcp_oob
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use dftd3_cutoff, only : realspace_cutoff
   use dftd3_gcp_param, only : gcp_param, get_gcp_param
   use dftd3_gcp, only : get_geometric_counterpoise
   implicit none
   type(structure_type) :: mol
   type(gcp_param) :: param
   type(realspace_cutoff) :: cutoff
   real(wp) :: energy, grad(3, 2), sigma(3, 3)
   character(len=16) :: arg

   call get_command_argument(1, arg)
   mol%nat = 2; mol%nid = 2; mol%id = [1, 2]; mol%periodic = .false.
   allocate(mol%xyz(3, 2))
   mol%xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   mol%xyz(:, 2) = [0.0_wp, 0.0_wp, 4.0_wp]
   cutoff = realspace_cutoff()
   energy = 0.0_wp; grad = 0.0_wp; sigma = 0.0_wp

   select case (trim(arg))
   case ("Fr")
      mol%num = [87, 8]
      call get_gcp_param(param, mol, method="r2scan3c")
      call get_geometric_counterpoise(mol, param, cutoff, energy, grad, sigma)
      write(*,'(A,2I4,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') 'Fr-O r2scan3c: zeff=', param%zeff, &
         ' emiss1=', param%emiss(1), ' slater1=', param%slater(1), ' E=', energy, ' gz=', grad(3,2)
   case ("Ra")
      mol%num = [88, 8]
      call get_gcp_param(param, mol, method="r2scan3c")
      call get_geometric_counterpoise(mol, param, cutoff, energy, grad, sigma)
      write(*,'(A,2I4,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') 'Ra-O r2scan3c: zeff=', param%zeff, &
         ' emiss1=', param%emiss(1), ' slater1=', param%slater(1), ' E=', energy, ' gz=', grad(3,2)
   case ("Ac")
      mol%num = [89, 8]
      call get_gcp_param(param, mol, method="r2scan3c")
      call get_geometric_counterpoise(mol, param, cutoff, energy, grad, sigma)
      write(*,'(A,2I4,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') 'Ac-O r2scan3c: zeff=', param%zeff, &
         ' emiss1=', param%emiss(1), ' slater1=', param%slater(1), ' E=', energy, ' gz=', grad(3,2)
   case ("Kr")
      mol%num = [36, 8]
      call get_gcp_param(param, mol, method="hf", basis="minis")
      call get_geometric_counterpoise(mol, param, cutoff, energy, grad, sigma)
      write(*,'(A,2I4,A,ES14.6,A,ES14.6)') 'Kr-O hf/minis: zeff=', param%zeff, ' emiss1=', param%emiss(1), ' E=', energy
   case ("Am")
      mol%num = [95, 8]
      call get_gcp_param(param, mol, method="r2scan3c")
      call get_geometric_counterpoise(mol, param, cutoff, energy, grad, sigma)
      write(*,'(A,2I4,A,ES14.6,A,ES14.6,A,ES14.6)') 'Am-O r2scan3c: zeff=', param%zeff, &
         ' emiss1=', param%emiss(1), ' slater1=', param%slater(1), ' E=', energy
   end select
end program repro_gcp_oob
