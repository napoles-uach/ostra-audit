program probe_tblite_fortran
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, new
   use tblite_context_type, only : context_type
   use tblite_results, only : results_type
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_xtb_h0, only : get_selfenergy
   use tblite_ceh_ceh, only : get_effective_qat
   use tblite_container, only : container_cache, container_type
   use tblite_solvation, only : new_solvation, solvation_type, solvation_input, alpb_input
   use tblite_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell, get_wignerseitz_weights
   implicit none

   real(wp), parameter :: SENTINEL = 999.0_wp
   type(structure_type) :: mol, mol_p, mol_hal
   type(context_type) :: ctx
   type(xtb_calculator) :: calc2, calc1, calc2_solv
   type(wavefunction_type) :: wfn
   type(error_type), allocatable :: err
   type(container_cache) :: cache
   class(solvation_type), allocatable :: solv
   class(container_type), allocatable :: cont
   type(solvation_input) :: sinput
   integer :: num(3), nat, iat, idir, a, b, ih
   real(wp) :: xyz(3, 3), lat(3, 3), xyz_tmp(3, 3), lat_tmp(3, 3), eps_mat(3, 3)
   logical :: pbc(3)
   real(wp) :: e_ref, g_ref(3, 3), s_ref(3, 3), e, g(3, 3), s(3, 3)
   real(wp) :: ep1, em1, ep2, em2, fd_val, max_fd_g, max_fd_s, h
   real(wp) :: hs(3) = [1.0e-3_wp, 1.0e-4_wp, 1.0e-5_wp]
   real(wp) :: en_ref(3), en(3)
   real(wp) :: cn_en(3), dcn_endr(3, 3, 3), dcn_endL(3, 3, 3)
   real(wp) :: qat_ref(3, 1), dqatdr_ref(3, 3, 3, 1), dqatdL_ref(3, 3, 3, 1)
   real(wp) :: qat(3, 1), dqatdr(3, 3, 3, 1), dqatdL(3, 3, 3, 1)
   real(wp), allocatable :: se_ref(:), dsedcn_ref(:), dsedcn_en(:)
   type(wignerseitz_cell) :: ws
   real(wp), allocatable :: w_ref(:), dwdr_ref(:, :), dwdL_ref(:, :, :), w(:), dwdr(:, :), dwdL(:, :, :)
   integer :: num_hal(3)
   real(wp) :: xyz_hal(3, 3)

   nat = 3
   num = [8, 1, 1]
   xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   xyz(:, 2) = [1.80_wp, 0.35_wp, 0.15_wp]
   xyz(:, 3) = [-0.45_wp, 1.75_wp, -0.25_wp]
   lat(:, 1) = [12.0_wp, 0.2_wp, 0.1_wp]
   lat(:, 2) = [0.3_wp, 12.2_wp, 0.15_wp]
   lat(:, 3) = [0.1_wp, 0.25_wp, 12.5_wp]
   pbc = [.false., .false., .false.]
   call new(mol, num, xyz)

   call new_gfn2_calculator(calc2, mol, err)
   if (allocated(err)) error stop err%message
   call new_gfn1_calculator(calc1, mol, err)
   if (allocated(err)) error stop err%message

   print '(a)', "=== 1. tblite: xtb_singlepoint (GFN2-xTB, 2^2 combinations of gradient/sigma) ==="
   call new_wavefunction(wfn, mol%nat, calc2%bas%nsh, calc2%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   e_ref = SENTINEL; g_ref = SENTINEL; s_ref = SENTINEL
   call xtb_singlepoint(ctx, mol, calc2, wfn, 1.0_wp, e_ref, gradient=g_ref, sigma=s_ref, verbosity=0)
   if (ctx%failed()) error stop "GFN2 ref failed"
   print '(a,es16.8,a,es12.5,a,es12.5)', "REF GFN2 E=", e_ref, " max|g|=", maxval(abs(g_ref)), &
      & " max|s|=", maxval(abs(s_ref))

   ! (F, F)
   call new_wavefunction(wfn, mol%nat, calc2%bas%nsh, calc2%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   e = SENTINEL
   call xtb_singlepoint(ctx, mol, calc2, wfn, 1.0_wp, e, verbosity=0)
   print '(a,es12.5)', "  (F,F) |E-ref|=", abs(e - e_ref)

   ! (T, F): gradient ONLY
   call new_wavefunction(wfn, mol%nat, calc2%bas%nsh, calc2%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   e = SENTINEL; g = SENTINEL
   call xtb_singlepoint(ctx, mol, calc2, wfn, 1.0_wp, e, gradient=g, verbosity=0)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  (T,F) grad only:  |E-ref|=", abs(e - e_ref), &
      & " max|g-ref|=", maxval(abs(g - g_ref)), " min(g)=", minval(g)

   ! (F, T): sigma ONLY
   call new_wavefunction(wfn, mol%nat, calc2%bas%nsh, calc2%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   e = SENTINEL; s = SENTINEL
   call xtb_singlepoint(ctx, mol, calc2, wfn, 1.0_wp, e, sigma=s, verbosity=0)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  (F,T) sigma only: |E-ref|=", abs(e - e_ref), &
      & " max|s-ref|=", maxval(abs(s - s_ref)), " min(s)=", minval(s)

   ! 4pt FD oracle for GFN2 gradient (non-periodic) & virial (periodic)
   do ih = 1, 3
      h = hs(ih)
      max_fd_g = 0.0_wp
      do iat = 1, nat
         do idir = 1, 3
            xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) + h
            call new(mol_p, num, xyz_tmp)
            call new_wavefunction(wfn, nat, calc2%bas%nsh, calc2%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
            call xtb_singlepoint(ctx, mol_p, calc2, wfn, 1.0e-3_wp, ep1, verbosity=0)
            xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) - h
            call new(mol_p, num, xyz_tmp)
            call new_wavefunction(wfn, nat, calc2%bas%nsh, calc2%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
            call xtb_singlepoint(ctx, mol_p, calc2, wfn, 1.0e-3_wp, em1, verbosity=0)
            xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) + 2.0_wp*h
            call new(mol_p, num, xyz_tmp)
            call new_wavefunction(wfn, nat, calc2%bas%nsh, calc2%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
            call xtb_singlepoint(ctx, mol_p, calc2, wfn, 1.0e-3_wp, ep2, verbosity=0)
            xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) - 2.0_wp*h
            call new(mol_p, num, xyz_tmp)
            call new_wavefunction(wfn, nat, calc2%bas%nsh, calc2%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
            call xtb_singlepoint(ctx, mol_p, calc2, wfn, 1.0e-3_wp, em2, verbosity=0)
            fd_val = (-ep2 + 8.0_wp*ep1 - 8.0_wp*em1 + em2) / (12.0_wp*h)
            max_fd_g = max(max_fd_g, abs(fd_val - g_ref(idir, iat)))
         end do
      end do
      print '(a,es8.1,a,es11.3)', "  4pt-FD GFN2 h=", h, " max|g_ref-FD|=", max_fd_g
   end do

   print '(a)', ""
   print '(a)', "=== 2. tblite: xtb_singlepoint with ALPB/CDS solvation (gradient only -> adds to uninit!) ==="
   call new_gfn2_calculator(calc2_solv, mol, err)
   allocate(sinput%alpb)
   sinput%alpb%solvent = "water"
   call new_solvation(solv, mol, sinput, err, method="gfn2")
   if (allocated(err)) error stop err%message
   call move_alloc(solv, cont)
   call calc2_solv%push_back(cont)
   call new_wavefunction(wfn, nat, calc2_solv%bas%nsh, calc2_solv%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   g_ref = SENTINEL; s_ref = SENTINEL
   call xtb_singlepoint(ctx, mol, calc2_solv, wfn, 1.0_wp, e_ref, gradient=g_ref, sigma=s_ref, verbosity=0)
   call new_wavefunction(wfn, nat, calc2_solv%bas%nsh, calc2_solv%bas%nao, 1, 300.0_wp * 3.166808578545117e-6_wp)
   g = SENTINEL
   call xtb_singlepoint(ctx, mol, calc2_solv, wfn, 1.0_wp, e, gradient=g, verbosity=0)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  GFN2+ALPB (grad only): max|g-ref|=", &
      & maxval(abs(g - g_ref)), " min(g)=", minval(g), " max(g)=", maxval(g)

   print '(a)', ""
   print '(a)', "=== 3. tblite: individual containers get_engrad (repulsion, d3, d4, halogen) ==="
   ! Repulsion (from calc2)
   call calc2%repulsion%update(mol, cache)
   en_ref = 0.0_wp; g_ref = 0.0_wp; s_ref = 0.0_wp
   ! Note: container%get_engrad has intent(inout) (accumulates into gradient/sigma).
   ! Initialize to 0.0_wp so exact comparison against ref (also from 0.0_wp) has zero cancellation.
   call calc2%repulsion%get_engrad(mol, cache, en_ref, g_ref, s_ref)
   en = 0.0_wp; g = 0.0_wp
   call calc2%repulsion%get_engrad(mol, cache, en, gradient=g)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  repulsion%get_engrad (grad only): ref=", &
      & maxval(abs(g_ref)), " max|g-ref|=", maxval(abs(g - g_ref)), " max|g|=", maxval(abs(g))
   en = 0.0_wp; s = 0.0_wp
   call calc2%repulsion%get_engrad(mol, cache, en, sigma=s)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  repulsion%get_engrad (sigma only):ref=", &
      & maxval(abs(s_ref)), " max|s-ref|=", maxval(abs(s - s_ref)), " max|s|=", maxval(abs(s))

   ! D4 dispersion (from calc2)
   call calc2%dispersion%update(mol, cache)
   en_ref = 0.0_wp; g_ref = 0.0_wp; s_ref = 0.0_wp
   call calc2%dispersion%get_engrad(mol, cache, en_ref, g_ref, s_ref)
   en = 0.0_wp; g = 0.0_wp
   call calc2%dispersion%get_engrad(mol, cache, en, gradient=g)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  d4_disp%get_engrad   (grad only): ref=", &
      & maxval(abs(g_ref)), " max|g-ref|=", maxval(abs(g - g_ref)), " max|g|=", maxval(abs(g))
   en = 0.0_wp; s = 0.0_wp
   call calc2%dispersion%get_engrad(mol, cache, en, sigma=s)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  d4_disp%get_engrad   (sigma only):ref=", &
      & maxval(abs(s_ref)), " max|s-ref|=", maxval(abs(s - s_ref)), " max|s|=", maxval(abs(s))

   ! D3 dispersion (from calc1)
   call calc1%dispersion%update(mol, cache)
   en_ref = 0.0_wp; g_ref = 0.0_wp; s_ref = 0.0_wp
   call calc1%dispersion%get_engrad(mol, cache, en_ref, g_ref, s_ref)
   en = 0.0_wp; g = 0.0_wp
   call calc1%dispersion%get_engrad(mol, cache, en, gradient=g)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  d3_disp%get_engrad   (grad only): ref=", &
      & maxval(abs(g_ref)), " max|g-ref|=", maxval(abs(g - g_ref)), " max|g|=", maxval(abs(g))
   en = 0.0_wp; s = 0.0_wp
   call calc1%dispersion%get_engrad(mol, cache, en, sigma=s)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  d3_disp%get_engrad   (sigma only):ref=", &
      & maxval(abs(s_ref)), " max|s-ref|=", maxval(abs(s - s_ref)), " max|s|=", maxval(abs(s))

   ! Halogen correction (from calc1 on C-Br...N system)
   num_hal = [6, 35, 7]
   xyz_hal(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   xyz_hal(:, 2) = [3.6_wp, 0.0_wp, 0.0_wp]
   xyz_hal(:, 3) = [9.0_wp, 0.2_wp, 0.0_wp]
   call new(mol_hal, num_hal, xyz_hal)
   call new_gfn1_calculator(calc1, mol_hal, err)
   call calc1%halogen%update(mol_hal, cache)
   en_ref = 0.0_wp; g_ref = 0.0_wp; s_ref = 0.0_wp
   call calc1%halogen%get_engrad(mol_hal, cache, en_ref, g_ref, s_ref)
   en = 0.0_wp; g = 0.0_wp
   call calc1%halogen%get_engrad(mol_hal, cache, en, gradient=g)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  halogen%get_engrad   (grad only): ref=", &
      & maxval(abs(g_ref)), " max|g-ref|=", maxval(abs(g - g_ref)), " max|g|=", maxval(abs(g))
   en = 0.0_wp; s = 0.0_wp
   call calc1%halogen%get_engrad(mol_hal, cache, en, sigma=s)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  halogen%get_engrad   (sigma only):ref=", &
      & maxval(abs(s_ref)), " max|s-ref|=", maxval(abs(s - s_ref)), " max|s|=", maxval(abs(s))

   print '(a)', ""
   print '(a)', "=== 4. tblite: get_effective_qat (CEH) & get_selfenergy (h0) ==="
   cn_en = [1.5_wp, 0.8_wp, 0.8_wp]
   dcn_endr = 0.25_wp; dcn_endL = 0.35_wp
   call get_effective_qat(mol, cn_en, qat_ref, dcn_endr, dcn_endL, dqatdr_ref, dqatdL_ref)
   dqatdr = SENTINEL
   call get_effective_qat(mol, cn_en, qat, dcn_endr=dcn_endr, dqatdr=dqatdr)
   print '(a,es12.5,a,es12.5,a,es12.5)', "  get_effective_qat (dqatdr only): ref=", &
      & maxval(abs(dqatdr_ref)), " max|err|=", maxval(abs(dqatdr - dqatdr_ref)), " max|val|=", maxval(abs(dqatdr))

   allocate(se_ref(calc2%bas%nsh), dsedcn_ref(calc2%bas%nsh), dsedcn_en(calc2%bas%nsh))
   dsedcn_en = SENTINEL
   call get_selfenergy(calc2%h0, mol%id, calc2%bas%ish_at, calc2%bas%nsh_id, &
      & cn=cn_en, selfenergy=se_ref, dsedcn=dsedcn_ref, dsedcn_en=dsedcn_en)
   print '(a,es12.5)', "  get_selfenergy (dsedcn_en w/o cn_en): min(dsedcn_en)=", minval(dsedcn_en)

   print '(a)', ""
   print '(a)', "=== 5. tblite CONTROL: get_wignerseitz_weights ==="
   call new(mol_p, num, xyz, lattice=lat, periodic=[.true., .true., .true.])
   call new_wignerseitz_cell(ws, mol_p)
   allocate(w_ref(maxval(ws%nimg)), dwdr_ref(3, maxval(ws%nimg)), dwdL_ref(3, 3, maxval(ws%nimg)))
   allocate(w(maxval(ws%nimg)), dwdr(3, maxval(ws%nimg)), dwdL(3, 3, maxval(ws%nimg)))
   call get_wignerseitz_weights(ws, 2, 1, xyz(:, 2) - xyz(:, 1), w_ref, dwdr_ref, dwdL_ref)
   dwdr = SENTINEL
   call get_wignerseitz_weights(ws, 2, 1, xyz(:, 2) - xyz(:, 1), w, dwdr=dwdr)
   print '(a,es12.5)', "  wignerseitz (dwdr only): max|err|=", maxval(abs(dwdr - dwdr_ref))
end program probe_tblite_fortran
