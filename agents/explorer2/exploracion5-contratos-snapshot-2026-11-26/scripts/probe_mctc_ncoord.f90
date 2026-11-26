program probe_mctc_ncoord
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, new
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   use mctc_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell, get_wignerseitz_weights
   use mctc_cutoff, only : get_lattice_points
   use mctc_csrlist, only : csr_list, new_csr_list
   implicit none

   real(wp), parameter :: SENTINEL = 999.0_wp
   type(structure_type) :: mol, mol_p
   type(error_type), allocatable :: err
   class(ncoord_type), allocatable :: ncoord
   real(wp), allocatable :: lattr(:, :), cn_ref(:), dcndr_ref(:, :, :), dcndL_ref(:, :, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: dcndrlist_ref(:, :), dcndrlist(:, :)
   type(csr_list) :: list
   integer :: num(3)
   real(wp) :: xyz(3, 3), lat(3, 3)
   logical :: pbc(3)
   integer :: ic, nat, iat, idir, a, b, ih
   real(wp) :: h, fd_val, max_fd_err_dr, max_fd_err_dL
   real(wp) :: hs(4) = [1.0e-2_wp, 1.0e-3_wp, 1.0e-4_wp, 1.0e-5_wp]
   real(wp) :: cn_p1(3), cn_m1(3), cn_p2(3), cn_m2(3)
   real(wp) :: xyz_tmp(3, 3), lat_tmp(3, 3), eps_mat(3, 3)
   character(len=16) :: names(3) = ["exp             ", "dftd4           ", "erf_cut6        "]
   integer :: ids(3) = [cn_count%exp, cn_count%dftd4, cn_count%erf]
   type(wignerseitz_cell) :: ws
   real(wp), allocatable :: w_ref(:), dwdr_ref(:, :), dwdL_ref(:, :, :)
   real(wp), allocatable :: w(:), dwdr(:, :), dwdL(:, :, :)
   integer :: max_img

   ! Asymmetric 3-atom H2O-like periodic cell so all dcndr and dcndL components are non-zero
   num = [8, 1, 1]
   xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   xyz(:, 2) = [1.80_wp, 0.35_wp, 0.15_wp]
   xyz(:, 3) = [-0.45_wp, 1.75_wp, -0.25_wp]
   lat(:, 1) = [8.0_wp, 0.2_wp, 0.1_wp]
   lat(:, 2) = [0.3_wp, 8.2_wp, 0.15_wp]
   lat(:, 3) = [0.1_wp, 0.25_wp, 8.5_wp]
   pbc = [.true., .true., .true.]
   nat = 3

   call new(mol, num, xyz, lattice=lat, periodic=pbc)
   allocate(cn_ref(nat), dcndr_ref(3, nat, nat), dcndL_ref(3, 3, nat))
   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))

   print '(a)', "=== 1. mctc-lib: ncoord%get_cn / get_coordination_number (2^2 combinations) ==="
   do ic = 1, 3
      if (allocated(ncoord)) deallocate(ncoord)
      if (ic == 3) then
         call new_ncoord(ncoord, mol, ids(ic), err, cut=6.0_wp)
      else
         call new_ncoord(ncoord, mol, ids(ic), err)
      end if
      if (allocated(err)) error stop err%message

      ! Reference (T, T)
      cn_ref = SENTINEL; dcndr_ref = SENTINEL; dcndL_ref = SENTINEL
      call ncoord%get_cn(mol, cn_ref, err, dcndr=dcndr_ref, dcndL=dcndL_ref)
      if (allocated(err)) error stop err%message
      print '(a,a,3(1x,es12.5))', "Model ", trim(names(ic)), cn_ref
      print '(a,es12.5,a,es12.5)', "  ref max|dcndr|=", maxval(abs(dcndr_ref)), &
         & " max|dcndL|=", maxval(abs(dcndL_ref))

      ! (F, F)
      cn = SENTINEL
      call ncoord%get_cn(mol, cn, err)
      print '(a,es12.5)', "  (F,F) max|cn-ref|     = ", maxval(abs(cn - cn_ref))

      ! (T, F): dcndr ONLY
      cn = SENTINEL; dcndr = SENTINEL
      call ncoord%get_cn(mol, cn, err, dcndr=dcndr)
      print '(a,es12.5,a,es12.5,a,es12.5)', "  (T,F) dcndr only: max|err|=", &
         & maxval(abs(dcndr - dcndr_ref)), " min(dcndr)=", minval(dcndr), " max(dcndr)=", maxval(dcndr)

      ! (F, T): dcndL ONLY
      cn = SENTINEL; dcndL = SENTINEL
      call ncoord%get_cn(mol, cn, err, dcndL=dcndL)
      print '(a,es12.5,a,es12.5,a,es12.5)', "  (F,T) dcndL only: max|err|=", &
         & maxval(abs(dcndL - dcndL_ref)), " min(dcndL)=", minval(dcndL), " max(dcndL)=", maxval(dcndL)

      ! Physical oracle: 4-point FD step scan for dcndr and dcndL
      do ih = 1, 4
         h = hs(ih)
         max_fd_err_dr = 0.0_wp
         do iat = 1, nat
            do idir = 1, 3
               xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) + h
               call new(mol_p, num, xyz_tmp, lattice=lat, periodic=pbc)
               call ncoord%get_cn(mol_p, cn_p1, err)
               xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) - h
               call new(mol_p, num, xyz_tmp, lattice=lat, periodic=pbc)
               call ncoord%get_cn(mol_p, cn_m1, err)
               xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) + 2.0_wp*h
               call new(mol_p, num, xyz_tmp, lattice=lat, periodic=pbc)
               call ncoord%get_cn(mol_p, cn_p2, err)
               xyz_tmp = xyz; xyz_tmp(idir, iat) = xyz(idir, iat) - 2.0_wp*h
               call new(mol_p, num, xyz_tmp, lattice=lat, periodic=pbc)
               call ncoord%get_cn(mol_p, cn_m2, err)
               do a = 1, nat
                  fd_val = (-cn_p2(a) + 8.0_wp*cn_p1(a) - 8.0_wp*cn_m1(a) + cn_m2(a)) / (12.0_wp*h)
                  max_fd_err_dr = max(max_fd_err_dr, abs(fd_val - dcndr_ref(idir, iat, a)))
               end do
            end do
         end do
         max_fd_err_dL = 0.0_wp
         do a = 1, 3
            do b = 1, 3
               eps_mat = 0.0_wp; eps_mat(a, b) = h
               xyz_tmp = xyz + matmul(eps_mat, xyz); lat_tmp = lat + matmul(eps_mat, lat)
               call new(mol_p, num, xyz_tmp, lattice=lat_tmp, periodic=pbc)
               call ncoord%get_cn(mol_p, cn_p1, err)
               eps_mat = 0.0_wp; eps_mat(a, b) = -h
               xyz_tmp = xyz + matmul(eps_mat, xyz); lat_tmp = lat + matmul(eps_mat, lat)
               call new(mol_p, num, xyz_tmp, lattice=lat_tmp, periodic=pbc)
               call ncoord%get_cn(mol_p, cn_m1, err)
               eps_mat = 0.0_wp; eps_mat(a, b) = 2.0_wp*h
               xyz_tmp = xyz + matmul(eps_mat, xyz); lat_tmp = lat + matmul(eps_mat, lat)
               call new(mol_p, num, xyz_tmp, lattice=lat_tmp, periodic=pbc)
               call ncoord%get_cn(mol_p, cn_p2, err)
               eps_mat = 0.0_wp; eps_mat(a, b) = -2.0_wp*h
               xyz_tmp = xyz + matmul(eps_mat, xyz); lat_tmp = lat + matmul(eps_mat, lat)
               call new(mol_p, num, xyz_tmp, lattice=lat_tmp, periodic=pbc)
               call ncoord%get_cn(mol_p, cn_m2, err)
               do iat = 1, nat
                  fd_val = (-cn_p2(iat) + 8.0_wp*cn_p1(iat) - 8.0_wp*cn_m1(iat) + cn_m2(iat)) / (12.0_wp*h)
                  max_fd_err_dL = max(max_fd_err_dL, abs(fd_val - dcndL_ref(a, b, iat)))
               end do
            end do
         end do
         print '(a,es8.1,a,es11.3,a,es11.3)', "  4pt-FD h=", h, " max|dcndr-FD|=", max_fd_err_dr, &
            & " max|dcndL-FD|=", max_fd_err_dL
      end do
   end do

   print '(a)', ""
   print '(a)', "=== 2. mctc-lib: get_coordination_number with CSR list (dcndrlist, dcndL) ==="
   if (allocated(ncoord)) deallocate(ncoord)
   call new_ncoord(ncoord, mol, cn_count%exp, err)
   call get_lattice_points(mol%periodic, mol%lattice, ncoord%cutoff, lattr, err)
   call new_csr_list(list, mol, err, cutoff=ncoord%cutoff)
   if (allocated(err)) error stop err%message
   allocate(dcndrlist_ref(3, size(list%nlat)), dcndrlist(3, size(list%nlat)))
   cn_ref = SENTINEL; dcndrlist_ref = SENTINEL; dcndL_ref = SENTINEL
   call ncoord%get_coordination_number(mol, lattr, cn_ref, dcndL=dcndL_ref, list=list, dcndrlist=dcndrlist_ref)
   print '(a,es12.5,a,es12.5)', "  ref max|dcndrlist|=", maxval(abs(dcndrlist_ref)), &
      & " max|dcndL|=", maxval(abs(dcndL_ref))
   cn = SENTINEL; dcndrlist = SENTINEL
   call ncoord%get_coordination_number(mol, lattr, cn, list=list, dcndrlist=dcndrlist)
   print '(a,es12.5,a,es12.5)', "  (T,F) dcndrlist only: max|err|=", &
      & maxval(abs(dcndrlist - dcndrlist_ref)), " min(val)=", minval(dcndrlist)
   cn = SENTINEL; dcndL = SENTINEL
   call ncoord%get_coordination_number(mol, lattr, cn, dcndL=dcndL, list=list)
   print '(a,es12.5,a,es12.5)', "  (F,T) dcndL only:     max|err|=", &
      & maxval(abs(dcndL - dcndL_ref)), " min(val)=", minval(dcndL)

   print '(a)', ""
   print '(a)', "=== 3. mctc-lib CONTROL: get_wignerseitz_weights (dwdr, dwdL) ==="
   call new_wignerseitz_cell(ws, mol)
   max_img = maxval(ws%nimg)
   allocate(w_ref(max_img), dwdr_ref(3, max_img), dwdL_ref(3, 3, max_img))
   allocate(w(max_img), dwdr(3, max_img), dwdL(3, 3, max_img))
   w_ref = SENTINEL; dwdr_ref = SENTINEL; dwdL_ref = SENTINEL
   call get_wignerseitz_weights(ws, 2, 1, xyz(:, 2) - xyz(:, 1), w_ref, dwdr=dwdr_ref, dwdL=dwdL_ref)
   w = SENTINEL; dwdr = SENTINEL
   call get_wignerseitz_weights(ws, 2, 1, xyz(:, 2) - xyz(:, 1), w, dwdr=dwdr)
   print '(a,es12.5,a,es12.5)', "  (T,F) dwdr only: max|w-ref|=", maxval(abs(w - w_ref)), &
      & " max|dwdr-ref|=", maxval(abs(dwdr - dwdr_ref))
   w = SENTINEL; dwdL = SENTINEL
   call get_wignerseitz_weights(ws, 2, 1, xyz(:, 2) - xyz(:, 1), w, dwdL=dwdL)
   print '(a,es12.5,a,es12.5)', "  (F,T) dwdL only: max|w-ref|=", maxval(abs(w - w_ref)), &
      & " max|dwdL-ref|=", maxval(abs(dwdL - dwdL_ref))
end program probe_mctc_ncoord
