module mctc_env
   use iso_fortran_env, only : sp => real32, dp => real64, wp => real64, i4 => int32, i8 => int64
   implicit none
   public :: sp, dp, wp, i4, i8, error_type, fatal_error
   type :: error_type
      character(len=:), allocatable :: message
   end type error_type
contains
   subroutine fatal_error(error, message)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: message
      allocate(error)
      error%message = message
   end subroutine fatal_error
end module mctc_env

module mctc_data
   use mctc_env, only : wp
   implicit none
contains
   elemental function get_covalent_rad(num) result(rad)
      integer, intent(in) :: num
      real(wp) :: rad
      rad = 1.5_wp
   end function get_covalent_rad
   elemental function get_pauling_en(num) result(en)
      integer, intent(in) :: num
      real(wp) :: en
      en = 2.0_wp
   end function get_pauling_en
   elemental function get_vdw_rad(iz, jz) result(rad)
      integer, intent(in) :: iz, jz
      real(wp) :: rad
      rad = 3.0_wp
   end function get_vdw_rad
end module mctc_data

module mctc_io_constants
   use mctc_env, only : wp
   implicit none
   real(wp), parameter :: pi = 3.1415926535897932384626433832795_wp
end module mctc_io_constants

module mctc_io_convert
   use mctc_env, only : wp
   implicit none
   real(wp), parameter :: autoaa = 0.529177210903_wp
   real(wp), parameter :: aatoau = 1.0_wp / autoaa
end module mctc_io_convert

module mctc_io_symbols
   implicit none
contains
   pure function to_number(sym) result(num)
      character(len=*), intent(in) :: sym
      integer :: num
      num = 0
   end function to_number
end module mctc_io_symbols

module mctc_io_math
   use mctc_env, only : wp
   implicit none
contains
   pure function matdet_3x3(a) result(d)
      real(wp), intent(in) :: a(3, 3)
      real(wp) :: d
      d = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
        - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
        + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
   end function matdet_3x3
   pure function matinv_3x3(a) result(inv)
      real(wp), intent(in) :: a(3, 3)
      real(wp) :: inv(3, 3)
      inv = 0.0_wp
      inv(1,1) = 1.0_wp; inv(2,2) = 1.0_wp; inv(3,3) = 1.0_wp
   end function matinv_3x3
   pure function outerprod(a, b) result(c)
      real(wp), intent(in) :: a(:), b(:)
      real(wp) :: c(size(a), size(b))
      integer :: i, j
      do j = 1, size(b)
         do i = 1, size(a)
            c(i, j) = a(i) * b(j)
         end do
      end do
   end function outerprod
end module mctc_io_math

module mctc_io
   use mctc_env, only : wp
   implicit none
   type :: structure_type
      integer :: nat = 0
      integer :: nid = 0
      integer, allocatable :: id(:)
      integer, allocatable :: num(:)
      real(wp), allocatable :: xyz(:, :)
      real(wp) :: lattice(3, 3) = 0.0_wp
      logical :: periodic(3) = .false.
      real(wp) :: charge = 0.0_wp
   end type structure_type
end module mctc_io

module mctc_cutoff
   use mctc_env, only : wp, error_type
   implicit none
   interface get_lattice_points
      module procedure glp_cut, glp_rep
   end interface
contains
   subroutine glp_cut(periodic, lattice, cutoff, lattr, error)
      logical, intent(in) :: periodic(3)
      real(wp), intent(in) :: lattice(3, 3)
      real(wp), intent(in) :: cutoff
      real(wp), allocatable, intent(out) :: lattr(:, :)
      type(error_type), allocatable, intent(out), optional :: error
      allocate(lattr(3, 1))
      lattr(:, 1) = 0.0_wp
   end subroutine glp_cut
   subroutine glp_rep(lattice, rep, origin, lattr)
      real(wp), intent(in) :: lattice(3, 3)
      integer, intent(in) :: rep(3)
      logical, intent(in) :: origin
      real(wp), allocatable, intent(out) :: lattr(:, :)
      allocate(lattr(3, 1))
      lattr(:, 1) = 0.0_wp
   end subroutine glp_rep
end module mctc_cutoff

module mctc_wignerseitz
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   implicit none
   type :: wignerseitz_cell
      integer, allocatable :: nimg(:, :)
      integer, allocatable :: tridx(:, :, :)
      real(wp), allocatable :: trans(:, :)
   end type wignerseitz_cell
contains
   subroutine new_wignerseitz_cell(wsc, mol, error)
      type(wignerseitz_cell), intent(out) :: wsc
      class(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out), optional :: error
      allocate(wsc%nimg(mol%nat, mol%nat), source=1)
      allocate(wsc%tridx(1, mol%nat, mol%nat), source=1)
      allocate(wsc%trans(3, 1), source=0.0_wp)
   end subroutine new_wignerseitz_cell
end module mctc_wignerseitz

module mctc_ncoord
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   implicit none
   type :: cn_count_type
      integer :: exp = 1
      integer :: erf = 2
      integer :: erf_en = 3
   end type cn_count_type
   type(cn_count_type), parameter :: cn_count = cn_count_type()
   type :: ncoord_type
      real(wp) :: cutoff = 25.0_wp
   contains
      procedure :: get_coordination_number => stub_get_cn
      procedure :: add_coordination_number_hessian => stub_add_cn_hess
   end type ncoord_type
contains
   subroutine new_ncoord(self, mol, cn_type, error, cutoff, kcn, rcov, en, cut, shift, norm_exp)
      class(ncoord_type), allocatable, intent(out) :: self
      class(structure_type), intent(in) :: mol
      integer, intent(in) :: cn_type
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: cutoff, kcn, cut, shift, norm_exp
      real(wp), intent(in), optional :: rcov(:), en(:)
      allocate(self)
      if (present(cutoff)) self%cutoff = cutoff
   end subroutine new_ncoord
   subroutine stub_get_cn(self, mol, trans, cn, dcndr, dcndL)
      class(ncoord_type), intent(in) :: self
      class(structure_type), intent(in) :: mol
      real(wp), intent(in) :: trans(:, :)
      real(wp), intent(out) :: cn(:)
      real(wp), intent(out), optional :: dcndr(:, :, :)
      real(wp), intent(out), optional :: dcndL(:, :, :)
      integer :: i, j
      real(wp) :: r, r0, ex, dcn_dr, vec(3)
      cn = 0.0_wp
      if (present(dcndr)) dcndr = 0.0_wp
      if (present(dcndL)) dcndL = 0.0_wp
      r0 = 2.5_wp
      do i = 1, mol%nat
         do j = i + 1, mol%nat
            vec = mol%xyz(:, i) - mol%xyz(:, j)
            r = norm2(vec)
            ex = exp(-16.0_wp * (r0 / r - 1.0_wp))
            cn(i) = cn(i) + 1.0_wp / (1.0_wp + ex)
            cn(j) = cn(j) + 1.0_wp / (1.0_wp + ex)
            if (present(dcndr)) then
               dcn_dr = (16.0_wp * r0 / (r*r)) * ex / (1.0_wp + ex)**2
               dcndr(:, i, i) = dcndr(:, i, i) + dcn_dr * vec / r
               dcndr(:, j, i) = dcndr(:, j, i) - dcn_dr * vec / r
               dcndr(:, i, j) = dcndr(:, i, j) + dcn_dr * vec / r
               dcndr(:, j, j) = dcndr(:, j, j) - dcn_dr * vec / r
            end if
         end do
      end do
   end subroutine stub_get_cn
   subroutine stub_add_cn_hess(self, mol, trans, dEdcn, hessian)
      class(ncoord_type), intent(in) :: self
      class(structure_type), intent(in) :: mol
      real(wp), intent(in) :: trans(:, :), dEdcn(:)
      real(wp), intent(inout) :: hessian(:, :)
   end subroutine stub_add_cn_hess
end module mctc_ncoord
