module mctc_env
   use iso_fortran_env, only : wp => real64, i8 => int64
   implicit none
   public :: wp, i8, error_type, fatal_error
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

module mctc_io
   use mctc_env, only : wp
   implicit none
   public :: structure_type
   type :: structure_type
      integer :: nat = 0, nid = 0
      integer, allocatable :: id(:), num(:)
      real(wp), allocatable :: xyz(:, :)
      real(wp) :: lattice(3, 3) = 0.0_wp
      logical :: periodic(3) = .false.
   end type structure_type
end module mctc_io

module dftd3_data
   use mctc_env, only : wp
   implicit none
   public :: get_vdw_rad
contains
   pure function get_vdw_rad(iz, jz) result(rvdw)
      integer, intent(in) :: iz, jz
      real(wp) :: rvdw
      rvdw = 3.0_wp
   end function get_vdw_rad
end module dftd3_data
