module multicharge_ewald
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_wignerseitz, only : wignerseitz_cell
   implicit none
contains
   subroutine get_alpha(lattice, alpha, error)
      real(wp), intent(in) :: lattice(3, 3)
      real(wp), intent(out) :: alpha
      type(error_type), allocatable, intent(out), optional :: error
      alpha = 0.0_wp
   end subroutine get_alpha
   subroutine get_ewald_amat(mol, wsc, alpha, rad, amat)
      class(structure_type), intent(in) :: mol
      type(wignerseitz_cell), intent(in) :: wsc
      real(wp), intent(in) :: alpha, rad(:)
      real(wp), intent(out) :: amat(:, :)
      amat = 0.0_wp
   end subroutine get_ewald_amat
   subroutine get_ewald_damat(mol, wsc, alpha, rad, qvec, dadr, dadL, atrace)
      class(structure_type), intent(in) :: mol
      type(wignerseitz_cell), intent(in) :: wsc
      real(wp), intent(in) :: alpha, rad(:), qvec(:)
      real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
      dadr = 0.0_wp; dadL = 0.0_wp; atrace = 0.0_wp
   end subroutine get_ewald_damat
   subroutine add_ewald_hessian(mol, wsc, alpha, rad, qvec, adj, hessian, diagonal)
      class(structure_type), intent(in) :: mol
      type(wignerseitz_cell), intent(in) :: wsc
      real(wp), intent(in) :: alpha, rad(:), qvec(:), adj(:)
      real(wp), intent(inout) :: hessian(:, :)
      real(wp), intent(in), optional :: diagonal(:, :, :)
   end subroutine add_ewald_hessian
   subroutine get_ewald_amat_bc(mol, wsc, alpha, rad, amat)
      class(structure_type), intent(in) :: mol
      type(wignerseitz_cell), intent(in) :: wsc
      real(wp), intent(in) :: alpha, rad(:)
      real(wp), intent(out) :: amat(:, :)
      amat = 0.0_wp
   end subroutine get_ewald_amat_bc
end module multicharge_ewald
