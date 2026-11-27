program p29f
  use mctc_env, only: wp, error_type
  use mctc_io, only: structure_type, new
  use dftd4, only: d4_model, new_d4_model, rational_damping_param, get_rational_damping, &
     & get_dispersion, realspace_cutoff, damping_param
  implicit none
  type(structure_type) :: mol
  type(d4_model) :: d4
  type(error_type), allocatable :: err
  class(damping_param), allocatable :: par
  real(wp) :: xyz(3,4), e0, e1, g0(3,4), s0(3,3), g1(3,4)
  xyz = reshape([0.0_wp,0.0_wp,0.22_wp, 0.0_wp,1.77_wp,-0.51_wp, 1.53_wp,-0.885_wp,-0.51_wp, -1.53_wp,-0.885_wp,-0.51_wp],[3,4])
  call new(mol, [7,1,1,1], xyz)
  call new_d4_model(err, d4, mol)
  call get_rational_damping("pbe0", par, s9=1.0_wp)
  call get_dispersion(mol, d4, par, realspace_cutoff(), e0, g0, s0)
  g1 = 999.0_wp
  call get_dispersion(mol, d4, par, realspace_cutoff(), e1, gradient=g1)
  print '(a,es10.3)', 'Fortran gradient only: max|g - g(both)| = ', maxval(abs(g1-g0))
end program
