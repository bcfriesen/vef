! NLTE statistical equilibrium equations for H I with levels 1 and 2 in
! NLTE, and levels 3 through 16 in LTE. taken from Auer & Mihalas, ApJ,
! 156, 157 (1969)
subroutine solve_statistical_equations(layer)
  use precision_mod
  use global, only: n_e, n_i_star
  use interfaces, only: c_ij, c_ji, c_ik, c_ki, &
                        r_ij, r_ik, r_ki, &
                        z_ji
  implicit none
  real(kind=dp), dimension(3,3) :: rate_matrix
  real(kind=dp), dimension(3) :: rhs
  ! # density of atoms in LTE levels (3-16)
  real(kind=dp), dimension(layer) :: n_u
  integer :: layer

  ! components of RHS are:
  ! 1.) n_1
  ! 2.) n_2
  ! 3.) n_H

  rate_matrix(1,1) = r_ik(1) + n_e(layer) * (c_ik(1) + c_ij(1,2))
  rate_matrix(1,2) = -(z_ji(2,1) + n_e(layer) * c_ji(2,1))
  rate_matrix(1,3) = 0.0
  rhs(1) = n_i_star(layer, 1) * (r_ki(1) + n_e(layer) * c_ki(1))

  rate_matrix(2,1) = -n_e(layer) * c_ij(1,2)
  rate_matrix(2,2) = z_ji(2,1) + n_e(layer) * c_ji(2,1) + r_ik(2) + n_e(layer) * c_ik(2)
  rate_matrix(2,3) = 0.0
  rhs(2) = n_i_star(layer, 2) * (r_ki(2) + n_e(layer) * c_ki(2))

  ! have to move stuff around since they didn't write this equation with
  ! unknowns on LHS and knowns on RHS
  rate_matrix(3,1) = 1.0
  rate_matrix(3,2) = 1.0
  rate_matrix(3,3) = -1.0
  rhs(3) = -n_u(layer) - n_e(layer)

end subroutine solve_statistical_equations
