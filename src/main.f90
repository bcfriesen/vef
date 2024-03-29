! A code for solving the radiative transfer equation in a plane-parallel
! atmosphere using variable Eddington factors. See article "Difference
! Equations and Linearization Methods" by Lawrence Auer
! in "Methods in Radiative Transfer" ed. Wolfgang Kalkofen (1984).

program main
  use precision_mod
  use interfaces
  use global
  use const, only: a2cm
  implicit none

  integer :: i1, i2, i3

  ! In order to sample the radiation field most efficently the optical depth
  ! grid should be spaced logarithmically, such that
  !
  !                           \tau_i = \tau_{i-1} * f, 
  !
  ! where f is this constant, to be calculated below.
  real(kind=dp) :: tau_grid_ratio

  ! We need to set the first two tau points by hand. The outermost layer should
  ! be
  !
  !                             \tau_1 = 0
  !
  ! and the next outermost layer should be
  !
  !                         \tau_2 = \tau_{min}
  !
  ! Then each deeper point in the grid will just be some constant multiplicative
  ! factor larger than the previous point.
  tau_grid( 1 ) = 0.0d+0
  tau_grid( 2 ) = tau_min
  tau_grid_ratio = ( tau_max / tau_grid( 2 ) )**(1.0d+0 / ( n_depth_pts - 2 ) )
  ! Set up the rest of the optical depth grid.
  do i1 = 3, n_depth_pts
    tau_grid( i1 ) = tau_grid( i1 - 1 ) * tau_grid_ratio
  end do

  ! Set up finite-difference-related grid values.
  do i1 = 1, n_depth_pts - 1
    dtau( i1 ) = calc_dtau( i1 )
  end do
  do i1 = 2, n_depth_pts - 1
    dtau_tilde( i1 ) = calc_dtau_tilde( i1 )
  end do

  ! Set up direction cosine grid. These should be evenly spaced.
  do i1 = 1, n_mu_pts
    mu_grid( i1 ) = real( i1 ) * ( 1.0d+0 / real( n_mu_pts ) )
  end do

  ! Set up wavelength grid.
  do i1 = 1, n_wl_pts
    wl_grid( i1 ) = wlmin + real( i1 - 1 ) * ( wlmax - wlmin ) / &
    real( n_wl_pts - 1 )
  end do

  open( unit = 22, file = 'source_fn.dat' )
  open( unit = 23, file = 'moments.dat' )
  open( unit = 24, file = 'f_K.dat' )
  open( unit = 25, file = 'f_H.dat' )
  open( unit = 30, file = 'spectrum.dat' )

!----------------------START LTE RUN TO GET FIRST GUESS FOR J-------------------
  write( *, * ) 'STARTING LTE MODE'
  write( *, * )

  ! set S = B for LTE
  do i1 = 1, n_depth_pts
    do i2 = 1, n_wl_pts
      source_fn( i1, i2 ) = planck_fn( wl_grid( i2 ) * a2cm, temp )
    end do
  end do
  !call write_source_fn

  ! When we switch to NLTE mode we will need a decent starting guess for S in
  ! order to solve the scattering problem for J. Starting with S = B is a
  ! terrible guess. A better guess would be S_NLTE = S(J_LTE). So in LTE we will
  ! calculate J according to its definition, that is, by integrating the
  ! Feautrier variable j.
  call solve_rte

  call calc_0th_moment_j

  write( *, * )
  write( *, * ) 'LTE COMPLETE!'
  write( *, * ) 

!--------------------------------END LTE RUN------------------------------------



!-------------------------------START NLTE RUN ---------------------------------
  write( *, * )
  write( *, * ) 'STARTING NLTE MODE'
  write( *, * )

  ! Eddington approximation works well as a first guess for the VEFs.
  vef_f_k( :, : ) = 1.0d+0 / 3.0d+0
  vef_f_h( :, : ) = 1.0d+0 / sqrt( 3.0d+0 )

  !call write_vefs

  ! keep previous results so we can test for convergence
  vef_f_k_old( :, : ) = vef_f_k( :, : )
  vef_f_h_old( :, : ) = vef_f_h( :, : )

  ! As a first guess for S_NLTE, use J_LTE. (This guess is waaaay better than
  ! S_NLTE = B.)
  call calc_source_fn
  !call write_source_fn

  ! TODO: iterate this loop until VEFs converge
  do i3 = 1, 20

    write( *, '(a20, 2x, i3)') 'NLTE ITERATION #: ', i3

    ! solve scattering problem to get J
    call solve_scatt_prob

    ! use J to calculate S(J)
    call calc_source_fn

    ! use S to calculate little_j (Feautrier variable)
    call solve_rte

    ! use little_j to calculate 2nd moment K
    call calc_2nd_moment_k

    ! use new values of J and K to get new values of f_K and f_H
    call calc_vefs
    !call write_vefs
    write(*, '(a20, 1x, es11.3e2)') 'f_K RMS change: ', &
    calc_rmsd( vef_f_k( :, 1 ), vef_f_k_old( :, 1 ) )
    vef_f_k_old( :, : ) = vef_f_k( :, : )
    vef_f_h_old( :, : ) = vef_f_h( :, : )

  end do
!-------------------------------END NLTE RUN------------------------------------

  ! use little_j to calculate flux
  call calc_flux
  call write_moments

  call write_spectrum

  close( 22 )
  close( 23 )
  close( 24 )
  close( 25 )
  close( 30 )

  stop
end program main
