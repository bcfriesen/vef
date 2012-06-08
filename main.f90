! A code for solving the radiative transfer equation in a plane-parallel
! atmosphere using variable Eddington factors. See article "Difference
! Equations and Linearization Methods" by Lawrence Auer
! in "Methods in Radiative Transfer" ed. Wolfgang Kalkofen (1984).

program main
  use interfaces
  use global
  implicit none

  integer :: i1, i2, i3

  ! The optical depth grid will be logarithmic, such that
  ! tau_i = tau_(i-1) * f, 
  ! where f is this constant, to be calculated below.
  real :: tau_grid_ratio

  ! Need to set first two tau points by hand because the logarithmic spacing
  ! won't do much good multiplying zero by some number over and over again.
  tau_grid( 1 ) = 0.0
  tau_grid( 2 ) = tau_min
  tau_grid_ratio = ( tau_max / tau_grid( 2 ) )**(1.0 / ( n_depth_pts - 2 ) )
  ! Set up optical depth grid.
  do i1 = 3, n_depth_pts
    tau_grid( i1 ) = tau_grid( i1 - 1 ) * tau_grid_ratio
  end do

  ! Set up direction cosine grid.
  do i1 = 1, n_mu_pts
    mu_grid( i1 ) = -1.0 + real( i1 - 1 ) * ( 2.0 / real( n_mu_pts - 1 ) )
  end do

  ! Set up wavelength grid.
  do i1 = 1, n_wl_pts
    ! TODO: add wavelength dependence later
    wl_grid( i1 ) = 5000.0
  end do

  ! Set up boundary conditions.
  i_lambda( :, :, : ) = 0.0
  do i1 = 1, n_wl_pts
    do i2 = 1, n_mu_pts
      ! At the surface, only downward rays have a boundary condition.
      if ( mu_grid ( i2 ) < 0.0 ) then
        i_lambda( 1, i2, i1 ) = 0.0
      ! At depth, only upward rays have a boundary condition.
      else
        i_lambda( n_depth_pts, i2, i1 ) = planck_fn( wl_grid( i1 ), temp )
      end if
    end do
  end do

!----------------------START LTE RUN TO GET FIRST GUESS FOR J-------------------
  write( *, * ) 'starting LTE run...'
  write( *, * ) 'setting S = B...'
  do i1 = 1, n_depth_pts
    do i2 = 1, n_wl_pts
      source_fn( i1, i2 ) = planck_fn( wl_grid( i2 ), temp )
    end do
  end do
  write( *, * ) 'finished setting S = B!'

  write( *, * ) 'calculating formal solution...'
  call formal_soln
  write( *, * ) 'formal solution calculated successfully!'

  ! In LTE we don't need J to find I since we can write down the formal solution
  ! immediately. However we need a not-completely-terrible guess for J (i.e.,
  ! better than J = B) when we start solving the scattering equation in NLTE. So
  ! we calculate J here according to its definition.
  write( *, * ) 'calulating J, H, and K...'
  call calc_moments
  write( *, * ) 'finished calulating J, H, and K!'

  write( *, * ) 'finished LTE run!'
  write( *, * )
!--------------------------------END LTE RUN------------------------------------



!-------------------------------START NLTE RUN ---------------------------------
  write( *, * ) 'starting NLTE run...'
  ! Eddington approximation works well as a first guess for the VEFs.
  write( *, * ) 'applying Eddington approx. for f_H and f_K...'
  vef_f_k( :, : ) = 1.0 / 3.0
  vef_f_h( :, : ) = 1.0 / sqrt( 3.0 )
  write( *, * ) 'finished applying Eddington approx. for f_H and f_K!'

  vef_f_k_old( :, : ) = vef_f_k
  vef_f_h_old( :, : ) = vef_f_h

  ! iterate this loop until VEFs converge
  do i3 = 1, 3

    write( *, * ) 'NLTE ITERATION #: ', i3

    do i1 = 1, n_depth_pts
      write( 11, * ) tau_grid( i1 ), vef_f_k( i1, 1 )
    end do

    ! calculate NLTE source function, given J (if this the first iteration in
    ! the NLTE loop then J will be J_LTE)
    write( *, * ) 'calculating NLTE source function...'
    call calc_source_fn
    write( *, * ) 'finished calculating NLTE source function!'

    write( *, * ) 'calculating formal solution...'
    call formal_soln
    write( *, * ) 'finished calculating formal solution!'

    write( *, * ) 'calculating J, H, and K...'
    call calc_moments
    write( *, * ) 'finished calculating J, H, and K!'

    write( *, * ) 'calculating Feautrier variables...'
    call calc_feautrier_vars
    write( *, * ) 'finished calculating Feautrier variables!'

    write( *, * ) 'calculating new VEFs...'
    call calc_vefs
    write( *, * ) 'finished calculating new VEFs!'

    ! solve scattering problem to find J
    write( *, * ) 'solving scattering problem...'
    call solve_scatt_prob
    write( *, * ) 'finished solving scattering problem!'

  end do
!-------------------------------END NLTE RUN------------------------------------


  stop
end program main



! Calculates the difference in optical depth between points d+1 and d.
function dtau( d )
  use global
  implicit none

  real :: dtau
  integer :: d

  if ( d < 1 ) then
    write( *, * ) 'ERROR: tau grid index must be at least 1'
    write( *, * ) 'requested index for dtau: ', d
    stop
  else if ( d > ( size( tau_grid ) - 1 ) ) then
    write( *, * ) 'ERROR: dtau is hanging off the edge of the grid!'
    write( *, * ) 'requested index for dtau: ', d
    stop
  else
    dtau = tau_grid( d + 1 ) - tau_grid( d )
  end if

end function dtau


function dtau_mu( depth_idx, mu_idx )
  use global
  use interfaces, only: dtau
  implicit none

  real :: dtau_mu
  integer :: depth_idx, mu_idx

  dtau_mu = dtau( depth_idx ) / abs( mu_grid( mu_idx ) )

end function dtau_mu



! Average change in optical depth between points d and d-1, and points d+1 and d.
function dtau_tilde( d )
  use interfaces, only: dtau
  use global
  implicit none

  real :: dtau_tilde
  integer :: d

  if ( d < 1 ) then
    write( *, * ) 'ERROR: tau grid index must be at least 1'
    write( *, * ) 'requested index for dtau: ', d
    stop
  else if ( d > ( size( tau_grid ) - 1 ) ) then
    write( *, * ) 'ERROR: dtau is hanging off the edge of the grid!'
    write( *, * ) 'requested index for dtau: ', d
    stop
  else
    dtau_tilde = ( dtau( d - 1 ) + dtau( d ) ) / 2.0
  end if

end function dtau_tilde

! Interpolate I in mu-space so we can integrate to find J.
function interp_i( depth_idx, mu, wl_idx )
  use global
  use interfaces, only: twerp
  implicit none
  real :: interp_i
  real :: mu
  integer :: depth_idx, wl_idx

  call twerp( mu_grid( : ), i_lambda( depth_idx, : , wl_idx ), &
              mu, interp_i, n_mu_pts )
end function interp_i
