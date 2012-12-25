! Calculates the difference in optical depth between points d+1 and d.
function calc_dtau( d )
  use precision_mod
  use global
  implicit none

  real(kind=dp) :: calc_dtau
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
    calc_dtau = tau_grid( d + 1 ) - tau_grid( d )
  end if

end function calc_dtau


function dtau_mu( depth_idx, mu_idx )
  use precision_mod
  use global
  implicit none

  real(kind=dp) :: dtau_mu
  integer :: depth_idx, mu_idx

  dtau_mu = dtau( depth_idx ) / abs( mu_grid( mu_idx ) )

end function dtau_mu



! Average change in optical depth between points d and d-1, and points d+1 and d.
function calc_dtau_tilde( d )
  use precision_mod
  use interfaces, only: calc_dtau
  use global
  implicit none

  real(kind=dp) :: calc_dtau_tilde
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
    calc_dtau_tilde = ( dtau( d - 1 ) + dtau( d ) ) / 2.0d+0
  end if

end function calc_dtau_tilde
