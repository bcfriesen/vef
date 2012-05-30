! a code for solving the radiative transfer equation in a plane-parallel
! atmosphere using variable Eddington factors (see article "Difference
! Equations and Linearization Methods" by Lawrence Auer
! in "Methods in Radiative Transfer" ed. Wolfgang Kalkofen (1984)

module interfaces
  interface

  function dtau( d )
    implicit none

    real :: dtau
    integer :: d
  end function dtau

  function dtau_tilde( d )
    implicit none

    real :: dtau_tilde
    integer :: d
  end function dtau_tilde
  
  end interface

end module interfaces



module global

  ! number of optical depth points
  integer, parameter :: n_depth_pts = 100
  ! number of direction cosine points
  integer, parameter :: n_mu_pts = 10

  ! maximum optical depth to consider
  real, parameter :: tau_max = 1.0e+4
  ! minimum non-zero optical depth to consider
  real, parameter :: tau_min = 1.0e-8

  ! thermalization parameter for source function in isotropic, monochromatic
  ! scattering (Milne-Eddington problem). eps = 1 means pure LTE; eps = 0 means
  ! pure scattering (like SYNOW)
  real, parameter :: eps = 1.0e-4
  ! optical depth grid
  real, dimension( n_depth_pts ) :: tau_grid
  ! Eddington factor f_K = K / J
  real, dimension( n_depth_pts, n_mu_pts ) :: f_k
  ! Eddington factor f_H = \int_0^1 j(\mu) \mu d\mu / J
  real :: f_h
  ! source function
  real, dimension( n_depth_pts, n_mu_pts ) :: s
  ! mean intensity
  real, dimension( n_depth_pts, n_mu_pts ) :: j

end module global




program main
  use interfaces
  use global
  implicit none

  integer :: i1
  ! the optical depth grid will be logarithmic, such that tau_i = tau_(i-1) * f, 
  ! where f is this constant, to be calculated below
  real :: tau_grid_ratio

  ! need to set first two tau points by hand because the logarithmic spacing
  ! won't do much good multiplying zero by some number over and over again
  tau_grid( 1 ) = 0.0
  tau_grid( 2 ) = tau_min
  tau_grid_ratio = ( tau_max / tau_grid( 2 ) )**(1.0 / ( n_depth_pts - 2 ) )
  ! set up optical depth grid
  do i1 = 3, n_depth_pts
    tau_grid( i1 ) = tau_grid( i1 - 1 ) * tau_grid_ratio
  end do

  stop
end program main



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
