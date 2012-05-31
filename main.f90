! A code for solving the radiative transfer equation in a plane-parallel
! atmosphere using variable Eddington factors (see article "Difference
! Equations and Linearization Methods" by Lawrence Auer
! in "Methods in Radiative Transfer" ed. Wolfgang Kalkofen (1984).

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

  function planck_fn( lambda, temp )
    implicit none

    real :: planck_fn
    real :: lambda, temp
  end function planck_fn
  
  end interface

end module interfaces



module global

  ! # of optical depth points
  integer, parameter :: n_depth_pts = 100
  ! # of direction cosine points
  integer, parameter :: n_mu_pts = 11
  ! # of wavelength points
  integer, parameter :: n_wl_pts = 1

  ! maximum optical depth to consider
  real, parameter :: tau_max = 1.0e+4
  ! minimum non-zero optical depth to consider
  real, parameter :: tau_min = 1.0e-8

  ! Thermalization parameter for source function in isotropic, monochromatic
  ! scattering (Milne-Eddington problem). eps = 1 means pure LTE; eps = 0 means
  ! pure scattering (like SYNOW).
  real, parameter :: me_therm_parm = 1.0e-4
  ! optical depth grid
  real, dimension( n_depth_pts ) :: tau_grid
  ! direction cosine grid
  real, dimension( n_mu_pts ) :: mu_grid
  ! wavelength grid
  real, dimension( n_wl_pts ) :: wl_grid
  ! Eddington factor f_K = K / J
  real, dimension( n_depth_pts, n_wl_pts ) :: vef_f_k
  ! Eddington factor f_H = \int_0^1 j(\mu) \mu d\mu / J
  real, dimension( 2, n_wl_pts ) :: vef_f_h
  ! source function
  real, dimension( n_depth_pts, n_wl_pts ) :: source_fn
  ! mean intensity
  real, dimension( n_depth_pts, n_wl_pts ) :: j_lambda
  ! specific intensity
  real, dimension( n_depth_pts, n_mu_pts, n_wl_pts ) :: i_lambda

  ! blackbody temperature at depth
  ! TODO: add temperature dependence later
  real, parameter :: temp = 1.0e4

end module global




program main
  use interfaces
  use global
  implicit none

  integer :: i1, i2, i3

  ! The optical depth grid will be logarithmic, such that tau_i = tau_(i-1) * f, 
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
    mu_grid( i1 ) = -1.0 + real( i1 ) * ( 2.0 / real( n_mu_pts ) )
  end do

  ! Set up wavelength grid.
  do i1 = 1, n_wl_pts
    ! TODO: add wavelength dependence later
    wl_grid( i1 ) = 5000.0
  end do

  ! LTE loop first
  do i1 = 1, n_depth_pts
    do i2 = 1, n_wl_pts
      source_fn( i1, i2 ) = planck_fn( wl_grid( i2 ), temp )
    end do
  end do

  ! Eddington approximation works well as a first guess for the VEFs.
  vef_f_k( :, : ) = 1.0 / 3.0
  vef_f_h( :, : ) = 1.0 / sqrt( 3.0 )

  ! solve scattering problem
  call solve_scatt_prob

  do i1 = 1, n_wl_pts
    !write( *, * ) 'lambda = ', wl_grid( i1 )
      do i2 = 1, n_depth_pts
        write( *, * ) tau_grid( i2 ), j_lambda( i2, i1 )
      end do
  end do


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


! Planck function
function planck_fn( lambda, temp )
  implicit none

  real :: planck_fn
  real :: lambda, temp
  real, parameter :: h_planck = 6.626e-27, &
                     c_light  = 2.998e+10, &
                     k_boltz  = 1.380e-16

  planck_fn = ( 2.0 * h_planck * c_light**2 / lambda**5 ) / &
              ( exp( ( h_planck * c_light / lambda ) / k_boltz * temp ) - 1.0 )
end function planck_fn



! Given the mean intensity J, calculate the source function. In this case we
! use the simple Milne-Eddington source function which is linear in J.
subroutine calc_source_fn
  use interfaces, only: planck_fn
  use global
  implicit none

  integer :: i1, i2
  do i1 = 1, n_depth_pts
    do i2 = 1, n_wl_pts
      ! Milne-Eddington source function
      source_fn( i1, i2 ) = me_therm_parm * &
      planck_fn( wl_grid( n_wl_pts ), temp ) * ( 1.0 - me_therm_parm ) * &
      j_lambda( i1, i2 )
    end do
  end do
end subroutine calc_source_fn



! Solve scattering problem, i.e., find mean intensity so that we can calculate
! source function. In LTE we don't need to do this at all because we already
! known S = B. So this is just for NLTE.
subroutine solve_scatt_prob
  use global
  use interfaces, only: dtau, dtau_tilde, planck_fn
  implicit none

  ! NOTE: this method assumes static media, i.e., wavelengths aren't coupled.
  ! If they are we have to solve the coupled scattering problems simultaneously.

  integer :: i1, i2

  ! Set up machinery for matrix equation. (matrix * unknowns = rhs)
  ! TODO: change from generic matrix solver to tridiagonal solver
  real, dimension( n_depth_pts, n_depth_pts ) :: matrix
  real, dimension( n_depth_pts ) :: rhs
  ! pivot indices; LAPACK needs these
  integer, dimension( n_depth_pts ) :: ipiv
  ! LAPACK driver info flag
  integer :: info

  ! Since wavelengths aren't coupled we can solve the scattering problem for
  ! each wavelength separately. This outer loop would parallelize perfectly.
  do i1 = 1, n_wl_pts

    ! boundary conditions at surface
    matrix( 1, 1 ) = ( -vef_f_k( 1, i1 ) / dtau( 1 ) ) - vef_f_h( 1, i1 )
    matrix( 1, 2 ) = vef_f_k( 2, i1 ) / dtau( 1 )
    ! No external illumination means this term is zero.
    rhs( 1 ) = 0.0

    ! boundary conditions at depth
    matrix( n_depth_pts, n_depth_pts - 1 ) = -vef_f_k( n_depth_pts - 1, i1 ) / &
    dtau( n_depth_pts - 1 )
    matrix( n_depth_pts, n_depth_pts ) = ( vef_f_k( n_depth_pts, i1 ) / &
    dtau( n_depth_pts - 1 ) ) + vef_f_h( 2, i1 )
    ! Blackbody at depth. Factor of 1/2 comes from integration of I_+ to get H_+.
    rhs ( n_depth_pts ) = planck_fn( wl_grid( i1 ), temp ) / 2.0

    do i2 = 2, n_depth_pts - 1
      matrix( i2, i2 - 1 ) = ( 2.0 / ( dtau( i2 - 1 ) * ( dtau( i2 - 1 ) + &
      dtau( i2 ) ) ) ) * vef_f_k( i2 - 1, i1 )
      matrix( i2, i2 ) = ( -( 2.0 / ( dtau( i2 - 1 ) * &
      dtau( i2 ) ) ) * vef_f_k( i2, i1 ) ) - 1.0
      matrix( i2, i2 + 1 ) = ( 2.0 / ( dtau( i2 ) * ( dtau( i2 - 1 ) + &
      dtau( i2 ) ) ) ) * vef_f_k( i2 + 1, i1 )
      rhs( i2 ) = -source_fn( i2, i1 )
    end do

    ! Solve for mean intensity.
    call sgesv( n_depth_pts, 1, matrix, n_depth_pts, ipiv, rhs, &
    n_depth_pts, info)

    j_lambda( :, i1 ) = rhs( : )

  end do

end subroutine solve_scatt_prob


! Solves the radiative transfer equation along each ray. Returns the specific
! intensity I_lambda as a function of wavelength, direction cosine, and optical
! depth.
subroutine solve_rte
  use global
  use interfaces, only: dtau, planck_fn
  implicit none
  
  real, dimension( n_depth_pts, n_depth_pts ) :: matrix
  real, dimension( n_depth_pts ) :: rhs
  integer, dimension( n_depth_pts ) :: ipiv

  integer :: i1, i2, i3
  integer :: info

  ! create first-order derivative stencil
  do i1 = 1, n_wl_pts
    do i2 = 1, n_mu_pts
      write( *, * ) 'solving RTE for wl = ', wl_grid( i1 ), &
      'and mu = ', mu_grid( i2 )
      do i3 = 2, n_depth_pts - 1
        matrix( i3, i3 - 1 ) = -1.0 / ( 2.0 * dtau( i3 - 1 ) )
        matrix( i3, i3 ) = ( ( dtau( i3 ) - dtau( i3 - 1 ) ) / &
        ( 2.0 * dtau( i3 ) * dtau( i3 - 1 ) ) ) - ( 1.0 / mu_grid( i2 ) )
        matrix( i3, i3 + 1 ) = 1.0 / ( 2.0 * dtau( i3 ) )
        rhs( i3 ) = -source_fn( i3, i1 ) / mu_grid( i2 )
      end do
      ! The boundary conditions for the RTE depend on the direction cosine of
      ! the ray.
      if ( mu_grid( i2 ) > 0.0 ) then
        matrix( 1, 1 ) = -( 1.0 / dtau( 1 ) ) - ( 1.0 / mu_grid( i2 ) )
        matrix( 1, 2 ) = 1.0 / dtau ( 1 )
        matrix( n_depth_pts, n_depth_pts ) = 1.0
        rhs( 1 ) = -source_fn( 1, i1 ) / mu_grid( i2 )
        rhs( n_depth_pts ) = planck_fn( wl_grid( i1 ), temp )
      else
        ! In principle this number can be anything because this BC is just
        ! I_-(tau = 0) = 0, which is the same as 5.0 * I_-(tau = 0) = 0.
        matrix( 1, 1 ) = 1.0
        matrix( n_depth_pts, n_depth_pts - 1 ) = - 1.0 / dtau( n_depth_pts - 1 )
        matrix( n_depth_pts, n_depth_pts ) = &
        ( 1.0 / dtau( n_depth_pts - 1 ) ) - ( 1.0 / mu_grid( i2 ) )
        rhs( 1 ) = 0.0
        rhs( n_depth_pts ) = -source_fn( n_depth_pts, i1 ) / mu_grid( i2 )
      end if
      call sgesv( n_depth_pts, 1, matrix, n_depth_pts, ipiv, rhs, &
                  n_depth_pts, info)
      write( *, * ) 'INFO = ', info
      i_lambda( :, i2, i1 ) = rhs( : )
    end do
  end do

end subroutine solve_rte
