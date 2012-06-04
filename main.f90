! A code for solving the radiative transfer equation in a plane-parallel
! atmosphere using variable Eddington factors. See article "Difference
! Equations and Linearization Methods" by Lawrence Auer
! in "Methods in Radiative Transfer" ed. Wolfgang Kalkofen (1984).

program main
  use interfaces
  use global
  implicit none

  integer :: i1, i2, i3
  integer, parameter :: n_int_pts = 20
  real, dimension( n_int_pts ) :: int_pts, int_wghts

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
    mu_grid( i1 ) = -1.0 + real( i1 ) * ( 2.0 / real( n_mu_pts ) )
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
  do i1 = 1, n_depth_pts
    do i2 = 1, n_wl_pts
      source_fn( i1, i2 ) = planck_fn( wl_grid( i2 ), temp )
    end do
  end do

  call formal_soln

  !do i1 = 1, n_wl_pts
  !  do i2 = 1, n_depth_pts
  !    write( *, * ) tau_grid( i2 ), i_lambda( i2, 1, i1 )
  !  end do
  !end do

  ! In LTE we don't need J to find I since we can write down the formal solution
  ! immediately. However we need a not-completely-terrible guess for J (i.e.,
  ! better than J = B) when we start solving the scattering equation in NLTE. So
  ! we calculate J here according to its definition.
  j_lambda( :, : ) = 0.0
  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      call gauleg( minval( mu_grid( : ) ), maxval( mu_grid( : ) ), int_pts, &
                   int_wghts, n_int_pts )
      do i3 = 1, n_int_pts
        j_lambda( i2, i1 ) = j_lambda( i2, i1 ) &
        + interp_i( i2, int_pts(i3), i1 ) * int_wghts( i3 )
      end do
      j_lambda( i2, i1 ) = ( 1.0 / 2.0 ) * j_lambda( i2, i1 )
    end do
  end do

  do i1 = 1, n_depth_pts
    write( *, * ) tau_grid( i1 ), j_lambda( i1, 1 )
  end do
!--------------------------------END LTE RUN------------------------------------



!-------------------------------START NLTE RUN ---------------------------------
  ! Eddington approximation works well as a first guess for the VEFs.
  vef_f_k( :, : ) = 1.0 / 3.0
  vef_f_h( :, : ) = 1.0 / sqrt( 3.0 )

  vef_f_k_old( :, : ) = vef_f_k
  vef_f_h_old( :, : ) = vef_f_h

  ! use LTE value of J to calculate first guess at NLTE source function
  call calc_source_fn

  ! solve scattering problem to find J
  call solve_scatt_prob

  do i3 = 1, 2
    ! calculate NLTE source function, given J (if this the first iteration in
    ! the NLTE loop then J will be J_LTE)
    call calc_source_fn

    !do i1 = 1, n_wl_pts
    !  do i2 = 1, n_depth_pts
    !    write( *, * ) tau_grid( i2 ), j_lambda( i2, i1 )
    !  end do
    !end do

    call solve_scatt_prob

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
    ! Blackbody at depth. Factor of 1/2 comes from integration of I_+ to get
    ! H_+.
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


subroutine formal_soln
  use global
  use interfaces, only: dtau_mu
  implicit none

  real, dimension( 2 : n_depth_pts, n_mu_pts ) :: e0i, e1i, e2i

  real, dimension( 2 : n_depth_pts, n_mu_pts ) :: ami, bmi, gmi
  real, dimension( 2 : n_depth_pts, n_mu_pts, n_wl_pts ) :: delta_i_m

  real, dimension( 1 : n_depth_pts - 1, n_mu_pts) :: api, bpi, gpi
  real, dimension( 1 : n_depth_pts - 1, n_mu_pts, n_wl_pts ) :: delta_i_p

  integer :: i1, i2, i3

  ! Once S is known everywhere, we use the formal solution to the RTE to
  ! calculate I everywhere.

  do i1 = 2, n_depth_pts
    do i2 = 1, n_mu_pts
      e0i( i1, i2 ) = 1.0 - exp( -dtau_mu( i1 - 1, i2 ) )
      e1i( i1, i2 ) = dtau_mu( i1 - 1, i2 ) - e0i( i1, i2 )
      e2i( i1, i2 ) = dtau_mu( i1 - 1, i2 )**2.0 - 2.0 * e1i( i1, i2 )
    end do
  end do

  ! These are parabolic interpolation coefficients. See Olson & Kunasz (1987),
  ! Eqs. 17 - 21.
! do i1 = 2, n_depth_pts - 1
!   do i2 = 1, n_mu_pts
!     ami( i1, i2 ) = e0i( i1, i2 ) + ( e2i( i1, i2 ) &
!                         - ( dtau_mu( i1, i2 ) &
!                         + 2.0 * dtau_mu( i1 - 1, i2 ) &
!                         * e1i( i1, i2 ) ) ) &
!                         / ( dtau_mu( i1 - 1, i2 ) * ( dtau_mu( i1, i2 ) &
!                         + dtau_mu( i1 - 1, i2 ) ) )
!     bmi( i1, i2 ) = ( ( dtau_mu( i1, i2 ) + dtau_mu( i1 - 1, i2 ) ) &
!                         * e1i( i1, i2 ) - e2i( i1, i2 ) ) &
!                         / ( dtau_mu( i1 - 1, i2 ) * dtau_mu( i1, i2 ) )
!     gmi( i1, i2 ) = ( e2i( i1, i2 ) - dtau_mu( i1 - 1, i2 ) &
!                         * e1i( i1, i2 ) ) / ( dtau_mu( i1, i2 ) &
!                         * ( dtau_mu( i1, i2 ) + dtau_mu( i1 - 1, i2 ) ) )
!     api( i1, i2 ) = ( e2i( i1 + 1, i2 ) - dtau_mu( i1, i2 ) &
!                         * e1i( i1 + 1, i2 ) ) &
!                         / ( dtau_mu( i1 - 1, i2 ) &
!                         * ( dtau_mu( i1, i2 ) +  dtau_mu( i1 - 1, i2 ) ) )
!     bpi( i1, i2 ) = ( ( dtau_mu( i1, i2 ) + dtau_mu( i1 - 1, i2 ) ) &
!                         * e1i( i1 + 1, i2 ) - e2i( i1 + 1, i2 ) ) &
!                         / ( dtau_mu( i1 - 1, i2 ) * dtau_mu( i1, i2 ) )
!     gpi( i1, i2 ) = e0i( i1 + 1, i2 ) + ( e2i( i1 + 1, i2 ) &
!                         - ( dtau_mu( i1 - 1, i2 ) &
!                         + 2.0 * dtau_mu( i1, i2 ) ) &
!                         * e1i( i1 + 1, i2 ) ) / ( dtau_mu( i1, i2 ) &
!                         * ( dtau_mu( i1, i2 ) + dtau_mu( i1 - 1, i2 ) ) )
!   end do
! end do

! These are linear interpolation coefficients.
  do i2 = 1, n_mu_pts
    do i1 = 2, n_depth_pts
      ami( i1, i2 ) = e0i( i1, i2 ) - e1i( i1, i2 ) / dtau_mu( i1 - 1, i2 )
      bmi( i1, i2 ) = e1i( i1, i2 ) / dtau_mu( i1 - 1, i2 )
      gmi( i1, i2 ) = 0.0
    end do
    do i1 = 1, n_depth_pts - 1
      api( i1, i2 ) = 0.0
      bpi( i1, i2 ) = e1i( i1 + 1, i2 ) / dtau_mu( i1, i2 )
      gpi( i1, i2 ) = e0i( i1 + 1, i2 ) - e1i( i1 + 1, i2 ) / dtau_mu( i1, i2 )
    end do
  end do

  do i1 = 1, n_wl_pts
    do i2 = 1, n_mu_pts
      if ( mu_grid( i2 ) > 0.0 ) then
        do i3 = 1, n_depth_pts - 1
          delta_i_p( i3, i2, i1 ) = api( i3, i2 ) * source_fn( i3 - 1, i1 ) &
                                  + bpi( i3, i2 ) * source_fn( i3    , i1 ) &
                                  + gpi( i3, i2 ) * source_fn( i3 + 1, i1 )
        end do
      else
        do i3 = 2, n_depth_pts
        delta_i_m( i3, i2, i1 ) = ami( i3, i2 ) * source_fn( i3 - 1, i1 ) &
                                + bmi( i3, i2 ) * source_fn( i3    , i1 ) &
                                + gmi( i3, i2 ) * source_fn( i3 + 1, i1 )
        end do
      end if
    end do
  end do

  ! Interpolate to calculate I everywhere.
  do i1 = 1, n_wl_pts
    do i2 = 1, n_mu_pts
      if ( mu_grid( i2 ) < 0.0 ) then
        do i3 = 2, n_depth_pts
          i_lambda( i3, i2, i1 ) = i_lambda( i3 - 1, i2, i1 ) &
                                   * exp( -dtau_mu( i3 - 1, i2 ) ) &
                                   + delta_i_m( i3, i2, i1 )
        end do
      else
        do i3 = n_depth_pts - 1, 1, -1
          i_lambda( i3, i2, i1 ) = i_lambda( i3 + 1, i2, i1 ) &
                                   * exp( -dtau_mu( i3, i2 ) ) &
                                   + delta_i_p( i3, i2, i1 )
        end do
      end if
    end do
  end do

end subroutine formal_soln


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
