subroutine formal_soln
  use global
  use interfaces, only: dtau_mu
  implicit none

  real, dimension( 2 : n_depth_pts, n_mu_pts ) :: e0i, e1i, e2i

  ! linear interpolation coefficients for constructing downward rays
  real, dimension( 2 : n_depth_pts, n_mu_pts ) :: amil, bmil, gmil
  ! parabolic interpolation coefficients for constructing downward rays
  real, dimension( 2 : n_depth_pts - 1, n_mu_pts ) :: amip, bmip, gmip
  real, dimension( 2 : n_depth_pts, n_mu_pts, n_wl_pts ) :: delta_i_m

  ! linear interpolation coefficients for constructing upward rays
  real, dimension( 1 : n_depth_pts - 1, n_mu_pts) :: apil, bpil, gpil
  ! parabolic interpolation coefficients for constructing upward rays
  real, dimension( 2 : n_depth_pts - 1, n_mu_pts) :: apip, bpip, gpip
  real, dimension( 1 : n_depth_pts - 1, n_mu_pts, n_wl_pts ) :: delta_i_p

  integer :: i1, i2, i3

  ! Once S is known everywhere, we use the formal solution to the RTE to
  ! calculate I everywhere.

  do i3 = 2, n_depth_pts
    do i2 = 1, n_mu_pts
      e0i( i3, i2 ) = 1.0 - exp( -dtau_mu( i3 - 1, i2 ) )
      e1i( i3, i2 ) = dtau_mu( i3 - 1, i2 ) - e0i( i3, i2 )
      e2i( i3, i2 ) = dtau_mu( i3 - 1, i2 )**2 - 2.0 * e1i( i3, i2 )
    end do
  end do

  ! These are parabolic interpolation coefficients. See Olson & Kunasz (1987),
  ! Eqs. 17 - 21.
  ! TODO: fix these! For some reason I get negative values for intensity when I
  ! use parabolic interpolation. I get reasonable answers for linear
  ! interpolation though. Something ain't right here.
  do i3 = 2, n_depth_pts - 1
    do i2 = 1, n_mu_pts
      amip( i3, i2 ) = e0i( i3, i2 ) + ( e2i( i3, i2 ) &
                          - ( dtau_mu( i3, i2 ) &
                          + 2.0 * dtau_mu( i3 - 1, i2 ) &
                          * e1i( i3, i2 ) ) ) &
                          / ( dtau_mu( i3 - 1, i2 ) * ( dtau_mu( i3, i2 ) &
                          + dtau_mu( i3 - 1, i2 ) ) )
      bmip( i3, i2 ) = ( ( dtau_mu( i3, i2 ) + dtau_mu( i3 - 1, i2 ) ) &
                          * e1i( i3, i2 ) - e2i( i3, i2 ) ) &
                          / ( dtau_mu( i3 - 1, i2 ) * dtau_mu( i3, i2 ) )
      gmip( i3, i2 ) = ( e2i( i3, i2 ) - dtau_mu( i3 - 1, i2 ) &
                          * e1i( i3, i2 ) ) / ( dtau_mu( i3, i2 ) &
                          * ( dtau_mu( i3, i2 ) + dtau_mu( i3 - 1, i2 ) ) )
      apip( i3, i2 ) = ( e2i( i3 + 1, i2 ) - dtau_mu( i3, i2 ) &
                          * e1i( i3 + 1, i2 ) ) &
                          / ( dtau_mu( i3 - 1, i2 ) &
                          * ( dtau_mu( i3, i2 ) + dtau_mu( i3 - 1, i2 ) ) )
      bpip( i3, i2 ) = ( ( dtau_mu( i3, i2 ) + dtau_mu( i3 - 1, i2 ) ) &
                          * e1i( i3 + 1, i2 ) - e2i( i3 + 1, i2 ) ) &
                          / ( dtau_mu( i3 - 1, i2 ) * dtau_mu( i3, i2 ) )
      gpip( i3, i2 ) = e0i( i3 + 1, i2 ) + ( e2i( i3 + 1, i2 ) &
                          - ( dtau_mu( i3 - 1, i2 ) &
                          + 2.0 * dtau_mu( i3, i2 ) ) &
                          * e1i( i3 + 1, i2 ) ) / ( dtau_mu( i3, i2 ) &
                          * ( dtau_mu( i3, i2 ) + dtau_mu( i3 - 1, i2 ) ) )
    end do
  end do

! These are linear interpolation coefficients. At the boundaries we always
! interpolate linearly because parabolic interpolation coefficients hang off the
! edge of the grid.
  do i2 = 1, n_mu_pts
    do i3 = 2, n_depth_pts
      amil( i3, i2 ) = e0i( i3, i2 ) - e1i( i3, i2 ) / dtau_mu( i3 - 1, i2 )
      bmil( i3, i2 ) = e1i( i3, i2 ) / dtau_mu( i3 - 1, i2 )
      gmil( i3, i2 ) = 0.0
    end do
    do i3 = 1, n_depth_pts - 1
      apil( i3, i2 ) = 0.0
      bpil( i3, i2 ) = e1i( i3 + 1, i2 ) / dtau_mu( i3, i2 )
      gpil( i3, i2 ) = e0i( i3 + 1, i2 ) - e1i( i3 + 1, i2 ) / dtau_mu( i3, i2 )
    end do
  end do

  do i1 = 1, n_wl_pts
    do i2 = 1, n_mu_pts
      if ( mu_grid( i2 ) > 0.0 ) then
        ! For all interior points we can interpolate parabolically. At
        ! boundaries we interpolate linearly.
        delta_i_p( 1, i2, i1 ) = bpil( 1, i2 ) * source_fn ( 1, i1 ) &
        + gpil( 1, i2) * source_fn( 2, i1 )
        do i3 = 2, n_depth_pts - 1
          delta_i_p( i3, i2, i1 ) = apil( i3, i2 ) * source_fn( i3 - 1, i1 ) &
                                  + bpil( i3, i2 ) * source_fn( i3    , i1 ) &
                                  + gpil( i3, i2 ) * source_fn( i3 + 1, i1 )
        end do
      else
        do i3 = 2, n_depth_pts - 1
        delta_i_m( i3, i2, i1 ) = amil( i3, i2 ) * source_fn( i3 - 1, i1 ) &
                                + bmil( i3, i2 ) * source_fn( i3    , i1 ) &
                                + gmil( i3, i2 ) * source_fn( i3 + 1, i1 )
        end do
        delta_i_m( n_depth_pts, i2, i1 ) = amil( n_depth_pts, i2 ) * &
        source_fn( n_depth_pts - 1, i1 ) + bmil( n_depth_pts, i2 ) * &
        source_fn( n_depth_pts, i1 )
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
          if ( i_lambda( i3, i2, i1 ) < 0.0 ) &
            call stop_exit( 1, 'ERROR: I_lambda < 0!' )
        end do
      else
        do i3 = n_depth_pts - 1, 1, -1
          i_lambda( i3, i2, i1 ) = i_lambda( i3 + 1, i2, i1 ) &
                                   * exp( -dtau_mu( i3, i2 ) ) &
                                   + delta_i_p( i3, i2, i1 )
          if ( i_lambda( i3, i2, i1 ) < 0.0 ) &
            call stop_exit( 1, 'ERROR: I_lambda < 0!' )
        end do
      end if
    end do
  end do

end subroutine formal_soln
