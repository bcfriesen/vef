subroutine calc_moments
  use global
  use interfaces, only: interp_i, gauleg
  implicit none

  integer :: i1, i2, i3
  integer, parameter :: n_int_pts = 20
  real, dimension( n_int_pts ) :: int_pts, int_wghts

  j_lambda( :, : ) = 0.0
  h_lambda( :, : ) = 0.0
  k_lambda( :, : ) = 0.0
  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      call gauleg( -1.0, 1.0, int_pts, int_wghts, n_int_pts )
      do i3 = 1, n_int_pts
        j_lambda( i2, i1 ) = j_lambda( i2, i1 ) &
        + interp_i( i2, int_pts(i3), i1 ) * int_wghts( i3 )
        h_lambda( i2, i1 ) = h_lambda( i2, i1 ) &
        + interp_i( i2, int_pts(i3), i1 ) * int_pts( i3 ) * int_wghts( i3 )
        k_lambda( i2, i1 ) = k_lambda( i2, i1 ) &
        + interp_i( i2, int_pts(i3), i1 ) * int_pts( i3 )**2 * int_wghts( i3 )
      end do

      ! even moments of I can't be negative
      if ( j_lambda( i2, i1 ) < 0.0 ) then
        call stop_exit( 1, 'ERROR: J < 0!' )
      else if ( k_lambda( i2, i1 ) < 0.0 ) then
        call stop_exit( 1, 'ERROR: K < 0!' )
      end if

      j_lambda( i2, i1 ) = ( 1.0 / 2.0 ) * j_lambda( i2, i1 )
      h_lambda( i2, i1 ) = ( 1.0 / 2.0 ) * h_lambda( i2, i1 )
      k_lambda( i2, i1 ) = ( 1.0 / 2.0 ) * k_lambda( i2, i1 )

    end do
  end do

end subroutine calc_moments
