subroutine calc_0th_moment_j
  use global
  implicit none

  character(len=*), parameter :: whoami = 'calc_0th_moment_j'
  integer :: i1, i2, i3

  j_lambda( :, : ) = 0.0d+0

  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      do i3 = 1, n_mu_pts - 1
        ! Romberg integration
        j_lambda( i2, i1 ) = j_lambda( i2, i1 ) + &
        0.5d+0 * ( little_j( i2, i3, i1 ) + &
        little_j( i2, i3 + 1, i1 ) ) * &
        ( mu_grid( i3 + 1 ) - mu_grid( i3 ) )
      end do

      ! even moments of I can't be negative, but odd moments can be
      if ( j_lambda( i2, i1 ) < 0.0d+0 ) then
        call stop_exit( 1, whoami, 'J < 0' )
      end if

    end do
  end do

end subroutine calc_0th_moment_j
