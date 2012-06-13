subroutine calc_2nd_moment_k
  use precision_mod
  use global
  use interfaces, only: stop_exit
  implicit none

  integer :: i1, i2, i3
  character(len=*), parameter :: whoami = 'calc_2nd_moment_k'

  k_lambda( :, : ) = 0.0d+0

  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      do i3 = 1, n_mu_pts - 1
        ! Romberg integration
        k_lambda( i2, i1 ) = k_lambda( i2, i1 ) + &
        0.5d+0 * ( little_j( i2, i3, i1 ) * mu_grid( i3 )**2 + &
        little_j( i2, i3 + 1, i1 ) * mu_grid( i3 + 1 )**2 ) * &
        ( mu_grid( i3 + 1 ) - mu_grid( i3 ) )
      end do

      ! even moments of I can't be negative, but odd moments can be
      if ( k_lambda( i2, i1 ) < 0.0d+0 ) then
        call stop_exit( 1, whoami, 'K < 0' )
      end if

    end do
  end do

end subroutine calc_2nd_moment_k
