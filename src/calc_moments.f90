subroutine calc_0th_moment_j
  use global
  implicit none

  character(len=*), parameter :: whoami = 'calc_0th_moment_j'
  integer :: i1, i2, i3

  j_lambda( :, : ) = 0.0d+0

  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      ! trapezoid rule
      j_lambda( i2, i1 ) = little_j( i2, 1, i1 ) + little_j( i2, n_mu_pts, i1 )
      do i3 = 2, n_mu_pts - 1
        j_lambda( i2, i1 ) = j_lambda( i2, i1 ) + &
        2.0d+0 * little_j( i2, i3, i1 )
      end do
      j_lambda( i2, i1 ) = j_lambda( i2, i1 ) * ( mu_grid( n_mu_pts ) - &
      mu_grid( 1 ) ) / ( 2.0d+0 * real( n_mu_pts - 1 ) )

      ! even moments of I can't be negative, but odd moments can be
      if ( j_lambda( i2, i1 ) < 0.0d+0 ) then
        call stop_exit( 1, whoami, 'J < 0' )
      end if

    end do
  end do

end subroutine calc_0th_moment_j


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
      ! trapezoid rule
      k_lambda( i2, i1 ) = little_j( i2, 1, i1 ) * mu_grid( 1 )**2 + &
      little_j( i2, n_mu_pts, i1 ) * mu_grid( n_mu_pts )**2
      do i3 = 2, n_mu_pts - 1
        k_lambda( i2, i1 ) = k_lambda( i2, i1 ) + &
        2.0d+0 * little_j( i2, i3, i1 ) * mu_grid( i3 )**2
      end do
      k_lambda( i2, i1 ) = k_lambda( i2, i1 ) * ( mu_grid( n_mu_pts ) - &
      mu_grid( 1 ) ) / ( 2.0d+0 * real( n_mu_pts - 1 ) )

      ! even moments of I can't be negative, but odd moments can be
      if ( k_lambda( i2, i1 ) < 0.0d+0 ) then
        call stop_exit( 1, whoami, 'K < 0' )
      end if

    end do
  end do

end subroutine calc_2nd_moment_k
