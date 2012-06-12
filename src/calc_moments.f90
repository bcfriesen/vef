subroutine calc_moments
  use precision_mod
  use global
  use interfaces, only: stop_exit
  implicit none

  integer :: i1, i2, i3
  character(len=*), parameter :: whoami = 'calc_moments'

  j_lambda( :, : ) = 0.0d+0
  h_lambda( :, : ) = 0.0d+0
  k_lambda( :, : ) = 0.0d+0

  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      do i3 = 1, n_mu_pts - 1

        j_lambda( i2, i1 ) = j_lambda( i2, i1 ) + &
        0.5d+0 * ( little_j( i2, i3, i1 ) + little_j( i2, i3 + 1, i1 ) ) * &
        ( mu_grid( i3 + 1 ) - mu_grid( i3 ) )

        h_lambda( i2, i1 ) = h_lambda( i2, i1 ) + &
        0.5d+0 * ( little_h( i2, i3, i1 ) * mu_grid( i3 ) + &
        little_h( i2, i3 + 1, i1 ) * mu_grid( i3 + 1 ) ) * &
        ( mu_grid( i3 + 1 ) - mu_grid( i3 ) )

        k_lambda( i2, i1 ) = k_lambda( i2, i1 ) + &
        0.5d+0 * ( little_j( i2, i3, i1 ) * mu_grid( i3 )**2 + &
        little_j( i2, i3 + 1, i1 ) * mu_grid( i3 + 1 )**2 ) * &
        ( mu_grid( i3 + 1 ) - mu_grid( i3 ) )

      end do

      ! even moments of I (J and K) can't be negative, but odd moments (H) can
      if ( j_lambda( i2, i1 ) < 0.0d+0 ) then
        call stop_exit( 1, whoami, 'J < 0' )
      else if ( k_lambda( i2, i1 ) < 0.0d+0 ) then
        call stop_exit( 1, whoami, 'K < 0' )
      end if

    end do
  end do

end subroutine calc_moments
