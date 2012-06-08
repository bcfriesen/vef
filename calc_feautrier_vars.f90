subroutine calc_feautrier_vars
  use global
  implicit none

  integer :: i1, i2, i3

  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      do i3 = 1, n_mu_pts
        little_j( i2, i3, i1 ) = ( 1.0 / 2.0 ) * &
        ( i_lambda( i2, i3, i1 ) + i_lambda( i2, n_mu_pts - ( i3 - 1 ), i1 ) )
        little_h( i2, i3, i1 ) = ( 1.0 / 2.0 ) * &
        ( i_lambda( i2, i3, i1 ) - i_lambda( i2, n_mu_pts - ( i3 - 1 ), i1 ) )
      end do
    end do
  end do
end subroutine calc_feautrier_vars
