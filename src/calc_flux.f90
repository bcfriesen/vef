subroutine calc_flux
  use global
  use interfaces, only: dtau_mu
  implicit none

  integer :: i1, i2, i3

  do i1 = 1, n_wl_pts
    do i2 = 1, n_mu_pts
      ! one-sided, 1st order derivative at surface
      little_h( 1, i2, i1 ) = ( little_j( 2, i2, i1 ) - &
                                little_j( 1, i2, i1 ) ) / &
                                dtau_mu( 1, i2 )
      ! centered, 2nd order deriative in interior
      do i3 = 2, n_depth_pts - 1
        little_h( i3, i2, i1 ) = ( little_j( i3 + 1, i2, i1 ) - &
                                   little_j( i3 - 1, i2, i1 ) ) / &
                                   ( 2.0d+0 * dtau_mu( i3, i2 ) )
      end do
      ! one_sided, 1st order derivative at depth
      little_h( n_depth_pts, i2, i1 ) = ( little_j( n_depth_pts, i2, i1 ) - &
                                      little_j( n_depth_pts - 1, i2, i1 ) ) / &
                                      dtau_mu( n_depth_pts - 1, i2 )
    end do
  end do

end subroutine calc_flux
