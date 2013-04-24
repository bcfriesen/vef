subroutine calc_vefs
  use global
  use interfaces, only: stop_exit
  implicit none

  integer :: i1, i2
  character(len=*), parameter :: whoami = 'calc_vefs'

  vef_f_k( :, : ) = 0.0d+0
  vef_f_h( :, : ) = 0.0d+0

  do i1 = 1, n_wl_pts

    ! calculate f_K at all interior points
    do i2 = 1, n_depth_pts
      vef_f_k( i2, i1 ) = k_lambda( i2, i1 ) / j_lambda( i2, i1 )
      if ( vef_f_k( i2, i1 ) < 0.0d+0 ) call stop_exit( 1, whoami, &
      'f_K < 0' )
    end do

    ! calculate f_H at surface
    ! trapezoid rule
    vef_f_h( 1, i1 ) = little_j( 1, 1, i1 ) * mu_grid( 1 ) + &
    little_j( 1, n_mu_pts, i1 ) * mu_grid( n_mu_pts )
    do i2 = 2, n_mu_pts - 1
      vef_f_h( 1, i1 ) = vef_f_h( 1, i1 ) + &
      2.0d+0 * little_j( 1, i2, i1 ) * mu_grid( i2 )
    end do
    vef_f_h( 1, i1 ) = vef_f_h( 1, i1 ) * ( mu_grid( n_mu_pts ) - &
    mu_grid( 1 ) ) / ( 2.0d+0 * real( n_mu_pts - 1 ) )

    vef_f_h( 1, i1 ) = vef_f_h( 1, i1 ) / j_lambda( 1, i1 )

    ! TODO: not sure what to do with this one...
    ! calculate f_H at depth
    ! trapezoid rule
    vef_f_h( 2, i1 ) = little_j( n_depth_pts, 1, i1 ) * mu_grid( 1 ) + &
    little_j( n_depth_pts, n_mu_pts, i1 ) * mu_grid( n_mu_pts )
    do i2 = 2, n_mu_pts - 1
      vef_f_h( 2, i1 ) = vef_f_h( 2, i1 ) + &
      2.0d+0 * little_j( n_depth_pts, i2, i1 ) * mu_grid( i2 )
    end do
    vef_f_h( 2, i1 ) = vef_f_h( 2, i1 ) * ( mu_grid( n_mu_pts ) - &
    mu_grid( 1 ) ) / ( 2.0d+0 * real( n_mu_pts - 1 ) )

    vef_f_h( 2, i1 ) = vef_f_h( 2, i1 ) / j_lambda( n_depth_pts, i1 )

  end do

end subroutine calc_vefs
