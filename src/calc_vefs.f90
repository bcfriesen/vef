subroutine calc_vefs
  use global
  use interfaces, only: stop_exit
  implicit none

  integer :: i1, i2
  character(len=*), parameter :: whoami = 'calc_vefs'

  vef_f_k( :, : ) = 0.0d+0
  vef_f_h( : ) = 0.0d+0

  do i1 = 1, n_wl_pts

    ! calculate f_K at all interior points
    do i2 = 1, n_depth_pts
      vef_f_k( i2, i1 ) = k_lambda( i2, i1 ) / j_lambda( i2, i1 )
      if ( vef_f_k( i2, i1 ) < 0.0d+0 ) call stop_exit( 1, whoami, &
      'f_K < 0' )
    end do

    ! calculate f_H at surface
    do i2 = 1, n_mu_pts - 1
      ! Romberg integration from 0 < mu < 1
      vef_f_h( i1 ) = vef_f_h( i1 ) + 0.5d+0 * ( little_j( 1, i2, i1 ) &
      * mu_grid( i2 ) + little_j( 1, i2 + 1, i1 ) * mu_grid( i2 + 1 ) ) * &
      ( mu_grid( i2 + 1 ) - mu_grid( i2 ) )
    end do
    vef_f_h( i1 ) = vef_f_h( i1 ) / j_lambda( 1, i1 )
    if ( vef_f_h( i1 ) < 0.0d+0 ) call stop_exit( 1, whoami, 'f_H < 0' )

  end do

end subroutine calc_vefs
