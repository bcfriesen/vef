subroutine calc_vefs
  use global
  implicit none

  integer :: i1, i2

  vef_f_k( :, : ) = 0.0
  vef_f_h( :, : ) = 0.0

  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      vef_f_k( i2, i1 ) = k_lambda( i2, i1 ) / j_lambda( i2, i1 )
    end do

    ! calculate f_H at surface
    do i2 = 1, n_mu_pts - 1
      ! Romberg integration from 0 < mu < 1
      if ( mu_grid( i2 ) > 0.0 ) then
        vef_f_h( 1, i1 ) = vef_f_h( 1, i1 ) + 0.5 * ( little_j( 1, i2, i1 ) &
        * mu_grid( i2 ) + little_j( 1, i2 + 1, i1 ) * mu_grid( i2 + 1 ) ) * &
        ( mu_grid( i2 + 1 ) - mu_grid( i2 ) )
      end if
    end do
    vef_f_h( 1, i1 ) = vef_f_h( 1, i1 ) / j_lambda( 1, i1 )
    write( *, * ) 'vef_f_h(1) = ', vef_f_h( 1, i1 )

    ! calculate f_H at depth
    do i2 = 1, n_mu_pts - 1
      ! Romberg integration
      if ( mu_grid( i2 ) > 0.0 ) then
        vef_f_h( 2, i1 ) = vef_f_h( 2, i1 ) + 0.5 * &
        ( little_j( n_depth_pts, i2, i1 ) * mu_grid( i2 ) + &
        little_j( n_depth_pts, i2 + 1, i1 ) * mu_grid ( i2 + 1 ) ) * &
        ( mu_grid( i2 + 1 ) - mu_grid( i2 ) )
      end if
    end do
    ! This index notation is dumb, with 1 being at the surface (as usual) but 2
    ! being at depth. Maybe I should just make vef_f_h_surf and vef_f_h_depth
    ! separated variables...
    vef_f_h( 2, i1 ) = vef_f_h( 2, i1 ) / j_lambda( n_depth_pts, i1 )
    write( *, * ) 'vef_f_h(2) = ', vef_f_h( 2, i1 )

  end do

end subroutine calc_vefs
