! Solve scattering problem, i.e., find mean intensity so that we can calculate
! source function. In LTE we don't need to do this at all because we already
! known S = B. So this is just for NLTE.
subroutine solve_scatt_prob
  use precision_mod
  use global
  use interfaces, only: dtau, dtau_tilde, planck_fn, stop_exit
  implicit none

  ! NOTE: this method assumes static media, i.e., wavelengths aren't coupled.
  ! If they are we have to solve the coupled scattering problems simultaneously.

  integer :: i1, i2

  ! Set up machinery for matrix equation. (matrix * unknown_vector = rhs_vector)
  ! TODO: change from generic matrix solver to tridiagonal solver
  real(kind=dp), dimension( n_depth_pts, n_depth_pts ) :: matrix
  real(kind=dp), dimension( n_depth_pts ) :: rhs
  ! pivot indices; LAPACK needs these
  integer, dimension( n_depth_pts ) :: ipiv
  ! LAPACK driver info flag
  integer :: info

  ! Since wavelengths aren't coupled we can solve the scattering problem for
  ! each wavelength separately. This outer loop would parallelize perfectly.
  do i1 = 1, n_wl_pts

    matrix( :, : ) = 0.0d+0
    rhs( : ) = 0.0d+0

    ! boundary conditions at surface
    ! TODO: integrate the BCs for I to get these BCs rather than hard-coding
    ! them here.
    matrix( 1, 1 ) = ( -vef_f_k( 1, i1 ) / dtau( 1 ) ) - vef_f_h( 1, i1 )
    matrix( 1, 2 ) = vef_f_k( 2, i1 ) / dtau( 1 )
    ! No external illumination means this term is zero.
    rhs( 1 ) = 0.0d+0

    ! boundary conditions at depth
    matrix( n_depth_pts, n_depth_pts - 1 ) = -vef_f_k( n_depth_pts - 1, i1 ) / &
    dtau( n_depth_pts - 1 )
    matrix( n_depth_pts, n_depth_pts ) = ( vef_f_k( n_depth_pts, i1 ) / &
    dtau( n_depth_pts - 1 ) ) + vef_f_h( 2, i1 )
    ! Blackbody at depth. Factor of 1/2 comes from integration of I_+ to get
    ! H_+.
    ! TODO: integrate the BCs for I to get these BCs rather than hard-coding
    ! them here.
    rhs ( n_depth_pts ) = planck_fn( wl_grid( i1 ), temp ) / 2.0d+0

    do i2 = 2, n_depth_pts - 1
      matrix( i2, i2 - 1 ) = ( 1.0d+0 / ( dtau( i2 - 1 ) * dtau_tilde( i2 ) ) ) * vef_f_k( i2 - 1, i1 )
      matrix( i2, i2 ) = ( -2.0d+0 / ( dtau( i2 - 1 ) * dtau( i2 ) ) ) * vef_f_k( i2, i1 ) - 1.0d+0
      matrix( i2, i2 + 1 ) = ( 1.0d+0 / ( dtau( i2 ) * dtau_tilde( i2 ) ) ) * vef_f_k( i2 + 1, i1 )
      rhs( i2 ) = -source_fn( i2, i1 )
    end do

    ! Solve for mean intensity.
    call dgesv( n_depth_pts, 1, matrix, n_depth_pts, ipiv, rhs, &
    n_depth_pts, info)

    if ( info /= 0 ) then
      write( *, * ) 'DGESV returned info = ', info
      call stop_exit(1, 'could not solve scattering problem')
    end if

    j_lambda( :, i1 ) = rhs( : )

    do i2 = 1, n_depth_pts
      if ( j_lambda( i2, i1 ) < 0.0d+0 ) then
        write( *, * ) i2, j_lambda( i2, i1 )
        call stop_exit( 1, 'J < 0' )
      end if
    end do

  end do

end subroutine solve_scatt_prob