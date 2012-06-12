subroutine solve_rte
  use precision_mod
  use global
  use interfaces, only: dtau, dtau_tilde, planck_fn, stop_exit
  implicit none

  ! NOTE: this method assumes static media, i.e., wavelengths aren't coupled.
  ! If they are we have to solve the coupled scattering problems simultaneously.

  integer :: i1, i2, i3
  character(len=*), parameter :: whoami = 'solve_rte'

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

    ! loop over direction cosines
    do i2 = 1, n_mu_pts

      matrix( :, : ) = 0.0d+0
      rhs( : ) = 0.0d+0

      ! boundary conditions at surface
      matrix( 1, 1 ) = -( 1.0d+0 / dtau( 1 ) ) - ( 1.0d+0 / mu_grid( i2 ) )
      matrix( 1, 2 ) = 1.0d+0 / dtau( 1 )
      ! No external illumination means this term is zero.
      rhs( 1 ) = 0.0d+0

      ! boundary conditions at depth
      matrix( n_depth_pts, n_depth_pts - 1 ) = -1.0d+0 / &
                                               dtau( n_depth_pts - 1 )
      matrix( n_depth_pts, n_depth_pts ) = 1.0d+0 / dtau( n_depth_pts - 1 ) + &
                                           1.0d+0 / mu_grid( i2 )
      ! Blackbody at depth. (isothermal -> no dB/dtau term)
      rhs ( n_depth_pts ) = planck_fn( wl_grid( i1 ), temp ) / &
                            mu_grid( i2 )

      do i3 = 2, n_depth_pts - 1
        matrix( i3, i3 - 1 ) = 1.0d+0 / ( dtau( i3 - 1 ) * dtau_tilde( i3 ) )
        matrix( i3, i3 ) = -2.0d+0 / ( dtau( i3 - 1 ) * dtau( i3 ) ) - &
                           1.0d+0 / mu_grid( i2 )**2
        matrix( i3, i3 + 1 ) = 1.0d+0 / ( dtau( i3 ) * dtau_tilde( i3 ) )
        rhs( i3 ) = -1.0d+0 / mu_grid( i2 )**2 * source_fn( i3, i1 )
      end do

      ! Solve for Feautrier variable j
      call dgesv( n_depth_pts, 1, matrix, n_depth_pts, ipiv, rhs, &
                  n_depth_pts, info)

      if ( info /= 0 ) then
        write( *, * ) 'DGESV returned info = ', info
        call stop_exit(1, whoami, 'could not solve scattering problem')
      end if

      little_j( :, i2, i1 ) = rhs( : )

      do i3 = 1, n_depth_pts
        if (little_j( i3, i2, i1 ) < 0.0d+0 ) then
          write( *, * ) i3, i2, little_j( i3, i2, i1 )
          call stop_exit( 1, whoami, 'little_j < 0' )
        end if
      end do

    end do

  end do

end subroutine solve_rte
