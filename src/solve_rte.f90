subroutine solve_rte
  use precision_mod
  use global
  use interfaces, only: dtau, dtau_tilde, planck_fn, stop_exit
  implicit none

  ! NOTE: this method assumes static media, i.e., wavelengths aren't coupled.
  ! If they are we have to solve the coupled scattering problems simultaneously.

  integer :: i1, i2, i3
  character(len=*), parameter :: whoami = 'solve_rte'
  ! error message
  character(len=80) :: errmsg
  ! error message format
  character(len=80) :: wfmt

  ! Set up machinery for matrix equation. (matrix * unknown_vector = rhs_vector)
  ! Since the matrix is tridiagonal we only need 3 vectors to hold all the
  ! elements.
  real(kind=dp), dimension( n_depth_pts ) :: diag
  real(kind=dp), dimension( n_depth_pts - 1 ) :: ldiag, udiag
  real(kind=dp), dimension( n_depth_pts ) :: rhs
  ! LAPACK driver info flag
  integer :: info

  ! Since wavelengths aren't coupled we can solve the scattering problem for
  ! each wavelength separately. This outer loop would parallelize perfectly.
  do i1 = 1, n_wl_pts

    ! loop over direction cosines
    do i2 = 1, n_mu_pts

      diag( : ) = 0.0d+0
      udiag( : ) = 0.0d+0
      ldiag( : ) = 0.0d+0
      rhs( : ) = 0.0d+0

      ! boundary conditions at surface
      diag( 1 ) = -( 1.0d+0 / dtau( 1 ) ) - ( 1.0d+0 / mu_grid( i2 ) )
      udiag( 1 ) = 1.0d+0 / dtau( 1 )
      ! No external illumination means this term is zero.
      ! TODO: integrate the BCs for I to get these BCs rather than hard-coding
      ! them here.
      rhs( 1 ) = 0.0d+0

      ! boundary conditions at depth
      ldiag( n_depth_pts - 1 ) = -1.0d+0 / dtau( n_depth_pts - 1 )
      diag( n_depth_pts ) = 1.0d+0 / dtau( n_depth_pts - 1 ) + &
                            1.0d+0 / mu_grid( i2 )
      ! Blackbody at depth. (isothermal -> no dB/dtau term)
      ! TODO: integrate the BCs for I to get these BCs rather than hard-coding
      ! them here.
      rhs ( n_depth_pts ) = planck_fn( wl_grid( i1 ), temp ) / &
                            mu_grid( i2 )

      do i3 = 2, n_depth_pts - 1
        ldiag( i3 - 1 ) = 1.0d+0 / ( dtau( i3 - 1 ) * dtau_tilde( i3 ) )
        diag( i3 ) = -2.0d+0 / ( dtau( i3 - 1 ) * dtau( i3 ) ) - &
                     1.0d+0 / mu_grid( i2 )**2
        udiag( i3 ) = 1.0d+0 / ( dtau( i3 ) * dtau_tilde( i3 ) )
        rhs( i3 ) = -1.0d+0 / mu_grid( i2 )**2 * source_fn( i3, i1 )
      end do

      ! Solve for Feautrier variable j
      call dgtsv( n_depth_pts, 1, ldiag, diag, udiag, rhs, n_depth_pts, info)

      if ( info /= 0 ) then
        if ( info < 0 ) then
          wfmt = '(a36, 2x, i2 )'
          write(errmsg, wfmt) 'this argument to DGTSV was invalid: ', -info
        else if ( info > 0 ) then
          wfmt = '(a50)'
          write(errmsg, wfmt) 'DGTSV upper triangular matrix factor U is &
                              &singular'
        end if
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
