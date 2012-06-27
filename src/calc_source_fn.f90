! Given the mean intensity J, calculate the source function. In this case we
! use the simple Milne-Eddington source function which is linear in J.
subroutine calc_source_fn
  use interfaces, only: planck_fn, stop_exit, phi=>gaussian
  use global
  use const
  implicit none

  integer :: i1, i2, i3
  character(len=*), parameter :: whoami = 'calc_source_fn'
  real(kind=dp) :: j_bar
  ! rest wavelength of the fake line (Angstrom)
  real(kind=dp) :: lambda_0
  ! broadening parameter of line (Angstrom)
  real(kind=dp) :: sigma
  ! # of integration points
  integer, parameter :: n_int_pts = 30
  ! normalization constant on line profile
  real(kind=dp) :: gauss_norm


  do i1 = 1, n_wl_pts

    do i2 = 1, n_depth_pts

      if ( j_lambda( i2, i1 ) < 0.0d+0 ) &
        call stop_exit( 1, whoami, 'J < 0' )

      lambda_0 = 3000.d0
      sigma = 50.d0
      gauss_norm = 1.d0 / (sigma * sqrt(2.d0*pi))
      ! integrate J_lambda over profile function to get J_bar (independent of
      ! lambda)
      j_bar = 0.d0
      ! trapezoid integration
      j_bar = j_lambda(i2, 1) * phi(wl_grid(1), gauss_norm, lambda_0, sigma) + &
              j_lambda(i2, n_wl_pts) * phi(wl_grid(n_wl_pts), gauss_norm, &
              lambda_0, sigma)
      do i3 = 2, n_wl_pts - 1
        j_bar = j_bar + 2.d0 * j_lambda(i2, i3) * &
                phi(wl_grid(i3), gauss_norm, lambda_0, sigma)
      end do
      j_bar = j_bar * (wl_grid(n_wl_pts) - wl_grid(1)) / &
              (2.d0 * (n_wl_pts-1))

      ! Milne-Eddington source function
      source_fn( i2, i1 ) = &
      me_therm_parm * planck_fn( wl_grid( i1 ) * a2cm, temp ) + &
      ( 1.0d+0 - me_therm_parm ) * j_bar
      
      if ( source_fn( i2, i1 ) < 0.0d+0 ) &
        call stop_exit( 1, whoami, 'S < 0' )

    end do

  end do
end subroutine calc_source_fn
