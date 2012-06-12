! Given the mean intensity J, calculate the source function. In this case we
! use the simple Milne-Eddington source function which is linear in J.
subroutine calc_source_fn
  use interfaces, only: planck_fn, stop_exit
  use global
  implicit none

  integer :: i1, i2
  character(len=*), parameter :: whoami = 'calc_source_fn'

  do i1 = 1, n_wl_pts

    do i2 = 1, n_depth_pts

      if ( j_lambda( i2, i1 ) < 0.0d+0 ) &
        call stop_exit( 1, whoami, 'J < 0' )

      ! Milne-Eddington source function
      source_fn( i2, i1 ) = &
      me_therm_parm * planck_fn( wl_grid( i1 ), temp ) + &
      ( 1.0d+0 - me_therm_parm ) * j_lambda( i2, i1 )
      
      if ( source_fn( i2, i1 ) < 0.0d+0 ) &
        call stop_exit( 1, whoami, 'S < 0' )

    end do

  end do
end subroutine calc_source_fn
