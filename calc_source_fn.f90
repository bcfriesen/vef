! Given the mean intensity J, calculate the source function. In this case we
! use the simple Milne-Eddington source function which is linear in J.
subroutine calc_source_fn
  use interfaces, only: planck_fn
  use global
  implicit none

  integer :: i1, i2
  do i1 = 1, n_depth_pts
    do i2 = 1, n_wl_pts
      ! Milne-Eddington source function
      source_fn( i1, i2 ) = me_therm_parm * &
      planck_fn( wl_grid( n_wl_pts ), temp ) * ( 1.0 - me_therm_parm ) * &
      j_lambda( i1, i2 )
    end do
  end do
end subroutine calc_source_fn
