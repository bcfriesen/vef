subroutine write_vefs
  use global
  implicit none

  integer :: i1, i2

  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      write( 24, * ) tau_grid( i2 ), vef_f_k( i2, i1 )
    end do
    write( 25, * ) tau_grid( 1 ), vef_f_h( 1, i1 )
    write( 25, * ) tau_grid( n_depth_pts ), vef_f_h( 2, i1 )
  end do

end subroutine write_vefs
