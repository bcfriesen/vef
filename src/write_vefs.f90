subroutine write_vefs
  use global
  implicit none

  integer :: i1, i2
  character(len=*), parameter :: wfmt = '(2es12.4e2)'

  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      write( 24, wfmt ) tau_grid( i2 ), vef_f_k( i2, i1 )
    end do
    write( 25, wfmt ) tau_grid( 1 ), vef_f_h( i1 )
  end do

end subroutine write_vefs