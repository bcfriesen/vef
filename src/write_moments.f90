subroutine write_moments
  use global
  implicit none

  integer :: i1, i2
  character(len=*), parameter :: wfmt = '(4es12.4e2)'
  
  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      write( 23, wfmt ) tau_grid( i2 ), j_lambda( i2, i1 ), h_lambda( i2, i1 ), &
      k_lambda( i2, i1 )
    end do
  end do

end subroutine write_moments
