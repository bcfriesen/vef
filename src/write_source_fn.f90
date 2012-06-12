subroutine write_source_fn
  use global
  implicit none

  integer :: i1, i2
  character(len=*), parameter :: wfmt = '(2es12.4e2)'

  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      write( 22, wfmt ) tau_grid( i2 ), source_fn( i2, i1 )
    end do
  end do

end subroutine write_source_fn
