subroutine write_spectrum
  use global
  implicit none

  integer :: i1
  character(len=80) :: wfmt

  wfmt = '(2a12)'
  write( 30, wfmt ) 'lambda (A)', 'H'
  wfmt = '(3es12.4e2)'
  do i1 = 1, n_wl_pts
    write( 30, wfmt ) wl_grid( i1 ), h_lambda( 1, i1 )
  end do

end subroutine write_spectrum
