subroutine write_vefs
  use global
  implicit none

  integer :: i1, i2
  character(len=20) :: wfmt

  write( 24, * )
  write( 25, * )
  wfmt = '(a15)'
  write( 24, wfmt ) '# NEW ITERATION'
  write( 25, wfmt ) '# NEW ITERATION'
  write( 24, * )
  write( 25, * )
  wfmt = '(2a12)'
  write( 24, wfmt ) 'tau', 'f_K'
  write( 25, wfmt ) 'tau', 'f_H'
  wfmt = '(2es12.4e2)'
  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      write( 24, wfmt ) tau_grid( i2 ), vef_f_k( i2, i1 )
    end do
    write( 25, wfmt ) tau_grid( 1 ), vef_f_h( i1 )
  end do

end subroutine write_vefs
