subroutine write_moments
  use global
  use interfaces, only: planck_fn
  implicit none

  integer :: i1, i2
  character(len=20) :: wfmt
  
  write( 23, * )
  wfmt = '(a15)'
  write( 23, wfmt ) '# NEW ITERATION'
  write( 23, * )
  wfmt = '(5a12)'
  write( 23, wfmt ) 'tau', 'J', 'H', 'K', 'B'
  wfmt = '(5es12.4e2)'
  do i1 = 1, n_wl_pts
    do i2 = 1, n_depth_pts
      write( 23, wfmt ) tau_grid( i2 ), j_lambda( i2, i1 ), &
                        h_lambda( i2, i1 ), k_lambda( i2, i1 ), &
                        planck_fn( wl_grid( i1 ) * a2cm, temp )
    end do
  end do

end subroutine write_moments
