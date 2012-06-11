module interfaces
  interface


  function dtau( d )
    use precision_mod
    implicit none
    real(kind=dp) :: dtau
    integer :: d
  end function dtau


  function dtau_tilde( d )
    use precision_mod
    implicit none
    real(kind=dp) :: dtau_tilde
    integer :: d
  end function dtau_tilde


  function dtau_mu( depth_idx, mu_idx )
    use precision_mod
    implicit none
    real(kind=dp) :: dtau_mu
    integer :: depth_idx, mu_idx
  end function dtau_mu


  function planck_fn( lambda, temp )
    use precision_mod
    implicit none
    real(kind=dp) :: planck_fn
    real(kind=dp) :: lambda, temp
  end function planck_fn


  subroutine gauleg(x1,x2,x,w,n)
    use precision_mod
    implicit none
    integer, intent(in) :: n
    real(kind=dp), intent(in) :: x1,x2
    real(kind=dp), intent(out) :: x(:),w(:)
  end subroutine gauleg

  subroutine twerp(x0,y0,x,y,npts)                                  
    use precision_mod
    implicit none                                                     
    real(kind=dp), dimension(:), intent(in) :: x0,y0
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(out) :: y
    integer, intent(in) :: npts
  end subroutine twerp

  subroutine locate(xarray,n,x,j)
    use precision_mod
    integer, intent(in) :: n
    real(kind=dp), intent(in) :: x, xarray(n)
    integer, intent(out) :: j
  end subroutine locate

  subroutine stop_exit(status,message)
    integer :: status
    character(len=*) :: message
  end subroutine stop_exit

  function interp_i( depth_idx, mu, wl_idx )
    use precision_mod
    implicit none
    real(kind=dp) :: interp_i
    real(kind=dp) :: mu
    integer :: depth_idx, wl_idx
  end function interp_i
  
  end interface

end module interfaces
