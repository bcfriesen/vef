module interfaces
  interface

  function dtau( d )
    implicit none

    real :: dtau
    integer :: d
  end function dtau

  function dtau_tilde( d )
    implicit none

    real :: dtau_tilde
    integer :: d
  end function dtau_tilde

  function dtau_mu( depth_idx, mu_idx )
    implicit none

    real :: dtau_mu
    integer :: depth_idx, mu_idx
  end function dtau_mu

  function planck_fn( lambda, temp )
    implicit none

    real :: planck_fn
    real :: lambda, temp
  end function planck_fn

  subroutine gauleg(x1,x2,x,w,n)
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: x1,x2
  real, intent(out) :: x(:),w(:)
  end subroutine gauleg

  subroutine twerp(x0,y0,x,y,npts)                                  
      implicit none                                                     
      integer npts                                                      
      real :: x,y                                                     
      real :: x0,y0                                                   
      dimension x0(*),y0(*)                                             
  end subroutine twerp

      subroutine locate(xarray,n,x,j)
!     -------------------------------
      integer :: j,n
      real :: x,xarray(n)
      end subroutine locate

      subroutine stop_exit(status,message)
!     ------------------------------------
      integer status
      character*(*) message
      end subroutine stop_exit

function interp_i( depth_idx, mu, wl_idx )
  implicit none
  real :: interp_i
  real :: mu
  integer :: depth_idx, wl_idx
end function interp_i
  
  end interface

end module interfaces
