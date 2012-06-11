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


  subroutine stop_exit(status,message)
    integer :: status
    character(len=*) :: message
  end subroutine stop_exit


  end interface

end module interfaces
