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


  subroutine stop_exit(stat, who, message)
    implicit none
    integer :: stat
    character(len=*) :: who, message
  end subroutine stop_exit


  Function calc_rmsd(array1, array2)
    Use precision_mod
    Implicit None
    Real (Kind=dp) :: calc_rmsd
    Real (Kind=dp), Dimension (:) :: array1, array2
  End Function calc_rmsd

  function gaussian(lambda, norm, lambda_0, sigma)
    use precision_mod
    implicit none
    real(kind=dp) :: gaussian
    real(kind=dp) :: lambda, norm, lambda_0, sigma
  end function gaussian

  end interface

end module interfaces
