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

  function c_ij(i, j)
    use precision_mod
    implicit none
    real(kind=dp) :: c_ij
    integer :: i, j
  end function c_ij

  function c_ji(j, i)
    use precision_mod
    implicit none
    real(kind=dp) :: c_ji
    integer :: j, i
  end function c_ji

  function c_ik(i)
    use precision_mod
    implicit none
    real(kind=dp) :: c_ik
    integer :: i
  end function c_ik

  function c_ki(i)
    use precision_mod
    implicit none
    real(kind=dp) :: c_ki
    integer :: i
  end function c_ki

  function r_ij(i, j)
    use precision_mod
    implicit none
    real(kind=dp) :: r_ij
    integer :: i, j
  end function r_ij

  function z_ji(j, i)
    use precision_mod
    implicit none
    real(kind=dp) :: z_ji
    integer :: j, i
  end function z_ji

  function r_ik(i)
    use precision_mod
    implicit none
    real(kind=dp) :: r_ik
    integer :: i
  end function r_ik

  function r_ki(i)
    use precision_mod
    implicit none
    real(kind=dp) :: r_ki
    integer :: i
  end function r_ki

    FUNCTION f_ij(i, j, n_e, t)
      USE precision_mod
      IMPLICIT NONE
      REAL (KIND=dp) :: f_ij
      INTEGER :: i
      INTEGER :: j
      REAL (KIND=dp) :: n_e
      REAL (KIND=dp) :: t
    END FUNCTION f_ij

    FUNCTION part(i, j)
      USE precision_mod
      IMPLICIT NONE
      REAL (KIND=dp) :: part
      INTEGER :: i
      INTEGER :: j
    END FUNCTION part

    FUNCTION saha(i, j, n_e, t)
      USE precision_mod
      IMPLICIT NONE
      REAL (KIND=dp) :: saha
      INTEGER :: i
      INTEGER :: j
      REAL (KIND=dp) :: n_e
      REAL (KIND=dp) :: t
    END FUNCTION saha

    FUNCTION solvene(n_e)
      USE precision_mod
      IMPLICIT NONE
      REAL (KIND=dp) :: solvene
      REAL (KIND=dp) :: n_e
    END FUNCTION solvene

    double precision function zeroin(ax,bx,f,tol)
      double precision ax,bx,f,tol
      external f
    end function zeroin

  end interface

end module interfaces
