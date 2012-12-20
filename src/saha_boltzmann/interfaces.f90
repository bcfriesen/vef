MODULE INTERFACES
  IMPLICIT NONE

  INTERFACE

    FUNCTION f_ij(i, j, n_e, t)
      USE machine
      IMPLICIT NONE
      REAL (KIND=dp) :: f_ij
      INTEGER :: i
      INTEGER :: j
      REAL (KIND=dp) :: n_e
      REAL (KIND=dp) :: t
    END FUNCTION f_ij

    FUNCTION part(i, j)
      USE machine
      IMPLICIT NONE
      REAL (KIND=dp) :: part
      INTEGER :: i
      INTEGER :: j
    END FUNCTION part

    FUNCTION saha(i, j, n_e, t)
      USE machine
      IMPLICIT NONE
      REAL (KIND=dp) :: saha
      INTEGER :: i
      INTEGER :: j
      REAL (KIND=dp) :: n_e
      REAL (KIND=dp) :: t
    END FUNCTION saha

    FUNCTION solvene(n_e)
      USE machine
      IMPLICIT NONE
      REAL (KIND=dp) :: solvene
      REAL (KIND=dp) :: n_e
    END FUNCTION solvene

    double precision function zeroin(ax,bx,f,tol)
      double precision ax,bx,f,tol
      external f
    end function zeroin

  END INTERFACE

END MODULE INTERFACES
