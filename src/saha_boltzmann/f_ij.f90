!> Determines fraction of atom \f$i\f$ which is in ionization state \f$j\f$
!! @param i atom index
!! @param j ion index for atom \f$i\f$
!! @param n_e free electron number density
!! @param t temperature
FUNCTION f_ij(i, j, n_e, t)
  USE machine
  USE interfaces, ONLY: saha
  IMPLICIT NONE
  REAL (KIND=dp) :: f_ij
  INTEGER :: i
  INTEGER :: j
  REAL (KIND=dp) :: n_e
  REAL (KIND=dp) :: t
  INTEGER :: m !< loops
  INTEGER :: n !< loops
  REAL (KIND=dp) :: num !< numerator of \f$f_ij\f$ equation we derived in class
  REAL (KIND=dp) :: denom !< denominator of \f$f_ij\f$ equation we derived in class
  REAL (KIND=dp) :: term !< each term in summation in denominator

  ! the numerator of f_ij is a product of j ratios
  num = 1.0D0
  DO n = 0, j - 1
    num = num*saha(i, n, n_e, t)
  END DO
  ! the denominator of f_ij is a sum of i+1 terms. the 1st term is 1, and after
  ! that the nth term in the sum is a product of n-1 ratios (i.e., 2nd term is
  ! 1 ratio, 3rd term is 2 ratios, etc.)
  denom = 0.0D0
  term = 1.0D0
  DO n = 1, i + 1
    term = 1.0D0
    DO m = 1, n - 1
      IF (n==1) EXIT
      term = term*saha(i, m-1, n_e, t)
    END DO
    denom = denom + term
  END DO
  f_ij = num/denom
  RETURN
END FUNCTION f_ij
