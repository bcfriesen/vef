!> This is the function the root finder uses to calculate \f$n_e\f$.
!! @param n_e free electron number density
FUNCTION solvene(n_e)
  USE atomicdata
  USE globalvars
  USE machine
  IMPLICIT NONE
  REAL (KIND=dp) :: solvene
  REAL (KIND=dp) :: n_e

  REAL (KIND=dp) :: ng !< total number of free particles (free electrons + atoms)
  REAL (KIND=dp) :: denom !< denominator of 2nd term of equation we need to solve to get \f$n_e\f$
  REAL (KIND=dp) :: term1 !< summation of \f$Y_i\f$ in denominator
  REAL (KIND=dp) :: term2 !< summation of \f$f_ij\f$ in denominator
  REAL (KIND=dp), EXTERNAL :: saha, f_ij
  INTEGER :: i !< loops
  INTEGER :: j !< loops

  ! get particle density from ideal gas law
  ng = p_g / (k_b * t)
  denom = 0.0D0
  term1 = 0.0D0
  term2 = 0.0D0
  ! here I would sum over all elements 1 through 92, but the partition
  ! function
  ! crashes when i>1 since I only have data for hydrogen. I'd fix it but I'm
  ! running out of time on this assignment.
  DO i = 1, 1
    ! sum over all the ionization states of each element
    DO j = 0, i
      term2 = term2 + dble(j)*f_ij(i, j, n_e, t)
    END DO
    term1 = term1 + y(i)*term2
  END DO
  denom = term1
  solvene = ng - n_e*(1.0D0+(1.0D0/denom))
  RETURN
END FUNCTION solvene
