!> Calculates ratio of number density of ion \f$j+1\f$ to ion \f$j\f$ of atom
!! \f$i\f$.
!! @param i atom index
!! @param j ion index for atom \f$i\f$
!! @param n_e free electron number density
!! @param t temperature
FUNCTION saha(i, j, n_e, t)
  USE atomicdata
  USE machine
  IMPLICIT NONE
  REAL (KIND=dp) :: saha
  INTEGER :: i
  INTEGER :: j
  REAL (KIND=dp) :: n_e
  REAL (KIND=dp) :: t
  REAL (KIND=dp), EXTERNAL :: part
  ! this is the Saha equation with n_e lumped together with all the
  ! temperature-dependent stuff, so what it returns is \f$n_{i,j+1}/n_{i,j}\f$
  saha = (1.0D0 / n_e ) * (dble(g_e) / (2.0D0 * pi * hbar)**3) * (part(i, j+1) / &
         part(i, j)) * (2.0D0 * pi * m_e * k_b * t)**(1.5D0) * dexp(-(chiion(i, j) * &
         ev2erg) / (k_b * t))
  RETURN
END FUNCTION saha
