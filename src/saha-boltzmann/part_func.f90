!> Calculates the partition function for atom \f$i\f$ in ionization state
!! \f$j\f$.
!! @param i atom index
!! @param j ion index for atom \f$i\f$
FUNCTION part(i, j)
  USE atomicdata
  USE global
  USE precision_mod
  use const, only: ev2erg, k_boltz
  IMPLICIT NONE
  REAL (KIND=dp) :: part
  INTEGER :: i
  INTEGER :: j

  INTEGER :: m !< loops

  part = 0.0D0
  DO m = 1, nlvl(i, j)
    part = part + dble(g_i(i,j,m))*dexp(-(chilvl(i,j,m)*ev2erg)/(k_boltz*temp))
  END DO
  RETURN
END FUNCTION part

