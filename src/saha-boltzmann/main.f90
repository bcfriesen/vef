!> Calculates level and ion populations in thermal equilibrium using the
!! Saha-Boltzmann equation.
PROGRAM hw7
  USE atomicdata
  USE globalvars
  USE machine
  USE interfaces, ONLY: zeroin, solvene, f_ij
  IMPLICIT NONE

  REAL (KIND=dp) :: n_e !< free electron number density
  NAMELIST /params/ p_g, y

  ! set gas pressure to something reasonable for the sun (this will be
  ! overridden by the namelist)
  p_g = 1.0d5

  ! set composition to pure hydrogen by default (this will be overridden by the
  ! namelist)
  y(:) = 0.0d0
  y(1) = 1.0d0

  ! read pressure and abundances from stdin
  READ (*, params)

  ! set atomic data
  call fill_atomic_data

  OPEN (UNIT = 11, FILE = 'output.dat', STATUS='new')
  ! sample logarithmically a wide range in temperature. surely everything
  ! interesting that happens to hydrogen will happen within this temperature
  ! range...
  t = 10.0d0
  DO
    t = t * 1.01
    ! solve non-linear equation for n_e
    n_e = zeroin(1.0D-10, 1.0D20, solvene, 1.0D-5)
    ! write out number fractions of n_HI and n_HII
    WRITE (11, *) t, f_ij(1, 0, n_e, t), f_ij(1, 1, n_e, t)
    IF (t > 1.0d5) EXIT
  END DO
  CLOSE (11)
  STOP
END PROGRAM hw7

