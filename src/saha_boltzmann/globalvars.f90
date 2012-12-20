!> Contains global variables that everybody can use so I don't have pass 83
!! arguments to every function
MODULE globalvars
  USE machine
  IMPLICIT NONE
  REAL (KIND=dp) :: t !< temperature
  REAL (KIND=dp) :: p_g !< total gas pressure
  REAL (KIND=dp) :: y(92) !< number fraction of each element
END MODULE globalvars
