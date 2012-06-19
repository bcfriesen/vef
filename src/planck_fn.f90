! Planck function
function planck_fn( lambda, temp )
  use precision_mod
  implicit none

  real(kind=dp) :: planck_fn
  real(kind=dp) :: lambda, temp
  real(kind=dp), parameter :: h_planck = 6.626d-27, &
                              c_light  = 2.998d+10, &
                              k_boltz  = 1.380d-16

  planck_fn = ( 2.0d+0 * h_planck * c_light**2 / lambda**5 ) / &
              ( exp( ( h_planck * c_light / lambda ) / ( k_boltz * temp ) ) - &
              1.0d+0 )

end function planck_fn
