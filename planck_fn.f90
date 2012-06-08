! Planck function
function planck_fn( lambda, temp )
  implicit none

  real :: planck_fn
  real :: lambda, temp
  real, parameter :: h_planck = 6.626e-27, &
                     c_light  = 2.998e+10, &
                     k_boltz  = 1.380e-16

  planck_fn = ( 2.0 * h_planck * c_light**2 / lambda**5 ) / &
              ( exp( ( h_planck * c_light / lambda ) / k_boltz * temp ) - 1.0 )
end function planck_fn
