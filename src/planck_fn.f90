! Planck function
function planck_fn( lambda, temp )
  use precision_mod
  use const, only: h_planck, c_light, k_boltz
  implicit none

  real(kind=dp) :: planck_fn
  real(kind=dp) :: lambda, temp

  planck_fn = ( 2.0d+0 * h_planck * c_light**2 / lambda**5 ) / &
              ( exp( ( h_planck * c_light / lambda ) / ( k_boltz * temp ) ) - &
              1.0d+0 )

end function planck_fn
