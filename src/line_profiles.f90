function gaussian(lambda, norm, lambda_0, sigma)
  use precision_mod
  implicit none

  real(kind=dp) :: gaussian
  real(kind=dp) :: lambda, norm, lambda_0, sigma

  gaussian = norm * exp( -( lambda - lambda_0 )**2 / ( 2.0d+0 * sigma**2 ) )

end function gaussian
