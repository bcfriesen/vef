Subroutine gauleg(x1, x2, x, w)
  use precision_mod
  Implicit None
  Real(kind=dp), Intent (In) :: x1, x2
  Real(kind=dp), dimension(:), Intent (Out) :: x, w
  Integer :: m, n, i, j
  Real(kind=dp) :: xm, xl, z, p1, p2, p3, pp, z1
  real(kind=dp), Parameter :: eps=3.D-14

  n = size(x(:))
  m = (n+1)/2
  xm = 0.5D0*(x2+x1)
  xl = 0.5D0*(x2-x1)
  Do i = 1, m
    z = cos(3.141592654D0*(i-.25D0)/(n+.5D0))
100 Continue
    p1 = 1.D0
    p2 = 0.D0
    Do j = 1, n
      p3 = p2
      p2 = p1
      p1 = ((2.D0*j-1.D0)*z*p2-(j-1.D0)*p3)/j
    End Do
    pp = n*(z*p1-p2)/(z*z-1.D0)
    z1 = z
    z = z1 - p1/pp
    If (abs(z-z1)>eps) Go To 100
    x(i) = xm - xl*z
    x(n+1-i) = xm + xl*z
    w(i) = 2.D0*xl/((1.D0-z*z)*pp*pp)
    w(n+1-i) = w(i)
  End Do
  Return
End Subroutine gauleg
