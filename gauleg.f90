      subroutine gauleg(x1,x2,x,w,n)
      use precision_mod
      implicit none
      integer, intent(in):: n
      real(kind=dp), intent(in):: x1,x2
      real(kind=dp), intent(out):: x(:),w(:)
      integer :: M, I, J
      real(kind=dp) :: XM, EPS, XL, Z, P1, P2, P3, PP, Z1
      parameter (eps=3.d-14)

      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=dcos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(dabs(z-z1).gt.eps)go to 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      end subroutine gauleg
