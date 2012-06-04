!
!
!
      subroutine locate(xarray,n,x,j)
!     -------------------------------
      integer j,n
      real*8 x,xarray(n)
!***************************************************************
! search in an ordered table.
! adapted from numerical recipes
! 29/apr/95: changed first .gt. to .ge. to adapt for search 
!            in the line lists (equal wavelengths can happen,
!            the orginal version would get stuck if this occurs)
!***************************************************************
      integer jl,jm,ju
!
      jl=0
      ju=n+1
!-- search loop jump point
   10 continue
      if(ju-jl .gt. 1) then
       jm=(ju+jl)/2
!-- changed .gt. to .ge. (phh, 29/apr/95)
       if( (xarray(n) .ge. xarray(1)) &
           .eqv. (x .ge. xarray(jm))) then
        jl=jm
       else
        ju=jm
       endif
       goto 10
      endif
!--
      j=jl
      return
      END
