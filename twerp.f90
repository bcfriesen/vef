!
! Time-stamp: <19/12/2007 12:21 baron>
!
                                                                        
!                                                                       
!***********************************************************************
      subroutine twerp(x0,y0,x,y,npts)                                  
!***********************************************************************
      use precision_mod
      implicit none                                                     
      integer, intent(in) :: npts
      real(kind=dp), intent(in) :: x
      real(kind=dp), intent(out) :: y
      real(kind=dp), dimension(:), intent(in) :: x0,y0
      integer :: idx
      real(kind=dp) :: power,slope,xrat,xratx,yrat
      logical :: dolog                                                  
!***********************************************************************
!*****does a two point log interpolation if there is no sign change**** 
!*****does a two point lin interpolation if there is a  sign change**** 
!***********************************************************************
!***********************************************************************
      call locate(x0,npts,x,idx)                                        
      if(idx .eq. 0 .or. idx .eq. npts) then                            
       print *,x0(1:npts),npts,x,idx                                    
       call stop_exit(1,'error in interp')                              
      endif                                                             
      slope  = (y0(idx+1)-y0(idx)) / (x0(idx+1)-x0(idx)+1.d-70)         
      dolog = ( x0(idx)*x0(idx+1) .gt. 0.d0 .and. &
         y0(idx)*y0(idx+1) .gt. 0.d0 )                                  
      if(dolog) then                                                    
       yrat = y0(idx+1)/y0(idx)                                         
       xrat = x0(idx+1)/x0(idx)                                         
       xratx = x/x0(idx)                                                
         power  = dlog(yrat) / dlog(xrat)
         y      = y0(idx)*(xratx**power)                                
      else                                                              
         y      = y0(idx) + slope * (x-x0(idx))                         
      endif                                                             
!***********************************************************************
      return                                                            
      end subroutine twerp                                              
!***********************************************************************

