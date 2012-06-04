      subroutine stop_exit(status,message)
!     ------------------------------------
      integer status
      character*(*) message
!***********************************************************************
!  this subroutine is used to print a message to stderr and 
!  stop the code with a return value, if possible.
!        version 1.0 of 20/may/94 by eab
!-- notes:
!  this routine is machine dependent (call to exit-routine only)
!  and this vresion is by far not complete ...
! exit() is for IBM RS/6000
!***********************************************************************
      write(0,'(a)')message
      stop
      return
      end
