      subroutine stop_exit(stat, who, message)
      implicit none

      integer :: stat
      character(len=*) :: who, message

      write(*,*) 'ERROR'
      write(*,*) 'status ID: ', stat
      write(*,*) 'who called error handler: ', who
      write(*,*) 'message: ', message

      stop
      return
      end subroutine stop_exit
