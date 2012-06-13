      subroutine stop_exit(stat, who, message)
      implicit none

      integer :: stat
      character(len=*) :: who, message

      write(*,*)
      write(*,'(a80)') '-------------------------------------&
      &ERROR--------------------------------------'
      write(*,*) 'STATUS ID: ', stat
      write(*,*) 'LOCATION: ', who
      write(*,*) 'MESSAGE: ', message
      write(*,*)

      stop
      return
      end subroutine stop_exit
