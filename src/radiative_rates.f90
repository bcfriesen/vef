function r_ij(i, j)
      use precision_mod
      use global
      implicit none
      real(kind=dp) :: r_ij
      integer :: i, j
end function r_ij

function r_ik(i)
      use precision_mod
      use global
      implicit none
      real(kind=dp) :: r_ik
      integer :: i
end function r_ik

function r_ki(i)
      use precision_mod
      use global
      implicit none
      real(kind=dp) :: r_ki
      integer :: i
end function r_ki

function z_ji(j, i)
      use precision_mod
      use global
      implicit none
      real(kind=dp) :: z_ji
      integer :: j, i
end function z_ji


