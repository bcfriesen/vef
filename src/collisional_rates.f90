function c_ij(i, j)
      use precision_mod
      use global
      implicit none
      real(kind=dp) :: c_ij
      integer :: i, j
end function c_ij


function c_ji(j, i)
      use precision_mod
      use global
      implicit none
      real(kind=dp) :: c_ji
      integer :: j, i
end function c_ji


function c_ik(i)
      use precision_mod
      use global
      implicit none
      real(kind=dp) :: c_ik
      integer :: i
end function c_ik


function c_ki(i)
      use precision_mod
      use global
      implicit none
      real(kind=dp) :: c_ki
      integer :: i
end function c_ki
