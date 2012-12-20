!> Calculates all the necessary atomic data.
subroutine fill_atomic_data
  use atomicdata
  implicit none

  integer :: i !< loops

  ! for this homework, H I = 10 levels, H II = 1 level
  nlvl(1, 0) = 10
  nlvl(1, 1) = 1

  ! set ionization potentials (from ground)
  chiion(1, 0) = 13.6d0
  chiion(1, 1) = 0.0d0

  ! for H I level k has energy 13.6 eV - (13.6 eV / k^2)
  do i = 1, nlvl(1, 0)
    chilvl(1, 0, i) = chiion(1, 0) * (1.0d0 - (1.0d0 / dble(i)**2))
  end do
  ! for H II we only include ground
  chilvl(1, 1, 1) = 0.0D0

  ! for H I the degeneracy of level k is 2k^2
  do i = 1, nlvl(1,0)
    g_i(1, 0, i) = 2 * i**2
  end do
  ! H II is nondegenerate
  g_i(1, 1, 1) = 1

end subroutine fill_atomic_data
