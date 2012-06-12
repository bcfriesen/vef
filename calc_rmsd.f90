! calculates root mean square deviation between two vectors
Function calc_rmsd(array1, array2)
  Use precision_mod
  Use interfaces, Only: stop_exit
  Implicit None
  Real :: calc_rmsd
  Real (Kind=dp), Dimension (:) :: array1, array2
  Integer :: i1
  Character (Len=*), Parameter :: whoami = 'calc_rmsd'

  If (size(array1)/=size(array2)) Call stop_exit(1, whoami, 'cannot &
    &calculate RMSD because arrays have unequal size!')

  calc_rmsd = 0.0D+0
  Do i1 = 1, size(array1)
    calc_rmsd = calc_rmsd + (array1(i1)-array2(i1))**2
  End Do
  calc_rmsd = calc_rmsd/dble(size(array1))
  calc_rmsd = sqrt(calc_rmsd)

End Function calc_rmsd
