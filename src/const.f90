module const
  use precision_mod
  implicit none
  real(kind=dp), parameter :: a2cm = 1.0d-8 !< convert Angstrom to cm
  REAL (KIND=dp), PARAMETER :: pi = 3.14159D0 !< \f$\pi\f$
  REAL (KIND=dp), PARAMETER :: m_e = 9.109D-28 !< electron rest mass
  REAL (KIND=dp), PARAMETER :: k_boltz = 1.38065D-16 !< Boltzmann's constant
  REAL (KIND=dp), PARAMETER :: h_planck = 6.626D-27 !< Planck's constant
  REAL (KIND=dp), PARAMETER :: h_bar = h_planck / (2.0 * pi) !< reduced Planck's constant
  REAL (KIND=dp), PARAMETER :: hbar = h_planck/(2.0D0*pi) !< reduced Planck's constant
  REAL (KIND=dp), PARAMETER :: c_light = 3.0D10 !< speed of light
  REAL (KIND=dp), PARAMETER :: ev2erg = 1.602D-12 !< converts eV to erg

end module const
