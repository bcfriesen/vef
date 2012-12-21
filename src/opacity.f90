! bound-free opacity for H I. Taken from Auer & Mihalas, ApJ, 156, 157
! (1969)
function sigma_ik(i, nu)
      use precision_mod
      use interfaces
      implicit none
      real(kind=dp) :: sigma_ik
      integer :: i
      real(kind=dp) :: nu
      character(len=*), parameter :: whoami = 'sigma_ik'

      select case (i)
        case (1)
          sigma_ik = 2.815d29 / nu**3
        case (2)
          sigma_ik = 8.8d27 / nu**3
        case default
          call stop_exit(1, whoami, 'invalid bound-free opacity index')
        end select
end function sigma_ik


! bound-bound opacity for H I. Taken from Auer & Mihalas, ApJ, 156, 157
! (1969)
function sigma_ij(i, j, nu)
      use precision_mod
      use interfaces
      implicit none
      real(kind=dp) :: sigma_ij
      integer :: i, j
      real(kind=dp) :: nu
      character(len=*), parameter :: whoami = 'sigma_ij'
      ! Lyman-alpha frequency
      real(kind=dp), parameter :: nu_12 = 2.466d15
      ! natural line width
      real(kind=dp), parameter :: delta_nu = 1.0d11

      select case (i)
        case (1)
          select case (j)
            case (2)
              sigma_ij = 0.01105 * exp(-(nu - nu_12)**2 / delta_nu**2) / delta_nu
            case default
              call stop_exit(1, whoami, 'invalid bound-bound opacity index')
          end select
        case default
          call stop_exit(1, whoami, 'invalid bound-bound opacity index')
        end select
end function sigma_ij

! free-free opacity for H I. Taken from Auer & Mihalas, ApJ, 156, 157
! (1969)
function sigma_kk(nu)
      use precision_mod
      use interfaces
      use const
      use const, only: h_planck, k_boltz
      use global, only: temp
      implicit none
      real(kind=dp) :: sigma_kk
      real(kind=dp) :: nu, nu_min
      real(kind=dp), parameter :: nu_paschen_limit = 3.656d14
      character(len=*), parameter :: whoami = 'sigma_kk'

      nu_min = min(nu, nu_paschen_limit)

      sigma_kk = 3.69d8 / (nu**3 * sqrt(temp)) * exp((h_planck * nu_min) / (k_boltz * temp))
end function sigma_kk
