module global

  use precision_mod
  ! # of optical depth points
  integer, parameter :: n_depth_pts = 100
  ! # of direction cosine points.
  integer, parameter :: n_mu_pts = n_depth_pts
  ! # of wavelength points
  integer, parameter :: n_wl_pts = 100

  ! maximum optical depth to consider
  real(kind=dp), parameter :: tau_max = 1.0d+6
  ! minimum non-zero optical depth to consider
  real(kind=dp), parameter :: tau_min = 1.0d-8

  ! minimum wavelength
  real(kind=dp), parameter :: wlmin = 1000.0d+0
  ! maximum wavelength
  real(kind=dp), parameter :: wlmax = 10000.0d+0

  ! Thermalization parameter for source function in isotropic, monochromatic
  ! scattering (Milne-Eddington problem). eps = 1 means pure LTE; eps = 0 means
  ! pure scattering (like SYNOW).
  real(kind=dp), parameter :: me_therm_parm = 1.0d-2
  ! optical depth grid
  real(kind=dp), dimension( n_depth_pts ) :: tau_grid
  ! finite-difference points on optical depth grid
  real(kind=dp), dimension( 1 : n_depth_pts-1 ) :: dtau
  real(kind=dp), dimension( 2 : n_depth_pts-1 ) :: dtau_tilde
  ! direction cosine grid
  real(kind=dp), dimension( n_mu_pts ) :: mu_grid
  ! wavelength grid (Angstrom)
  real(kind=dp), dimension( n_wl_pts ) :: wl_grid
  ! Eddington factor f_K = K / J
  real(kind=dp), dimension( n_depth_pts, n_wl_pts ) :: vef_f_k
  ! Eddington factor f_H = \int_0^1 j(\mu) \mu d\mu / J. There's one of these at
  ! the top boundary and another at the bottom.
  real(kind=dp), dimension( 2, n_wl_pts ) :: vef_f_h
  ! VEF values from previous iteration. Need these so we can find out when
  ! they've converged.
  real(kind=dp), dimension( n_depth_pts, n_wl_pts ) :: vef_f_k_old
  real(kind=dp), dimension( 2, n_wl_pts ) :: vef_f_h_old
  ! source function
  real(kind=dp), dimension( n_depth_pts, n_wl_pts ) :: source_fn
  ! first 3 moments of I
  real(kind=dp), dimension( n_depth_pts, n_wl_pts ) :: j_lambda, h_lambda, k_lambda
  ! Feautrier variables. These are symmetric (j) and anti-symmetric (h) averages
  ! of specific intensities, so they only span from 0 < mu < 1, rather than from
  ! -1 < mu < 1
  real(kind=dp), dimension( n_depth_pts, n_mu_pts, n_wl_pts ) :: little_j, little_h
  ! blackbody temperature at depth
  ! TODO: add temperature dependence
  real(kind=dp), parameter :: temp = 1.0d+4
  ! free electron density
  real(kind=dp), dimension( n_depth_pts ) :: n_e
  ! # of levels in model atom for H I
  integer, parameter :: n_levels = 16
  ! LTE level populations
  real(kind=dp), dimension( n_depth_pts, n_levels ) :: n_i_star
  REAL (KIND=dp) :: p_g !< total gas pressure
  REAL (KIND=dp) :: y(92) !< number fraction of each element

end module global
