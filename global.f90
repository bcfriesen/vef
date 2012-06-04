module global

  ! # of optical depth points
  integer, parameter :: n_depth_pts = 100
  ! # of direction cosine points
  integer, parameter :: n_mu_pts = 51
  ! # of wavelength points
  integer, parameter :: n_wl_pts = 1

  ! maximum optical depth to consider
  real, parameter :: tau_max = 1.0e+4
  ! minimum non-zero optical depth to consider
  real, parameter :: tau_min = 1.0e-8

  ! Thermalization parameter for source function in isotropic, monochromatic
  ! scattering (Milne-Eddington problem). eps = 1 means pure LTE; eps = 0 means
  ! pure scattering (like SYNOW).
  real, parameter :: me_therm_parm = 1.0e-4
  ! optical depth grid
  real, dimension( n_depth_pts ) :: tau_grid
  ! direction cosine grid
  real, dimension( n_mu_pts ) :: mu_grid
  ! wavelength grid
  real, dimension( n_wl_pts ) :: wl_grid
  ! Eddington factor f_K = K / J
  real, dimension( n_depth_pts, n_wl_pts ) :: vef_f_k
  ! Eddington factor f_H = \int_0^1 j(\mu) \mu d\mu / J
  real, dimension( 2, n_wl_pts ) :: vef_f_h
  ! VEF values from previous iteration. Need these so we can find out when
  ! they've converged.
  real, dimension( n_depth_pts, n_wl_pts ) :: vef_f_k_old
  real, dimension( 2, n_wl_pts ) :: vef_f_h_old
  ! source function
  real, dimension( n_depth_pts, n_wl_pts ) :: source_fn
  ! mean intensity
  real, dimension( n_depth_pts, n_wl_pts ) :: j_lambda
  ! specific intensity
  real, dimension( n_depth_pts, n_mu_pts, n_wl_pts ) :: i_lambda

  ! blackbody temperature at depth
  ! TODO: add temperature dependence later
  real, parameter :: temp = 1.0e4

end module global
