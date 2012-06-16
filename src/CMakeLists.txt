set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(vef_SOURCES
calc_flux.f90
calc_moments.f90
calc_rmsd.f90
calc_source_fn.f90
calc_vefs.f90
global.f90
interfaces.f90
main.f90
planck_fn.f90
precision_mod.f90
solve_rte.f90
solve_scatt_prob.f90
stop_exit.f90
write_moments.f90
write_source_fn.f90
write_spectrum.f90
write_vefs.f90
)

set (EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(vef ${vef_SOURCES})

#------------------------------BEGIN ADD LIBRARIES------------------------------
find_package(BLAS)
if (BLAS_FOUND)
  include_directories(${BLAS_INCLUDE_DIRS})
  target_link_libraries(vef ${BLAS_LIBRARIES})
endif (BLAS_FOUND)
find_package(LAPACK)
if (LAPACK_FOUND)
  include_directories(${LAPACK_INCLUDE_DIRS})
  target_link_libraries(vef ${LAPACK_LIBRARIES})
endif (LAPACK_FOUND)
#-------------------------------END ADD LIBRARIES-------------------------------
