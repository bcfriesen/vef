SET (CMAKE_INCLUDE_CURRENT_DIR ON)
SET (vef_SOURCES
calc_flux.f90
calc_moments.f90
calc_rmsd.f90
calc_source_fn.f90
calc_vefs.f90
const.f90
global.f90
grid_funcs.f90
interfaces.f90
line_profiles.f90
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

SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

ADD_EXECUTABLE(vef ${vef_SOURCES})

#------------------------------BEGIN ADD LIBRARIES------------------------------
IF (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
  SET (CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -framework Accelerate")
ELSE ()
  FIND_PACKAGE (BLAS)
  IF (BLAS_FOUND)
    INCLUDE_DIRECTORIES (${BLAS_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES (vef ${BLAS_LIBRARIES})
  ENDIF (BLAS_FOUND)

  FIND_PACKAGE (LAPACK)
  IF (LAPACK_FOUND)
    INCLUDE_DIRECTORIES (${LAPACK_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES (vef ${LAPACK_LIBRARIES})
  ENDIF (LAPACK_FOUND)
ENDIF (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
#-------------------------------END ADD LIBRARIES-------------------------------
