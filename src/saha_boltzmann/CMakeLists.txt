SET (sb_solver_SOURCES
    saha_boltzmann/atomicdata.f90
    saha_boltzmann/f_ij.f90
    saha_boltzmann/fill_atomic_data.f90
    saha_boltzmann/part_func.f90
    saha_boltzmann/saha.f90
    saha_boltzmann/solvene.f90
    saha_boltzmann/zeroin.f
    )

include_directories(
    ${PROJECT_BINARY_DIR}/src
    )
