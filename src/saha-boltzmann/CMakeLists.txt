SET (sb_solver_SOURCES
    saha-boltzmann/atomicdata.f90
    saha-boltzmann/f_ij.f90
    saha-boltzmann/fill_atomic_data.f90
    saha-boltzmann/part_func.f90
    saha-boltzmann/saha.f90
    saha-boltzmann/solvene.f90
    saha-boltzmann/zeroin.f
    )

include_directories(
    ${PROJECT_BINARY_DIR}/src
    )
