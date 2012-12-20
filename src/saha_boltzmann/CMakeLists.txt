SET (sb_solver_SOURCES
    atomicdata.f90
    globalvars.f90
    f_ij.f90
    fill_atomic_data.f90
    interfaces.f90
    machine.f90
    main.f90
    part_func.f90
    saha.f90
    solvene.f90
    zeroin.f
    )

ADD_LIBRARY (saha_boltzmann ${sb_solver_SOURCES})
