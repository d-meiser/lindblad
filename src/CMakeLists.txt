include_directories(${PROJECT_SOURCE_DIR}/include/detail)
add_library(lindblad
    MasterEqn.cpp
    MasterEqnEvolution.cpp
    Integrator.cpp
    RK4.cpp
    SparseApply.cpp
    Coupling.cpp
    Decay.cpp
    )

