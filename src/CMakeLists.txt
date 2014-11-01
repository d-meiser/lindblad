include_directories(${PROJECT_SOURCE_DIR}/include/detail)
add_library(lindblad
    Lindblad.cpp
    Integrator.cpp
    RK4.cpp
    MasterEqnRhs.cpp
    SparseApply.cpp
    Coupling.cpp
    Decay.cpp
    )

