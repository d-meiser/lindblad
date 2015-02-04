include_directories(${PROJECT_SOURCE_DIR}/include/detail)
include_directories(${PROJECT_BINARY_DIR})
add_library(lindblad
    Coupling.cpp
    Decay.cpp
    Integrator.cpp
    MasterEqn.cpp
    MasterEqnEvolution.cpp
    RK4.cpp
    SourceSink.cpp
    RK45.cpp
    SparseApply.cpp
    )

