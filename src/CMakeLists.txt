include_directories(${PROJECT_SOURCE_DIR}/include/detail)
add_library(lindblad
    Coupling.cpp
    Decay.cpp
    Integrator.cpp
    MasterEqn.cpp
    MasterEqnEvolution.cpp
    Memory.cpp
    RK4.cpp
    SourceSink.cpp
    SparseApply.cpp
    )

