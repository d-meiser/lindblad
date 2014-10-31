include_directories(${PROJECT_SOURCE_DIR}/include/detail)
add_library(qsys
    qsys.cpp
    Integrator.cpp
    RK4.cpp
    MasterEqnRhs.cpp
    SparseApply.cpp
    )

