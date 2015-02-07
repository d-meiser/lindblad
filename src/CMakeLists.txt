include_directories(${PROJECT_SOURCE_DIR}/include/detail)
include_directories(${PROJECT_BINARY_DIR})
add_library(lindblad
    OBJECT
    Coupling.cpp
    Decay.cpp
    GeneralDecayOperator.cpp
    Integrator.cpp
    MasterEqn.cpp
    MasterEqnEvolution.cpp
    RK4.cpp
    SourceSink.cpp
    RK45.cpp
    SparseApply.cpp
    )
add_library(lindblad_static $<TARGET_OBJECTS:lindblad>)
add_library(lindblad_shared SHARED $<TARGET_OBJECTS:lindblad>)
set_target_properties(lindblad_shared PROPERTIES OUTPUT_NAME lindblad)

