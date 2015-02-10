include_directories(${PROJECT_SOURCE_DIR}/include/detail)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
set(LINDBLAD_SRCS
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
set(LINDBLAD_HDRS
    ${PROJECT_SOURCE_DIR}/include/Amplitude.hpp
    ${PROJECT_SOURCE_DIR}/include/Lindblad.hpp
    ${PROJECT_SOURCE_DIR}/include/MasterEqn.hpp
    ${PROJECT_SOURCE_DIR}/include/MasterEqnEvolution.hpp
    ${PROJECT_SOURCE_DIR}/include/SparseMatrix.hpp
    ${PROJECT_SOURCE_DIR}/include/detail/Coupling.hpp
    ${PROJECT_SOURCE_DIR}/include/detail/Decay.hpp
    ${PROJECT_SOURCE_DIR}/include/detail/GeneralDecayOperator.hpp
    ${PROJECT_SOURCE_DIR}/include/detail/Integrator.hpp
    ${PROJECT_SOURCE_DIR}/include/detail/RK4.hpp
    ${PROJECT_SOURCE_DIR}/include/detail/RK45.hpp
    ${PROJECT_SOURCE_DIR}/include/detail/SourceSink.hpp
    ${PROJECT_SOURCE_DIR}/include/detail/SparseApply.hpp
    )

add_library(lindblad_static STATIC ${LINDBLAD_SRCS} ${LINDBLAD_HDRS})
add_library(lindblad_shared SHARED ${LINDBLAD_SRCS} ${LINDBLAD_HDRS})
generate_export_header(lindblad_shared
    EXPORT_FILE_NAME LindbladExport.h
    BASE_NAME lindblad
    EXPORT_MACRO_NAME LINDBLAD_API
    )
set_target_properties(lindblad_static PROPERTIES
	OUTPUT_NAME lindblad
	COMPILE_FLAGS -DLINDBLAD_STATIC_DEFINE
	)

install(TARGETS lindblad_static lindblad_shared DESTINATION lib)

