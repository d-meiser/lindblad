set(EXAMPLES
    RabiOscillations
    EIT)
foreach(e ${EXAMPLES})
  add_executable(${e} ${e}.cpp)
  target_link_libraries(${e} lindblad)
  add_test(${e} ${e})
endforeach()

set(PETSC_EXAMPLES
    SteadyStateEIT
    SteadyStateEIA
    SteadyState)
find_package(DL)
find_path(PETSC_INCLUDE_DIR petsc.h ${PETSC_DIR}/include)
message(STATUS "PETSC_INCLUDE_DIR==${PETSC_INCLUDE_DIR}")
find_library(PETSC_LIBS petsc ${PETSC_DIR}/lib)
message(STATUS "PETSC_LIBS==${PETSC_LIBS}")
if(PETSC_INCLUDE_DIR AND PETSC_LIBS)
  set (PETSC_FOUND TRUE)
  include_directories(${PETSC_INCLUDE_DIR}/mpiuni ${PETSC_INCLUDE_DIR})
endif()
if(PETSC_FOUND)
foreach(e ${PETSC_EXAMPLES})
  add_executable(${e} ${e}.cpp)
  target_link_libraries(${e} lindblad ${PETSC_LIBS} lapack blas pthread ${DL_LIBRARIES})
  add_test(${e} ${e} -pc_type none)
endforeach()
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotRabiOscillationsData.py
    ${CMAKE_CURRENT_BINARY_DIR}/plotRabiOscillationsData.py)

