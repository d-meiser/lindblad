cmake_minimum_required(VERSION 2.6.2)
project(LINDBLAD)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)
enable_testing()
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")
include(CheckCXXCompilerFlag)
include(GenerateExportHeader)
set(LINDBLAD_VER_MAJOR 0)
set(LINDBLAD_VER_MINOR 0)
set(LINDBLAD_VER_PATCH 1)
set(LINDBLAD_VERSION
    "${LINDBLAD_VER_MAJOR}.${LINDBLAD_VER_MINOR}.${LINDBLAD_VER_PATCH}")
message(STATUS "LINDBLAD_VERSION == ${LINDBLAD_VERSION}")

option(LINDBLAD_BUILD_TESTS "Whether to build tests" OFF)
option(LINDBLAD_BUILD_DOCS "Whether to build documentation" OFF)
option(LINDBLAD_ENABLE_COVERAGE "Instrument tests for coverage" OFF)
option(LINDBLAD_WITH_ADDRESS_SANITIZER "Whether to build with address sanitizer enabled" OFF)

if(LINDBLAD_WITH_ADDRESS_SANITIZER)
  set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
  set(CMAKE_REQUIRED_FLAGS "-fsanitize=address")
  CHECK_CXX_COMPILER_FLAG(-fsanitize=address HAVE_ADDRESS_SANITIZER)
  if(HAVE_ADDRESS_SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
  endif()
  set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
endif()
if(LINDBLAD_BUILD_DOCS)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doxyfile @ONLY)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc")
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc"
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
  else()
    message(STATUS "Doxygen not found -- not building docs")
  endif()
endif()


set(CMAKE_REQUIRED_FLAGS "-fopenmp-simd")
CHECK_CXX_COMPILER_FLAG(-fopenmp-simd HAVE_PRAGMA_OMP_SIMD)
if(HAVE_PRAGMA_OMP_SIMD)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -fopenmp-simd")
endif()
set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})

if(CMAKE_COMPILER_IS_GNUCC OR "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  if(LINDBLAD_ENABLE_COVERAGE)
    add_definitions("-ftest-coverage -fprofile-arcs")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftest-coverage -fprofile-arcs")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftest-coverage -fprofile-arcs")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ftest-coverage -fprofile-arcs")
    set(COV_LIBRARIES gcov)
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
    ${PROJECT_BINARY_DIR}/src
    )
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(include)
if(LINDBLAD_BUILD_TESTS)
  add_subdirectory(gtest-1.7.0)
  include_directories(${gtest_SOURCE_DIR}/include)
  add_subdirectory(tests)
endif()

configure_file(LindbladConfig.h.in LindbladConfig.h)

