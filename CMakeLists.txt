cmake_minimum_required(VERSION 2.6.2)
project(LINDBLAD)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)
enable_testing()
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")
include(CheckCXXCompilerFlag)

option(LINDBLAD_BUILD_TESTS "Whether to build tests" OFF)
option(LINDBLAD_ENABLE_COVERAGE "Instrument tests for coverage" OFF)

set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -fsanitize=address")
CHECK_CXX_COMPILER_FLAG("-fsanitize=address" HAVE_ADDRESS_SANITIZER)
if(HAVE_ADDRESS_SANITIZER)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")
endif()

if(CMAKE_COMPILER_IS_GNUCC OR "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  if(LINDBLAD_ENABLE_COVERAGE)
    add_definitions("-ftest-coverage -fprofile-arcs")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftest-coverage -fprofile-arcs")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftest-coverage -fprofile-arcs")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ftest-coverage -fprofile-arcs")
    set(COV_LIBRARIES gcov)
  endif()
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
    )
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(include)
if(LINDBLAD_BUILD_TESTS)
  add_subdirectory(gtest-1.7.0)
  include_directories(${gtest_SOURCE_DIR}/include)
  add_subdirectory(tests)
endif()

try_compile(HAVE_PRAGMA_OMP_SIMD ${PROJECT_BINARY_DIR}/cmake
            ${PROJECT_SOURCE_DIR}/cmake/try_compile/pragma_omp_simd.cpp)
configure_file(config.h.in config.h)

include(GenerateExportHeader)
generate_export_header(lindblad
    EXPORT_FILE_NAME LindbladExport.h
    EXPORT_MACRO_NAME LINDBLAD_API
    )
