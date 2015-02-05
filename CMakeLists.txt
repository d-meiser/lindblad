cmake_minimum_required(VERSION 2.6.2)
project(LINDBLAD)
enable_testing()
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

option(LINDBLAD_BUILD_TESTS "Whether to build tests" OFF)

include_directories(${PROJECT_SOURCE_DIR}/include)
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
