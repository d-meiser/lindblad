cmake_minimum_required(VERSION 2.6.2)
project(LINDBLAD)
enable_testing()
include(GenerateExportHeader)

option(LINDBLAD_BUILD_TESTS "Whether to build tests" OFF)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(src)
if(LINDBLAD_BUILD_TESTS)
  add_subdirectory(gtest-1.7.0)
  include_directories(${gtest_SOURCE_DIR}/include)
  add_subdirectory(tests)
endif()

