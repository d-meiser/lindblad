cmake_minimum_required(VERSION 2.8.6)
project(QSYS)
enable_testing()
include(GenerateExportHeader)

option(QSYS_BUILD_TESTS "Whether to build tests" OFF)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(src)
if(QSYS_BUILD_TESTS)
  add_subdirectory(gtest-1.7.0)
  include_directories(${gtest_SOURCE_DIR}/include)
  add_subdirectory(tests)
endif()

