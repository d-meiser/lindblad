message(STATUS "GTEST_INCLUDE_DIRS=${GTEST_INCLUDE_DIRS}")
set(TESTS
    test_leftApply
   )

foreach(t ${TESTS})
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} gtest_main gtest)
  add_test(${t} ${t})
endforeach()

