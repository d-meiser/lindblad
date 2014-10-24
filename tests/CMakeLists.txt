message(STATUS "GTEST_INCLUDE_DIRS=${GTEST_INCLUDE_DIRS}")
set(TESTS
    test_leftApply
    test_rightApply
    test_Coupling
   )

foreach(t ${TESTS})
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} qsys gtest_main gtest)
  add_test(${t} ${t})
endforeach()

