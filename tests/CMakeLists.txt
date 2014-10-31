include_directories(${PROJECT_SOURCE_DIR}/include/detail)
set(TESTS
    test_leftApply
    test_rightApply
    test_Coupling
    test_Decay
    test_MasterEqnRhs
    test_MasterEqn
    test_Integrator
    test_RK4
   )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
foreach(t ${TESTS})
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} lindblad gtest_main gtest)
  add_test(${t} ${t})
endforeach()

