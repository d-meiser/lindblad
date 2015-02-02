include_directories(${PROJECT_SOURCE_DIR}/include/detail)
set(TESTS
    test_leftApply
    test_rightApply
    test_Coupling
    test_Decay
    test_MasterEqn
    test_MasterEqnEvolution
    test_Memory
    test_Integrator
    test_RK4
    test_SourceSink
   )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
foreach(t ${TESTS})
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} lindblad gtest_main gtest)
  add_test(${t} ${t})
endforeach()

