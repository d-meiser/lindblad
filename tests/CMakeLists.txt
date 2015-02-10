include_directories(${PROJECT_SOURCE_DIR}/include/detail)
set(TESTS
    test_Coupling
    test_Decay
    test_GeneralDecayOperator
    test_Integrator
    test_MasterEqn
    test_MasterEqnEvolution
    test_Memory
    test_RK4
    test_RK45
    test_SourceSink
    test_leftApply
    test_rightApply
    test_SparseMatrix
   )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
foreach(t ${TESTS})
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} lindblad_static gtest_main gtest)
  set_target_properties(${t} PROPERTIES COMPILE_FLAGS -DLINDBLAD_STATIC_DEFINE)
  add_test(${t} ${t})
endforeach()

