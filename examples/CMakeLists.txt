set(EXAMPLES
    RabiOscillations
    EIT)
foreach(e ${EXAMPLES})
  add_executable(${e} ${e}.cpp)
  target_link_libraries(${e} lindblad)
  add_test(${e} ${e})
endforeach()

