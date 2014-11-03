set(EXAMPLES
    RabiOscillations
    EIT)
foreach(e ${EXAMPLES})
  add_executable(${e} ${e}.cpp)
  target_link_libraries(${e} lindblad)
  add_test(${e} ${e})
endforeach()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plotRabiOscillationsData.py
    ${CMAKE_CURRENT_BINARY_DIR}/plotRabiOscillationsData.py)

