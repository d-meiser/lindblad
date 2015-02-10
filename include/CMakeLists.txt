set(LINDBLAD_HEADERS
  Amplitude.hpp
  Lindblad.hpp
  MasterEqnEvolution.hpp
  MasterEqn.hpp
  SparseMatrix.hpp
  )
install(FILES ${LINDBLAD_HEADERS} ${PROJECT_BINARY_DIR}/config.h ${PROJECT_BINARY_DIR}/src/LindbladExport.h
        DESTINATION include)

