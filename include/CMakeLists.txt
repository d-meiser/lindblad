set(LINDBLAD_HEADERS
  Amplitude.hpp
  Lindblad.hpp
  MasterEqnEvolution.hpp
  MasterEqn.hpp
  SparseMatrix.hpp
  )
install(FILES
        ${LINDBLAD_HEADERS}
        ${PROJECT_BINARY_DIR}/LindbladConfig.h
        ${PROJECT_BINARY_DIR}/LindbladExport.h
        DESTINATION include)

