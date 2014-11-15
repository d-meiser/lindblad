A few examples illustrating how to use the lindblad library.
The steady state examples [SteadyState.cpp](./SteadyState.cpp),
[SteadyStateEIT.cpp](./SteadyStateEIT.cpp), and
[SteadyStateEIA.cpp](./SteadyStateEIA.cpp) require
[PETSc](http://www.mcs.anl.gov/petsc/).

- [RabiOscillations.cpp](./RabiOscillations.cpp): Dynamics of a driven two
  level atom.
- [EIT.cpp](./EIT.cpp): Electromagnetically induced transparency in a four
  level system.
- [SteadyState.cpp](./SteadyState.cpp): Steady state solution for a
  driven damped two level atom.
- [SteadyStateEIT.cpp](./SteadyStateEIT.cpp): Steady state solution for
  a four level system exhibiting electromagnetically induced transparency.
- [SteadyStateEIA.cpp](./SteadyStateEIA.cpp): Steady state solution for
  a multi-levelsystem exhibiting electromagnetically induced absorption.
  This example is work in progress.
- [plotRabiOscillationsData.py](./plotRabiOscillationsData.py): Python
  script for plotting excited state population and coherence generated
  by [RabiOscillations.cpp](./RabiOscillations.cpp).
