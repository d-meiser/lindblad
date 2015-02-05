[![Build Status](https://travis-ci.org/d-meiser/lindblad.png?branch=master)](https://travis-ci.org/d-meiser/lindblad)
[![Build status](https://ci.appveyor.com/api/projects/status/nsun22swq5f2tn15/branch/master?svg=true)](https://ci.appveyor.com/project/d-meiser/lindblad/branch/master)
[![Coverage Status](https://coveralls.io/repos/d-meiser/lindblad/badge.svg?branch=master)](https://coveralls.io/r/d-meiser/lindblad?branch=master)

lindblad
========

Quantum optical master equation simulations.

The lindblad library provides functionality for numerical simulations
with open quantum systems.  The focus of the library is on medium sized
problems with on the order of tens of quantum levels.


## Example:

```
#include <Lindblad.hpp>
#include <vector>

static const int numIters = 1000;

int main() {
  MasterEqn meqn(2);
  double g = 1.0;
  meqn.addCoupling(0, 1, g);
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]);
  for (int i = 0; i < numIters; ++i) {
    evolution.takeStep();
  }
  return 0;
}
```


## Obtaining lindblad

The lindblad library is distributed in source form.  It can be
downloaded from the 
[lindblad github repository](https://github.com/d-meiser/lindblad).
Configuration is via cmake.  On linux the library is downloaded,
configured, built, and installed with the following commands:

```
git clone https://github.com/d-meiser/lindblad
cd lindblad
mkdir build
cd build
cmake ..
make -j4
make install
```

## Building PETSc

Several of the examples require the
[PETSc](http://www.mcs.anl.gov/petsc/) library.  To build PETSc, blas
and lapack headers and libraries are needed.  A minimal configuration of
PETSc can be created with the script
[./utilities/get_petsc.sh](./utilities/get_petsc.sh).  On ubuntu, the
following works for me:
```
sudo apt-get install libblas-dev libatlas-dev liblapack-dev
cd $LINDBLAD_ROOT_DIRECTORY
sh ./utilities/get_petsc.sh
```
`$LINDBLAD_ROOT_DIRECTORY` is the directory in which this README file
lives.  The `get_petsc.sh` script builds a minimal PETSc library and
installs it into the `${LINDBLAD_ROOT_DIRECTORY}/petsc-cpucmplx`
directory.

Once PETSc is built and installed lindblad can be configured with PETSc
support as follows:

```
cd ${LINDBLAD_ROOT_DIRECTORY}/build
cmake -DPETSC_DIR=../petsc-cpucmplx ..
make
```


## Overview of the source code

- [include](include): Public API of the lindblad library
- [examples](examples): Some examples.
- [include/detail](include/detail): Internal headers.
- [tests](tests): Tests for lindblad.
- [src](src): Implementation files of the library.


## Questions, comments, and issues

Questions, comments, and issues should be reported via the [github issue
tracker](https://github.com/d-meiser/lindblad/issues).

