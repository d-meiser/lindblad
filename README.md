[![Build Status](https://travis-ci.org/d-meiser/lindblad.png?branch=master)](https://travis-ci.org/d-meiser/lindblad)

lindblad
========

Quantum optical master equation simulations.

The lindblad library provides functionality for numerical simulations
with open quantum systems.  The focus of the library is on medium sized
problems with on the order of tens of quantum levels.


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

## Example:

```
#include <Lindblad.hpp>

static const int numIters = 10 * 1000;

int main() {
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEqn meqn(2, &rhoInitial[0]);
  double g = 1.0;
  meqn.addCoupling(Coupling(0, 1, g));
  for (int i = 0; i < numIters; ++i) {
    meqn.takeStep();
  }
  return 0;
}
```

## Questions, comments, and issues

Questions, comments, and issues should be reported via the [github issue
tracker](https://github.com/d-meiser/lindblad/issuesa).

