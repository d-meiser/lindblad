/*
Copyright 2014 Dominic Meiser

This file is part of lindblad.

lindblad is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

lindblad is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with lindblad.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <Lindblad.hpp>
#include <cmath>
#include <vector>
#include <cstdio>

static const int numIters = 1000;

struct SystemParameters {
  double OmegaR, OmegaB, Delta, gamma, Gamma, deltaB;
};

static int readParameter(int n, const char** argv, double* param,
                  const char* name);
static int getParameters(int argn, const char** argv, SystemParameters* parameters);
static int printDensityMatrix(const Amplitude* rho, int N);

int main(int argn, const char** argv) {
  int dim = 4;
  SystemParameters params;
  getParameters(argn, argv, &params);

  MasterEqn meqn(dim);
  meqn.addCoupling(0, 0, params.OmegaB / 2.0);
  meqn.addCoupling(0, 1, params.deltaB / 2.0 / sqrt(2.0));
  meqn.addCoupling(1, 2, params.deltaB / 2.0 / sqrt(2.0));
  meqn.addCoupling(2, 2, -params.OmegaB / 2.0);
  meqn.addCoupling(0, 3, 1.0 / 4.0 * params.OmegaR / sqrt(6.0));
  meqn.addCoupling(2, 3, -1.0 / 4.0 * params.OmegaR / sqrt(6.0));
  meqn.addCoupling(3, 3, -params.Delta);
  meqn.addDecay(0, 0, params.gamma);
  meqn.addDecay(1, 1, params.gamma);
  meqn.addDecay(2, 2, params.gamma);
  meqn.addDecay(3, 3, params.gamma);
  meqn.addDecay(0, 3, params.Gamma / 3.0);
  meqn.addDecay(1, 3, params.Gamma / 3.0);
  meqn.addDecay(2, 3, params.Gamma / 3.0);

  std::vector<Amplitude> rhoInitial(dim * dim, 0);
  rhoInitial[1 + 1 * dim] = 1.0;
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]);
  evolution.setTimeStep(1.0e-9);
  for (int i = 0; i < numIters; ++i) {
    evolution.takeStep();
    printDensityMatrix(evolution.getState(), dim);
  }
  return 0;
}

int getParameters(int argn, const char** argv, SystemParameters* parameters) {
  int ierr;

  if (argn == 1) {
    parameters->OmegaR = 1.25e6 * 2.0 * M_PI;
    parameters->OmegaB = 700.0e3 * 2.0 * M_PI * 0.01;
    parameters->Delta = 0.0;
    parameters->gamma = 1.0e3 * 2.0 * M_PI;
    parameters->Gamma = 6.0e6 * 2.0 * M_PI;
    parameters->deltaB = 0.0;
  } else if (argn < 7) {
    printf("Insufficient parameters provided.");
    return -1;
  } else {
    ierr = readParameter(1, argv, &parameters->OmegaR, "OmegaR");
    ierr = readParameter(2, argv, &parameters->OmegaB, "OmegaB");
    ierr = readParameter(3, argv, &parameters->Delta, "Delta");
    ierr = readParameter(4, argv, &parameters->gamma, "gamma");
    ierr = readParameter(5, argv, &parameters->Gamma, "Gamma");
    ierr = readParameter(6, argv, &parameters->deltaB, "deltaB");
  }
  return ierr;
}

int readParameter(int n, const char** argv, double* param,
                  const char* name) {
  if (sscanf(argv[n], "%lf", param) != 1) {
    printf("Failed to read argument %d: %s\n", n, name);
    return -1;
  }
  return 0;
}

int printDensityMatrix(const Amplitude* rho, int N) {
  int i, j;
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      printf("%+1.5lf %+1.5lf  ", rho[i * N + j].real(), rho[i * N + j].imag());
    }
    printf("\n");
  }
}
