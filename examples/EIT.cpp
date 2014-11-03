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

static const int numIters = 1000;

int main() {
  int dim = 4;
  double Delta = -0.1;
  double OmegaB = 0.01;
  double OmegaR = 1.0;

  MasterEqn meqn(dim);
  meqn.addCoupling(3, 0, -OmegaR / 2.0 / sqrt(6.0));
  meqn.addCoupling(3, 2, -OmegaR / 2.0 / sqrt(6.0));
  meqn.addCoupling(0, 0, OmegaB / 8.0);
  meqn.addCoupling(2, 2, -OmegaB / 8.0);
  meqn.addCoupling(3, 3, -Delta);

  std::vector<Amplitude> rhoInitial(dim * dim, 0);
  rhoInitial[1 + 1 * dim] = 1.0;
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]);
  for (int i = 0; i < numIters; ++i) {
    evolution.takeStep();
  }
  return 0;
}
