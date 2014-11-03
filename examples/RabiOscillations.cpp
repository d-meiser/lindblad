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
#include <vector>
#include <iostream>

static void dumpState(const MasterEqnEvolution& evolution);

int main() {
  MasterEqn meqn(2);
  double g = 2.0;
  meqn.addCoupling(0, 1, g);
  double delta = 0.1;
  meqn.addCoupling(1, 1, delta);
  double gamma = 1.0;
  meqn.addDecay(0, 1, gamma);

  std::vector<Amplitude> rhoInitial(4, 0);
  // Start in ground state
  rhoInitial[0] = 1.0;
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]);
  evolution.setTimeStep(1.0e-3);
  const int numSteps = 40;
  const int itersPerStep = 200;
  dumpState(evolution);
  for (int i = 0; i < numSteps; ++i) {
    for (int j = 0; j < itersPerStep; ++j){
      evolution.takeStep();
    }
    dumpState(evolution);
  }
  return 0;
}

void dumpState(const MasterEqnEvolution& evolution) {
  std::cout << evolution.getTime() << " ";
  const Amplitude* state = evolution.getState();
  std::cout << state[0].real() << "\n";
}
