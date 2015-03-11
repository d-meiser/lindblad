/*
Copyright 2014-2015 Dominic Meiser

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

/*
\par
This example illustrates how to use the Lindblad library to simulate the
time evolution of a driven, damped two level atom.

\par
First we need to include the Lindblad header and bring the Lindblad namespace into scope:
*/
#include <Lindblad.hpp>
using namespace Lindblad;

/*
We also need a couple of standard C++ headers:
*/
#include <vector>
#include <iostream>

/*
The <tt>dumpState</tt> function prints the real part of the ground state
amplitude and the modulus of the coherence between ground and excited
state.  Its definition is at the end of this source file.
*/
static void dumpState(const MasterEqnEvolution& evolution);

int main() {
/*
The main object needed for describing master equations is \ref MasterEqn.
To construct one we need to supply the dimension of the quantum system.
For a two level system we have a Hilbert space dimension of 2.  We use
the convention that 0 corresponds to the ground state and 1 corresponds
to the excited state.
*/
  MasterEqn meqn(2);
/*
Next we need to add the different physical processes to the master
equation.  The Hermitian parts are added using \ref addCoupling.
\f$g\f$ is the coupling constant between ground and excited state.
*/
  double g = 2.0;
  meqn.addCoupling(0, 1, g);
/*
\f$ \delta \f$ is the detuning between the external laser field and the
atomic transition.
*/
  double delta = 0.1;
  meqn.addCoupling(1, 1, delta);
/*
\f$ \gamma \f$ is the decay rate of the excited state.  Dissipative
processes like decay or decoherence are added using \ref addDecay.
*/
  double gamma = 1.0;
  meqn.addDecay(0, 1, gamma);

/*
Quantum mechanical states and density matrices are described by arrays
in row major order.  That means that \f$ \rho_{i,j} = \rho[i * dim + j]
\f$ where \f$dim\f$ is the Hilbert space dimension.  Here we use a
<tt>std::vector</tt> to manage states more easily.
*/
  std::vector<Amplitude> rhoInitial(4, 0);
/*
We set the ground state amplitude to 1 (the excited state amplitude was
set to zero during initialization in the <tt> std::vector </tt>
constructor).
*/
  rhoInitial[0] = 1.0;
/*
Next we create a MasterEqnEvolution object in order to compute the
evolution of the quantum system starting from the initial state
<tt>rhoInitial</tt>.
*/
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]);
/*
We choose a time step size of 1.0e-3.  For a fixed step Runge-Kutta
integrator of 4th order (the default in Lindblad) this results in an
estimated single step error of 1.0e-12 and a total error of 1.0e-3 for a
time interval of length 1.  In this example we dump the state of the
system every <tt>itersPerStep</tt> steps with a total number of
<tt>numSteps * itersPerStep</tt> time steps.
*/
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
  std::cout << state[0].real() << " " << abs(state[1]) << "\n";
}
