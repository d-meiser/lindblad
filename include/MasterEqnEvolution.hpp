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
#ifndef MASTER_EQN_EVOLUTION_HPP
#define MASTER_EQN_EVOLUTION_HPP

#include <MasterEqn.hpp>
#include <LindbladExport.h>

namespace Lindblad {

/**
@brief Numerical time evolution of a master equation
*/
class LINDBLAD_API MasterEqnEvolution {
 public:
/**
@brief Construct from a master equation and an initial state
*/
  MasterEqnEvolution(const MasterEqn& eqn, const Amplitude* initialState);
  ~MasterEqnEvolution();

/**
@brief Get the current time of the master equation evolution.
*/
  double getTime() const;

/**
@brief Take one explicit time step

Note that the density matrix might be advanced by more or less than
getTimeStep if an adaptive time step solver is used.  Use getTime to
query the actual time after the step has been completed.

@sa getTime, getTimeStep, setTimeStep
*/
  void takeStep();

/**
@brief Returns the current state of the quantum system

@return A view on the density matrix.  The density matrix is in row
major order and of size \f$dim\times dim\f$ where \f$dim\f$ is the
dimension of the underlying quantum system.  Note that the returned
density matrix does not correspond to a copy of the internal state.  It
is illegal to free this pointer or to manipulate the data in any way.
If the density matrix is to be stored or manipulated it is necessary to
make a copy.
*/
  const Amplitude* getState() const;

/**
@brief Set hint for the time step size

Note that the MasterEqnEvolution is not forced to use this exact time
step size.  When adaptive time step integration routines are used the
actual time step size is usually different.

@sa getTimeStep, getTime
*/
  void setTimeStep(double dt);

/**
@brief Returns the current guess for the next time step size.
*/
  double getTimeStep() const;

  friend struct MasterEqnEvolutionContext;

 private:
  struct Impl;
  Impl* impl;
};

}

#endif

