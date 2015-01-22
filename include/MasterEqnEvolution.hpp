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

class MasterEqnEvolution {
 public:
  MasterEqnEvolution(const MasterEqn& eqn, const Amplitude* initialState);
  ~MasterEqnEvolution();
  double getTime() const;
  void takeStep();
  const Amplitude* getState() const;
  void setTimeStep(double dt);
  double getTimeStep() const;

  friend struct MasterEqnEvolutionContext;

 private:
  struct Impl;
  Impl* impl;
};

#endif

