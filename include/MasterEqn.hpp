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
#ifndef MASTEREQN_HPP
#define MASTEREQN_HPP

#include <Amplitude.hpp>

class MasterEqn {
 public:
  MasterEqn(int dim);
  MasterEqn(const MasterEqn& other);
  MasterEqn& operator=(const MasterEqn& rhs);
  ~MasterEqn();
  void addCoupling(int m, int n, Amplitude a);
  void addDecay(int into, int outOf, double gamma);
  void apply(const Amplitude* A, Amplitude* B) const;
  int getDim() const;

 private:
  struct Impl;
  Impl* impl;
};

#endif

