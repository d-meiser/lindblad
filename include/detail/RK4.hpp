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
#ifndef RK4_HPP
#define RK4_HPP

#include <Integrator.hpp>
#include <vector>

class RK4 : public Integrator {
 public:
  RK4(int dim, double time, const double* state,
      void (*f)(double* x, double* y, double t, void* ctx));

 private:
  std::vector<double> k1, k2, k3, k4, y, work;
  virtual const double* getCurrentState() const;
  virtual void advance(double* t, double* dt, void* ctx);
  virtual RK4* makeCopy() const;
};

#endif

