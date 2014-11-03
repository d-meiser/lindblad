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
#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <cstddef>

class Integrator {
 public:
  Integrator(int dim, double time, const double* state,
             void (*f)(double* x, double* y, double t, void* ctx))
      : d(dim),
        t(time),
        dt(1.0e-3),
        rhs(f) {}
  void takeStep(void* ctx);
  double getTime() const { return t; }
  const double* getState() const;
  size_t getDim() const { return d; }
  void evaluateRHS(double* in, double* out, double t, void* ctx);
  Integrator* copy() const;
  void setTimeStep(double deltat) { dt = deltat; }
  double getTimeStep() const { return dt; }

 private:
  int d;
  double t;
  double dt;
  void (*rhs)(double* x, double* y, double t, void* ctx);
  virtual const double* getCurrentState() const = 0;
  virtual void advance(double* t, double* dt, void* ctx) = 0;
  virtual Integrator* makeCopy() const = 0;
};

#endif
