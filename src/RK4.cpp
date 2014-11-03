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
#include <RK4.hpp>
#include <iostream>

void RK4::buildIntegratorData(size_t dim, const double* state, double t) {
  k1.resize(dim);
  k2.resize(dim);
  k3.resize(dim);
  k4.resize(dim);
  y.resize(dim);
  work.resize(dim);
  std::copy(state, state + dim, y.begin());
}

const double* RK4::getCurrentState() const { return &y[0]; }

void RK4::advance(double* t, double* dt, void* ctx) {
  evaluateRHS(&y[0], &k1[0], *t, ctx);
  for (size_t i = 0; i < y.size(); ++i) {
    work[i] = y[i] + 0.5 * *dt * k1[i];
  }
  evaluateRHS(&work[0], &k2[0], *t + 0.5 * *dt, ctx);
  for (size_t i = 0; i < y.size(); ++i) {
    work[i] = y[i] + 0.5 * *dt * k2[i];
  }
  evaluateRHS(&work[0], &k3[0], *t + 0.5 * *dt, ctx);
  for (size_t i = 0; i < y.size(); ++i) {
    work[i] = y[i] + *dt * k3[i];
  }
  evaluateRHS(&work[0], &k4[0], *t + *dt, ctx);
  double prefactor = *dt / 6.0;
  for (size_t i = 0; i < y.size(); ++i) {
    y[i] += prefactor * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }
  *t += *dt;
}

RK4* RK4::makeCopy() const {
  return new RK4(*this);
}
