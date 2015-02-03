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
#include <RK45.hpp>

#include <algorithm>
#include <iostream>
#include <cmath>

const double RK45::cs[6] = {0.0,         1.0 / 4.0, 3.0 / 8.0,
                            12.0 / 13.0, 1.0,       1.0 / 2.0};
const double RK45::as[25] = {
  1.0 / 4.0,       0,                0,                0,               0,
  3.0 / 32.0,      9.0 / 32.0,       0,                0,               0,
  1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0,  0,               0,
  439.0 / 216.0,   -8.0,             3680.0 / 513.0,   -845.0 / 4104.0, 0,
  -8.0 / 27.0,     2.0,              -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0};
const double RK45::b5[6] = {16.0 / 135.0, 0, 6656.0 / 12825.0,
                            28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0};
const double RK45::b4[6] = {25.0 / 216.0, 0, 1408.0 / 2565.0, 2197.0 / 4104.0,
                            -1.0 / 5.0, 0};

RK45::RK45(int dim, double time, const double *state,
           void (*f)(double *x, double *y, double t, void *ctx))
    : Integrator(dim, time, state, f), y(dim), ks(6), work(dim), sol4(dim),
      sol5(dim) {
  std::copy(state, state + dim, &y[0]);
  for (int i = 0; i < 6; ++i) {
    ks[i].resize(dim);
  }
  for (int i = 0; i < dim; ++i) {
    y[i] = state[i];
  }
}

RK45::~RK45() {}
const double* RK45::getCurrentState() const {
  return &y[0];
}

void RK45::advance(double* t, double* dt, void* ctx) {
  double error;
  do {
    evaluateRHS(&y[0], &ks[0][0], *t, ctx);
  for (int row = 1; row < 6; ++row) {
    computeWork(row, *dt);
    evaluateRHS(&work[0], &ks[row][0], *t + cs[row], ctx);
  }

  for (std::size_t i = 0; i < y.size(); ++i) {
    sol5[i] = y[i];
    for (int j = 0; j < 6; ++j) {
      sol5[i] += b5[j] * *dt * ks[j][i];
    }
  }
  for (std::size_t i = 0; i < y.size(); ++i) {
    sol4[i] = y[i];
    for (int j = 0; j < 6; ++j) {
      sol4[i] += b4[j] * *dt * ks[j][i];
    }
  }

  error = errorEstimate(sol4, sol5);
  if (error > 1.0e-9) {
    *dt /= 2.0;
  }
  } while (error > 1.0e-9);

  for (std::size_t i = 0; i < y.size(); ++i) {
    y[i] = sol5[i];
  }

  *t += *dt;
  if (error < 1.0e-11) {
    *dt *= 2.0;
  }
}

void RK45::computeWork(int row, double dt) {
  for (std::size_t i = 0; i < y.size(); ++i) {
    work[i] = y[i];
    for (int col = 0; col < row; ++col) {
      work[i] += as[(row - 1) * 5 + col] * dt * ks[col][i];
    }
  }
}

double RK45::errorEstimate(const std::vector<double>& x, const std::vector<double>& y) const {
  double l2Error = 0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    l2Error += (x[i] - y[i]) * (x[i] - y[i]);
  }
  return sqrt(l2Error);
}

RK45* RK45::makeCopy() const {
  return new RK45(*this);
}
