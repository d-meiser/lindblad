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
#include <gtest/gtest.h>
#include <RK45.hpp>
using namespace Lindblad;
using namespace Lindblad::Detail;
#include <cmath>

struct FooCtx {
  size_t dim;
};
void foo(double* x, double* y, double t, void* ctx) {
  struct FooCtx* fooCtx = (struct FooCtx*)ctx;
  for (size_t i = 0; i < fooCtx->dim; ++i) {
    y[i] = x[i];
  }
}

TEST(RK45, Constructor) {
  std::vector<double> x(10);
  RK45 integrator(10, 0, &x[0], &foo);
}

struct DecayCtx {
  double gamma;
};
void exponentialDecay(double* x, double* y, double t, void* ctx) {
  struct DecayCtx *decCtx = (struct DecayCtx*)ctx;
  *y = -decCtx->gamma * *x;
}

TEST(RK45, Exp) {
  static const double gamma = 6.3;
  struct DecayCtx ctx = {gamma};
  double x = 1.7;
  RK45 integrator(1, 0, &x, &exponentialDecay);
  integrator.takeStep(&ctx);
  static const double EPS = 1.0e-10;
  double t = integrator.getTime();
  EXPECT_GE(*integrator.getState(), x * exp(-gamma * t) - EPS);
  EXPECT_LE(*integrator.getState(), x * exp(-gamma * t) + EPS);
}

