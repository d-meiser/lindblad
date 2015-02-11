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
#include <Integrator.hpp>
using namespace Lindblad;
using namespace Lindblad::Detail;

#include <algorithm>

class MockEuler : public Integrator {
 public:
  MockEuler(int dim, double time, double* state,
            void (*f)(double* x, double* y, double t, void* ctx))
      : Integrator(dim, time, state, f), y(dim), dt(1.0e-6) {
      y.resize(dim);
      std::copy(state, state + dim, y.begin());
      k1.resize(dim);
      }
 private:
  virtual const double* getCurrentState() const {
    return &y[0];
  }
  virtual void advance(double* t, double* dt, void* ctx) {
    evaluateRHS(&y[0], &k1[0], *t, ctx);
    for (size_t i = 0; i < y.size(); ++i) {
      y[i] += *dt * k1[i];
    }
    *t += *dt;
  }
  std::vector<double> y;
  std::vector<double> k1;
  double dt;
  virtual MockEuler* makeCopy() const {
    return new MockEuler(*this);
  }
};

struct FooCtx {
  size_t dim;
};
void foo(double* x, double* y, double t, void* ctx) {
  struct FooCtx* fooCtx = (struct FooCtx*)ctx;
  for (size_t i = 0; i < fooCtx->dim; ++i) {
    y[i] = x[i];
  }
}

TEST(Integrator, Constructor) {
  std::vector<double> y(10);
  MockEuler integ(10, 0, &y[0], &foo);
}

TEST(Integrator, TakeStep) {
  std::vector<double> x(10, 2.0);
  MockEuler integ(10, 0, &x[0], &foo);
  struct FooCtx ctx = {10};
  integ.takeStep(&ctx);
  EXPECT_TRUE(integ.getTime() > 0);
}

TEST(Integrator, GetState) {
  std::vector<double> x(10, 2.0);
  MockEuler integ(10, 0, &x[0], &foo);
  struct FooCtx ctx = {10};
  integ.takeStep(&ctx);
  std::vector<double> y(10);
  std::copy(integ.getState(), integ.getState() + 10, y.begin());
  EXPECT_TRUE(y[0] > x[0]);
}

TEST(Integrator, EvaluateRHS) {
  std::vector<double> x(10, 2.0);
  MockEuler integ(10, 0, &x[0], &foo);
  struct FooCtx ctx = {10};
  std::vector<double> y(10);
  integ.evaluateRHS(&x[0], &y[0], 0, &ctx);
  for (int i = 0; i < 10; ++i) {
    EXPECT_DOUBLE_EQ(y[i], x[i]);
  }
}

TEST(Integrator, GetTimeStep) {
  MockEuler integ(0, 0, 0, &foo);
  integ.setTimeStep(1.0e-5);
  double dt = integ.getTimeStep();
  EXPECT_DOUBLE_EQ(1.0e-5, dt);
}
