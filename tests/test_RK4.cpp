#include <gtest/gtest.h>
#include <RK4.hpp>
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

TEST(RK4, Constructor) {
  RK4 rk4(10, 0, 0, &foo, 1.0e-5);
}

struct DecayCtx {
  double gamma;
};
void exponentialDecay(double* x, double* y, double t, void* ctx) {
  struct DecayCtx *decCtx = (struct DecayCtx*)ctx;
  *y = -decCtx->gamma * *x;
}

TEST(RK4, Exp) {
  static const double gamma = 6.3;
  struct DecayCtx ctx = {gamma};
  double x = 1.7;
  static const double dt = 1.0e-3;
  RK4 rk4(1, 0, &x, &exponentialDecay, dt);
  rk4.takeStep(&ctx);
  static const double EPS = 1.0e-10;
  EXPECT_TRUE(*rk4.getState() > x * exp(-gamma * dt) - EPS);
  EXPECT_TRUE(*rk4.getState() < x * exp(-gamma * dt) + EPS);
}

TEST(RK4, ExpMultipleSteps) {
  static const double gamma = 6.3;
  struct DecayCtx ctx = {gamma};
  double x = 1.7;
  static const double dt = 1.0e-3;
  RK4 rk4(1, 0, &x, &exponentialDecay, dt);
  static const int N = 10;
  for (int i = 0; i < N; ++i) {
    rk4.takeStep(&ctx);
  }
  static const double EPS = 1.0e-10;
  EXPECT_TRUE(*rk4.getState() > x * exp(-gamma * N * dt) - EPS);
  EXPECT_TRUE(*rk4.getState() < x * exp(-gamma * N * dt) + EPS);
}

