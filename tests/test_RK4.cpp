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
  RK4 rk4(10, 0, 0, &foo);
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
  RK4 rk4(1, 0, &x, &exponentialDecay);
  rk4.takeStep(&ctx);
  static const double EPS = 1.0e-10;
  double t = rk4.getTime();
  EXPECT_GE(*rk4.getState(), x * exp(-gamma * t) - EPS);
  EXPECT_LE(*rk4.getState(), x * exp(-gamma * t) + EPS);
}

TEST(RK4, ExpMultipleSteps) {
  static const double gamma = 6.3;
  struct DecayCtx ctx = {gamma};
  double x = 1.7;
  RK4 rk4(1, 0, &x, &exponentialDecay);
  static const int N = 10;
  for (int i = 0; i < N; ++i) {
    rk4.takeStep(&ctx);
  }
  static const double EPS = 1.0e-10;
  double t = rk4.getTime();
  EXPECT_GE(*rk4.getState(), x * exp(-gamma * t) - EPS);
  EXPECT_LE(*rk4.getState(), x * exp(-gamma * t) + EPS);
}

