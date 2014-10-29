#include <gtest/gtest.h>
#include <RK4.hpp>

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

