#include <gtest/gtest.h>
#include <Integrator.hpp>

#include <algorithm>

class MockEuler : public Integrator {
 public:
  MockEuler(int dim, double time, double* state,
            void (*f)(double* x, double* y, double t, void* ctx))
      : Integrator(dim, time, state, f), y(dim), dt(1.0e-6) {}
 private:
  virtual void buildIntegratorData(size_t dim, const double* state, double t) {
    k1.resize(dim);
    y.resize(dim);
    std::copy(state, state + dim, y.begin());
  }
  virtual const double* getCurrentState() const {
    return &y[0];
  }
  virtual void advance(double* t, void* ctx) {
    evaluateRHS(&y[0], &k1[0], *t, ctx);
    for (size_t i = 0; i < y.size(); ++i) {
      y[i] += dt * k1[i];
    }
    *t += dt;
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
  MockEuler integ(10, 0, 0, &foo);
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
    EXPECT_FLOAT_EQ(y[i], x[i]);
  }
}


