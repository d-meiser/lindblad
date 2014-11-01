#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <cstddef>

class Integrator {
 public:
  Integrator(int dim, double time, const double* state,
             void (*f)(double* x, double* y, double t, void* ctx))
      : d(dim), t(time), state(state), rhs(f), engineInitialized(false) {}
  void takeStep(void* ctx);
  double getTime() const { return t; }
  const double* getState() const;
  size_t getDim() const { return d; }
  void evaluateRHS(double* in, double* out, double t, void* ctx) {
    rhs(in, out, t, ctx);
  }
  Integrator* copy() const;

 private:
  int d;
  double t;
  const double *state;
  void (*rhs)(double* x, double* y, double t, void* ctx);
  bool engineInitialized;
  virtual void buildIntegratorData(size_t dim, const double* state,
                                   double t) = 0;
  virtual const double* getCurrentState() const = 0;
  virtual void advance(double* t, void* ctx) = 0;
  virtual Integrator* makeCopy() const = 0;
};

#endif
