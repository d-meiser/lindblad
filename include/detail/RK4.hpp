#ifndef RK4_HPP
#define RK4_HPP

#include <Integrator.hpp>
#include <vector>

class RK4 : public Integrator {
 public:
  RK4(int dim, double time, const double* state,
      void (*f)(double* x, double* y, double t, void* ctx), double deltat)
      : Integrator(dim, time, state, f), dt(deltat) {}

 private:
  double dt;
  std::vector<double> k1, k2, k3, k4, y, work;
  virtual void buildIntegratorData(size_t dim, const double* state, double t);
  virtual const double* getCurrentState() const;
  virtual void advance(double* t, void* ctx);
  virtual RK4* makeCopy() const;
};

#endif

