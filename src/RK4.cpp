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

void RK4::advance(double* t, void* ctx) {
  evaluateRHS(&y[0], &k1[0], *t, ctx);
  for (size_t i = 0; i < y.size(); ++i) {
    work[i] = y[i] + 0.5 * dt * k1[i];
  }
  evaluateRHS(&work[0], &k2[0], *t + 0.5 * dt, ctx);
  for (size_t i = 0; i < y.size(); ++i) {
    work[i] = y[i] + 0.5 * dt * k2[i];
  }
  evaluateRHS(&work[0], &k3[0], *t + 0.5 * dt, ctx);
  for (size_t i = 0; i < y.size(); ++i) {
    work[i] = y[i] + dt * k3[i];
  }
  evaluateRHS(&work[0], &k4[0], *t + dt, ctx);
  double prefactor = dt / 6.0;
  for (size_t i = 0; i < y.size(); ++i) {
    y[i] += prefactor * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }
  *t += dt;
}
