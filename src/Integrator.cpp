#include <Integrator.hpp>

void Integrator::takeStep(void* ctx) {
  if (!engineInitialized) {
    buildIntegratorData(d, state, t);
    engineInitialized = true;
  }
  advance(&t, ctx);
}

const double* Integrator::getState() const {
  return getCurrentState();
}
void Integrator::evaluateRHS(double* in, double* out, double t, void* ctx) {
  rhs(in, out, t, ctx);
}

Integrator* Integrator::copy() const {
  makeCopy();
}
