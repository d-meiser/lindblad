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

Integrator* Integrator::copy() const {
  makeCopy();
}
