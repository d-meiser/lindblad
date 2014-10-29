#include <Integrator.hpp>

void Integrator::takeStep(void* ctx) {
  if (!engineInitialized) {
    buildIntegratorData(d, state, t);
  }
  advance(&t, ctx);
}

const double* Integrator::getState() const {
  return getCurrentState();
}

