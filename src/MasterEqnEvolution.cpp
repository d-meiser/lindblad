#include <MasterEqnEvolution.hpp>
#include <RK4.hpp>
#include <vector>

static void applyRhs(double* x, double* y, double t, void* ctx);

struct MasterEqnEvolutionContext {
  MasterEqnEvolutionContext() : impl(0) {}
  MasterEqnEvolutionContext(MasterEqnEvolution::Impl* p) : impl(p) {}
  MasterEqnEvolution::Impl *impl;
};

struct MasterEqnEvolution::Impl {
  Impl(const MasterEqn& eqn, const Amplitude* initialState)
      : meqn(eqn) {
    int d = eqn.getDim();
    integrator =
        new RK4(2 * d * d, 0, (const double*)initialState, &applyRhs, 1.0e-2);
    state.resize(d * d);
    std::copy(initialState, initialState + d * d, state.begin());
    ctx = MasterEqnEvolutionContext(this);
  }
  MasterEqn meqn;
  MasterEqnEvolutionContext ctx;
  std::vector<Amplitude> state;
  Integrator* integrator;
};


MasterEqnEvolution::MasterEqnEvolution(const MasterEqn& eqn,
                                       const Amplitude* initialState)
    : impl(new MasterEqnEvolution::Impl(eqn, initialState)) {
}

MasterEqnEvolution::~MasterEqnEvolution() {
  delete impl;
}

double MasterEqnEvolution::getTime() const {
  return impl->integrator->getTime();
}

void MasterEqnEvolution::takeStep() {
  impl->integrator->takeStep(&impl->ctx);
}

const Amplitude* MasterEqnEvolution::getState() const {
  return (const Amplitude*) impl->integrator->getState();
}

static void applyRhs(double* x, double* y, double t, void* ctx) {
  MasterEqnEvolutionContext* meCtx = (MasterEqnEvolutionContext*)ctx;
  meCtx->impl->meqn.apply(meCtx->impl->meqn.getDim(), (const Amplitude*)x,
                          (Amplitude*)y);
}

