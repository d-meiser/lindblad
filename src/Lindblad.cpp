#include <Lindblad.hpp>
#include <RK4.hpp>
#include <MasterEqnRhs.hpp>

struct MasterEqnRhsContext {
  int dim;
  const MasterEqnRhs* rhs;
};

static void applyRhs(double* x, double* y, double t, void* ctx) {
  MasterEqnRhsContext* meCtx = (MasterEqnRhsContext *)ctx;
  meCtx->rhs->apply(meCtx->dim, (const Amplitude*) x, (Amplitude*) y);
}

struct MasterEqn::MasterEqnImpl {
  MasterEqnImpl(int d, const Amplitude* A) {
    ctx.dim = d;
    ctx.rhs = &rhs;
    integrator = new RK4(2 * d * d, 0, (const double *)A, &applyRhs, 1.0e-2);
  }
  ~MasterEqnImpl() { delete integrator; }
  MasterEqnRhs rhs;
  MasterEqnRhsContext ctx;
  Integrator* integrator;
};

MasterEqn::MasterEqn(int d, const Amplitude *A)
     {
  impl = new MasterEqnImpl(d, A);
}

MasterEqn::~MasterEqn() {
  delete impl;
}

double MasterEqn::getTime() const {
  return impl->integrator->getTime();
}

void MasterEqn::takeStep() {
  impl->integrator->takeStep(&impl->ctx);
}

const Amplitude* MasterEqn::getState() const {
  return (const Amplitude*) impl->integrator->getState();
}

void MasterEqn::addCoupling(int m, int n, Amplitude a) {
  impl->rhs.addCoupling(Coupling(m, n, a));
}

void MasterEqn::addDecay(int into, int outOf, double gamma) {
  impl->rhs.addDecay(Decay(into, outOf, gamma));
}
