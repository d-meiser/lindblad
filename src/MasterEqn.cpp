#include <MasterEqn.hpp>
#include <RK4.hpp>
#include <MasterEqnRhs.hpp>

struct MasterEqnRhsContext {
  int dim;
  const MasterEqn::Impl* rhs;
};

static void applyRhs(double* x, double* y, double t, void* ctx);

struct MasterEqn::Impl {
  Impl(int d, const Amplitude* A) {
    ctx.dim = d;
    ctx.rhs = this;
    integrator = new RK4(2 * d * d, 0, (const double *)A, &applyRhs, 1.0e-2);
  }
  Impl(const Impl& other)
      : couplings(other.couplings),
        decays(other.decays),
        ctx(other.ctx),
        integrator(other.integrator->copy()) {
    ctx.rhs = this;
  }
  ~Impl() { delete integrator; }
  void apply(int dim, const Amplitude* A, Amplitude* B) const {
    std::fill(B, B + dim * dim, 0);
    for (std::vector<Coupling>::const_iterator c = couplings.begin();
         c != couplings.end(); ++c) {
      c->apply(dim, A, B);
    }
    for (std::vector<Decay>::const_iterator d = decays.begin();
         d != decays.end(); ++d) {
      d->apply(dim, A, B);
    }
  }

  std::vector<Coupling> couplings;
  std::vector<Decay> decays;
  MasterEqnRhsContext ctx;
  Integrator* integrator;
};

MasterEqn::MasterEqn(int d, const Amplitude* A) { impl = new Impl(d, A); }

MasterEqn::MasterEqn(const MasterEqn& other)
    : impl(new MasterEqn::Impl(*other.impl)) {}

MasterEqn& MasterEqn::operator=(const MasterEqn& rhs) {
  if (this != &rhs) {
    delete impl;
    impl = new MasterEqn::Impl(*rhs.impl);
  }
  return *this;
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
  impl->couplings.push_back(Coupling(m, n, a));
}

void MasterEqn::addDecay(int into, int outOf, double gamma) {
  impl->decays.push_back(Decay(into, outOf, gamma));
}

void MasterEqn::apply(int dim, const Amplitude* A, Amplitude *B) const {
  impl->apply(dim, A, B);
}


static void applyRhs(double* x, double* y, double t, void* ctx) {
  MasterEqnRhsContext* meCtx = (MasterEqnRhsContext*)ctx;
  meCtx->rhs->apply(meCtx->dim, (const Amplitude*)x, (Amplitude*)y);
}

