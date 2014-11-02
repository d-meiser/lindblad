#include <MasterEqn.hpp>
#include <Coupling.hpp>
#include <Decay.hpp>
#include <vector>

struct MasterEqn::Impl {
  Impl(int d) : dim(d) {}
  Impl(const Impl& other)
      : dim(other.dim), couplings(other.couplings), decays(other.decays) {}
  Impl& operator=(const Impl& other) {
    if (this != &other) {
      dim = other.dim;
      couplings = other.couplings;
      decays = other.decays;
    }
    return *this;
  }
  ~Impl() {}
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
  int getDim() const { return dim; }

  int dim;
  std::vector<Coupling> couplings;
  std::vector<Decay> decays;
};

MasterEqn::MasterEqn(int d) { impl = new Impl(d); }

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

void MasterEqn::addCoupling(int m, int n, Amplitude a) {
  impl->couplings.push_back(Coupling(m, n, a));
}

void MasterEqn::addDecay(int into, int outOf, double gamma) {
  impl->decays.push_back(Decay(into, outOf, gamma));
}

void MasterEqn::apply(int dim, const Amplitude* A, Amplitude *B) const {
  impl->apply(dim, A, B);
}

int MasterEqn::getDim() const {
  return impl->getDim();
}

