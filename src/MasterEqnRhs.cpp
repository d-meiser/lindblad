#include <MasterEqnRhs.hpp>

void MasterEqnRhs::addCoupling(Coupling c) {
  couplings.push_back(c);
}

void MasterEqnRhs::addDecay(Decay d) {
  decays.push_back(d);
}

void MasterEqnRhs::apply(int dim, const Amplitude *A, Amplitude *B) const {
  std::fill(B, B + dim * dim, 0);
  for (std::vector<Coupling>::const_iterator c = couplings.begin();
       c != couplings.end(); ++c) {
    c->apply(dim, A, B);
  }
  for (std::vector<Decay>::const_iterator d = decays.begin(); d != decays.end();
       ++d) {
    d->apply(dim, A, B);
  }
}

