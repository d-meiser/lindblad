#include <qsys.hpp>
#include <cmath>

static const int numIters = 1000;

int main() {
  int dim = 4;
  double Delta = -0.1;
  double OmegaB = 0.01;
  double OmegaR = 1.0;

  std::vector<Amplitude> rhoInitial(dim * dim, 0);
  rhoInitial[1 + 1 * dim] = 1.0;
  MasterEqn meqn(dim, &rhoInitial[0]);
  meqn.addCoupling(Coupling(3, 0, -OmegaR / 2.0 / sqrt(6.0)));
  meqn.addCoupling(Coupling(3, 2, -OmegaR / 2.0 / sqrt(6.0)));
  meqn.addCoupling(Coupling(0, 0, OmegaB / 8.0));
  meqn.addCoupling(Coupling(2, 2, -OmegaB / 8.0));
  meqn.addCoupling(Coupling(3, 3, -Delta));
  for (int i = 0; i < numIters; ++i) {
    meqn.takeStep();
  }
  return 0;
}
