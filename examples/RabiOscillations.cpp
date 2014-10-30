#include <qsys.hpp>

static const int numIters = 10 * 1000;

int main() {
  MasterEqnRhs rhs;
  double g = 1.0;
  rhs.addCoupling(Coupling(0, 1, g));
  std::vector<Amplitude> rhoInitial(4, 0);
  // Start in ground state
  rhoInitial[0] = 1.0;
  MasterEquation meqn(2, &rhoInitial[0], &rhs);
  for (int i = 0; i < numIters; ++i) {
    meqn.takeStep();
  }
  return 0;
}
