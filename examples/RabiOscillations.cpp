#include <Lindblad.hpp>

static const int numIters = 10 * 1000;

int main() {
  std::vector<Amplitude> rhoInitial(4, 0);
  // Start in ground state
  rhoInitial[0] = 1.0;
  MasterEqn meqn(2, &rhoInitial[0]);
  double g = 1.0;
  meqn.addCoupling(Coupling(0, 1, g));
  for (int i = 0; i < numIters; ++i) {
    meqn.takeStep();
  }
  return 0;
}
