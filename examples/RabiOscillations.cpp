#include <MasterEqn.hpp>
#include <MasterEqnEvolution.hpp>
#include <Amplitude.hpp>
#include <vector>

static const int numIters = 10 * 1000;

int main() {
  MasterEqn meqn(2);
  double g = 1.0;
  meqn.addCoupling(0, 1, g);

  std::vector<Amplitude> rhoInitial(4, 0);
  // Start in ground state
  rhoInitial[0] = 1.0;
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]);
  for (int i = 0; i < numIters; ++i) {
    evolution.takeStep();
  }
  return 0;
}
