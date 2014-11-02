#include <gtest/gtest.h>
#include <MasterEqnEvolution.hpp>
#include <MasterEqn.hpp>

TEST(MasterEqnEvolution, Constructor) { 
  int dim = 4;
  MasterEqn eqn(dim); 
  std::vector<Amplitude> psi(dim * dim);
  MasterEqnEvolution(eqn, &psi[0]);
}

TEST(MasterEqnEvolution, GetTime) {
  MasterEqn meqn(2);
  meqn.addCoupling(1, 0, 1.0);

  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]);
  EXPECT_FLOAT_EQ(0, evolution.getTime());
}

TEST(MasterEqnEvolution, TakeStep) {
  MasterEqn meqn(2);
  meqn.addCoupling(1, 0, 1.0);
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]); 
  evolution.takeStep();
  EXPECT_LE(0, evolution.getTime());
}

TEST(MasterEqnEvolution, SpontaneousDecay) {
  MasterEqn meqn(2);
  double gamma = 1.0;
  meqn.addDecay(0, 1, gamma);
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[3] = 1.0;
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]);

  while (evolution.getTime() < 2.0) {
    evolution.takeStep();
  }
  const Amplitude* finalState = evolution.getState();
  double trace = 0;
  for (int i = 0; i < 2; ++i) {
    trace += finalState[i + i * 2].real();
  }
  double EPS = 1.0e-8;
  EXPECT_GE(trace, 1.0 - EPS);
  EXPECT_LE(trace, 1.0 + EPS);
  double pe = finalState[1 + 1 * 2].real();
  double finalTime = evolution.getTime();
  EXPECT_GE(pe, 1.0 * exp(-gamma * finalTime) - EPS);
  EXPECT_LE(pe, 1.0 * exp(-gamma * finalTime) + EPS);
}

TEST(MasterEqn, RabiOscillations) {
  MasterEqn meqn(2);
  double g = 1.0;
  meqn.addCoupling(0, 1, g);

  std::vector<Amplitude> rhoInitial(4, 0);
  // Start in ground state
  rhoInitial[0] = 1.0;
  MasterEqnEvolution evolution(meqn, &rhoInitial[0]);

  while (evolution.getTime() < 2.0) {
    evolution.takeStep();
  }
  const Amplitude* finalState = evolution.getState();
  double trace = 0;
  for (int i = 0; i < 2; ++i) {
    trace += finalState[i + i * 2].real();
  }
  double EPS = 1.0e-8;
  EXPECT_GE(trace, 1.0 - EPS);
  EXPECT_LE(trace, 1.0 + EPS);
  double pe = finalState[1 + 1 * 2].real();
  double finalTime = evolution.getTime();
  EXPECT_GE(pe, pow(sin(0.5 * g * finalTime), 2) - EPS);
  EXPECT_LE(pe, pow(sin(0.5 * g * finalTime), 2) + EPS);
}

