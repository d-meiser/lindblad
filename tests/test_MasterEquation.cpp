#include <gtest/gtest.h>
#include <qsys.hpp>

TEST(MasterEquation, Constructor) {
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEquation* meqn = new MasterEquation(2, &rhoInitial[0]);
  meqn->addCoupling(Coupling(1, 0, 1.0));
  EXPECT_TRUE(meqn != 0);
  delete meqn;
}

TEST(MasterEquation, GetTime) {
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEquation meqn(2, &rhoInitial[0]);
  meqn.addCoupling(Coupling(1, 0, 1.0));
  EXPECT_FLOAT_EQ(0, meqn.getTime());
}

TEST(MasterEquation, TakeStep) {
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEquation meqn(2, &rhoInitial[0]);
  meqn.addCoupling(Coupling(1, 0, 1.0));
  meqn.takeStep();
  EXPECT_LE(0, meqn.getTime());
}

TEST(MasterEquation, SpontaneousDecay) {
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[3] = 1.0;
  MasterEquation meqn(2, &rhoInitial[0]);
  double gamma = 1.0;
  meqn.addDecay(Decay(0, 1, gamma));

  while (meqn.getTime() < 2.0) {
    meqn.takeStep();
  }
  const Amplitude* finalState = meqn.getState();
  double trace = 0;
  for (int i = 0; i < 2; ++i) {
    trace += finalState[i + i * 2].real();
  }
  double EPS = 1.0e-8;
  EXPECT_GE(trace, 1.0 - EPS);
  EXPECT_LE(trace, 1.0 + EPS);
  double pe = finalState[1 + 1 * 2].real();
  double finalTime = meqn.getTime();
  EXPECT_GE(pe, 1.0 * exp(-gamma * finalTime) - EPS);
  EXPECT_LE(pe, 1.0 * exp(-gamma * finalTime) + EPS);
}

TEST(MasterEquation, RabiOscillations) {
  std::vector<Amplitude> rhoInitial(4, 0);
  // Start in ground state
  rhoInitial[0] = 1.0;
  MasterEquation meqn(2, &rhoInitial[0]);
  double g = 1.0;
  meqn.addCoupling(Coupling(0, 1, g));

  while (meqn.getTime() < 2.0) {
    meqn.takeStep();
  }
  const Amplitude* finalState = meqn.getState();
  double trace = 0;
  for (int i = 0; i < 2; ++i) {
    trace += finalState[i + i * 2].real();
  }
  double EPS = 1.0e-8;
  EXPECT_GE(trace, 1.0 - EPS);
  EXPECT_LE(trace, 1.0 + EPS);
  double pe = finalState[1 + 1 * 2].real();
  double finalTime = meqn.getTime();
  EXPECT_GE(pe, pow(sin(0.5 * g * finalTime), 2) - EPS);
  EXPECT_LE(pe, pow(sin(0.5 * g * finalTime), 2) + EPS);
}
