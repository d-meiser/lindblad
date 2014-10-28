#include <gtest/gtest.h>
#include <qsys.hpp>

TEST(MasterEquation, Constructor) {
  MasterEqnRhs rhs;
  rhs.addCoupling(Coupling(1, 0, 1.0));
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEquation* meqn = new MasterEquation(2, &rhoInitial[0], &rhs);
  EXPECT_TRUE(meqn != 0);
  delete meqn;
}

TEST(MasterEquation, GetTime) {
  MasterEqnRhs rhs;
  rhs.addCoupling(Coupling(1, 0, 1.0));
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEquation meqn(2, &rhoInitial[0], &rhs);
  EXPECT_FLOAT_EQ(0, meqn.getTime());
}

TEST(MasterEquation, TakeStep) {
  MasterEqnRhs rhs;
  rhs.addCoupling(Coupling(1, 0, 1.0));
  std::vector<Amplitude> rhoInitial(4, 0);
  rhoInitial[0] = 1.0;
  MasterEquation meqn(2, &rhoInitial[0], &rhs);
  meqn.takeStep();
  EXPECT_LE(0, meqn.getTime());
}
