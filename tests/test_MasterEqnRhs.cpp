#include <gtest/gtest.h>
#include <MasterEqnRhs.hpp>
#include <testUtils.hpp>

TEST(MasterEqnRhs, Constructor) {
  MasterEqnRhs* meqn = new MasterEqnRhs();
  EXPECT_TRUE(meqn != 0);
  delete meqn;
}

TEST(MasterEqnRhs, AddCoupling) {
  MasterEqnRhs meqn;
  EXPECT_NO_THROW(meqn.addCoupling(Coupling(2, 3, Amplitude(1.0, 11.0))));
}

TEST(MasterEqnRhs, AddDecay) {
  MasterEqnRhs meqn;
  EXPECT_NO_THROW(meqn.addDecay(Decay(2, 3, 2.0)));
}

MasterEqnRhs buildRhs() {
  MasterEqnRhs meqn;
  meqn.addCoupling(Coupling(2, 1, Amplitude(2.0, -1.0)));
  meqn.addCoupling(Coupling(1, 1, Amplitude(1.0, -7.0)));
  meqn.addCoupling(Coupling(2, 0, Amplitude(3.0, -2.0)));
  meqn.addCoupling(Coupling(3, 4, Amplitude(4.0, -3.0)));
  meqn.addCoupling(Coupling(4, 2, Amplitude(7.0, -3.0)));
  meqn.addDecay(Decay(0, 0, 2.0));
  meqn.addDecay(Decay(2, 3, 1.0));
  meqn.addDecay(Decay(2, 1, 3.0));
  meqn.addDecay(Decay(5, 7, 4.0));
  return meqn;
}

TEST(MasterEqnRhs, TraceProperty) {
  MasterEqnRhs meqn = buildRhs();
  CheckLindbladTraceProperty(meqn, 8);
}

TEST(MasterEqnRhs, HermiticityProperty) {
  MasterEqnRhs meqn = buildRhs();
  CheckLindbladHermiticityProperty(meqn, 8);
}
