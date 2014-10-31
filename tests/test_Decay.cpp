#include <gtest/gtest.h>
#include <Lindblad.hpp>
#include <testUtils.hpp>

TEST(Decay, Constructor) {
  Decay d(0, 1, 0.2);
}

TEST(Decay, EToG) {
  Decay d(0, 1, 0.2);
  std::vector<Amplitude> rho1(4, 1);
  std::vector<Amplitude> rho2(4, 0);
  d.apply(2, &rho1[0], &rho2[0]);
}

TEST(Decay, EToGTrace) {
  Decay d(0, 1, 0.2);
  CheckLindbladTraceProperty(d, 2);
}

TEST(Decay, EToGHermiticity) {
  Decay d(0, 1, 0.2);
  CheckLindbladHermiticityProperty(d, 2);
}
