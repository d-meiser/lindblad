#include <gtest/gtest.h>
#include <qsys.hpp>
#include <vector>
#include <testUtils.hpp>

TEST(Coupling, Constructor) {
  Coupling a(1, 0, 1.0);
}

TEST(Coupling, ApplySpSm) {
  Coupling a(0, 0, 1.0);
  std::vector<Amplitude> rho1(4, 1);
  std::vector<Amplitude> rho2(4, 0);
  a.apply(2, &rho1[0], &rho2[0]);
  Amplitude expected[] = {0, Amplitude(0, -1), Amplitude(0, 1), 0};
  EXPECT_RANGES_EQ(expected, expected + 4, rho2.begin(), rho2.end());
}

TEST(Coupling, ApplySmSp) {
  Coupling a(1, 1, 1.0);
  std::vector<Amplitude> rho1(4, 1);
  std::vector<Amplitude> rho2(4, 0);
  a.apply(2, &rho1[0], &rho2[0]);
  Amplitude expected[] = {0, Amplitude(0, 1), Amplitude(0, -1), 0};
  EXPECT_RANGES_EQ(expected, expected + 4, rho2.begin(), rho2.end());
}

TEST(Coupling, ApplySm) {
  Coupling a(0, 1, 1.0);
  std::vector<Amplitude> rho1(4, 1);
  std::vector<Amplitude> rho2(4, 0);
  a.apply(2, &rho1[0], &rho2[0]);
  Amplitude expected[] = {0, 0, 0, 0};
  EXPECT_RANGES_EQ(expected, expected + 4, rho2.begin(), rho2.end());
}

TEST(Coupling, ApplySp) {
  Coupling a(1, 0, 1.0);
  std::vector<Amplitude> rho1(4, 1);
  std::vector<Amplitude> rho2(4, 0);
  a.apply(2, &rho1[0], &rho2[0]);
  Amplitude expected[] = {0, 0, 0, 0};
  EXPECT_RANGES_EQ(expected, expected + 4, rho2.begin(), rho2.end());
}

TEST(Coupling, TraceSp) {
  Coupling a(1, 0, 1.0);
  CheckLindbladTraceProperty(a, 2);
}

TEST(Coupling, HermiticitySp) {
  Coupling a(1, 0, 1.0);
  CheckLindbladHermiticityProperty(a, 2);
}

TEST(Coupling, TraceSpSm) {
  Coupling a(1, 1, 1.0);
  CheckLindbladTraceProperty(a, 2);
}

TEST(Coupling, HermiticitySpSm) {
  Coupling a(1, 1, 1.0);
  CheckLindbladHermiticityProperty(a, 2);
}

TEST(Coupling, HermiticityComplexMatElem) {
  Coupling a(0, 1, Amplitude(2.3, 4.8));
  CheckLindbladHermiticityProperty(a, 2);
}