#include <gtest/gtest.h>
#include <qsys.hpp>
#include <vector>

TEST(Coupling, Constructor) {
  Coupling a(1, 0, 1.0);
}

template <typename T>
void MY_EXPECT_EQ(T a, T b, size_t i) {
  EXPECT_EQ(a, b) << " i == " << i;
}

template <>
void MY_EXPECT_EQ(Amplitude a, Amplitude b, size_t i) {
  EXPECT_FLOAT_EQ(a.real(), b.real()) << "real, i == " << i;
  EXPECT_FLOAT_EQ(a.imag(), b.imag()) << "imag, i == " << i;
}

template <typename IterA, typename IterB>
void EXPECT_RANGES_EQ(IterA beginA, IterA endA, IterB beginB, IterB endB) {
  size_t i = 0;
  while (beginA != endA && beginB != endB) {
    MY_EXPECT_EQ(*beginA, *beginB, i);
    ++beginA;
    ++beginB;
    ++i;
  }
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

