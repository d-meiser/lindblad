#include <gtest/gtest.h>
#include <SparseApply.hpp>

#include <vector>

TEST(LeftApply, SigmaPlus) {
  std::vector<Amplitude> rho1(4, 1);
  std::vector<Amplitude> rho2(4, 0);
  rightApply(1, 0, 1.0, 2, &rho1[0], &rho2[0]);
  EXPECT_FLOAT_EQ(1, rho2[0].real());
  EXPECT_FLOAT_EQ(0, rho2[0].imag());
  EXPECT_FLOAT_EQ(0, rho2[1].real());
  EXPECT_FLOAT_EQ(0, rho2[1].imag());
  EXPECT_FLOAT_EQ(1, rho2[2].real());
  EXPECT_FLOAT_EQ(0, rho2[2].imag());
  EXPECT_FLOAT_EQ(0, rho2[3].real());
  EXPECT_FLOAT_EQ(0, rho2[3].imag());
}
