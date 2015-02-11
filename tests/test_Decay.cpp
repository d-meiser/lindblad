/*
Copyright 2014-2015 Dominic Meiser

This file is part of lindblad.

lindblad is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

lindblad is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with lindblad.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <gtest/gtest.h>
#include <Decay.hpp>
using namespace Lindblad;
using namespace Lindblad::Detail;
#include <testUtils.hpp>

TEST(Decay, Constructor) {
  Decay d(0, 1, 0.2);
}

TEST(Decay, EToG) {
  double gamma = 0.2;
  Decay d(0, 1, gamma);
  std::vector<Amplitude> rho1(4, 1);
  std::vector<Amplitude> rho2(4, 0);
  d.apply(2, &rho1[0], &rho2[0]);
  EXPECT_DOUBLE_EQ(gamma, rho2[0].real());
  EXPECT_DOUBLE_EQ(0.0, rho2[0].imag());
  EXPECT_DOUBLE_EQ(-gamma, rho2[3].real());
  EXPECT_DOUBLE_EQ(0.0, rho2[3].imag());
  EXPECT_DOUBLE_EQ(-0.5 * gamma, rho2[1].real());
  EXPECT_DOUBLE_EQ(0.0, rho2[1].imag());
  EXPECT_DOUBLE_EQ(-0.5 * gamma, rho2[2].real());
  EXPECT_DOUBLE_EQ(0.0, rho2[2].imag());
}

TEST(Decay, EToGTrace) {
  Decay d(0, 1, 0.2);
  CheckLindbladTraceProperty(d, 2);
}

TEST(Decay, EToGHermiticity) {
  Decay d(0, 1, 0.2);
  CheckLindbladHermiticityProperty(d, 2);
}
