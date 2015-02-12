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
#include <Coupling.hpp>
using namespace Lindblad;
using namespace Lindblad::Detail;
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

TEST(Coupling, Diagonal) {
  Coupling a(1, 0, 0.0);
  EXPECT_FALSE(a.isDiagonal());
  Coupling b(0, 0, 0.0);
  EXPECT_TRUE(b.isDiagonal());
}

TEST(Coupling, DiagonalEntryHasSameRowAsCol) {
  Coupling a(0, 0, 1.0);
  ASSERT_TRUE(a.isDiagonal());
  EXPECT_EQ(a.getRow(), a.getCol());
}
