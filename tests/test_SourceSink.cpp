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
#include <SourceSink.hpp>
using namespace Lindblad;
using namespace Lindblad::Detail;
#include <vector>
#include <cmath>
#include <testUtils.hpp>

TEST(SourceSink, Constructor) {
  SourceSink ss(1, 1.0);
}

TEST(SourceSink, WithZeroRateIsIdentity) {
  SourceSink ss(0, 0.0);
  std::vector<Amplitude> rho1(4, 1);
  std::vector<Amplitude> rho2(4, 0);
  ss.apply(2, &rho1[0], &rho2[0]);
  std::vector<Amplitude> expected(4, 0);
  EXPECT_RANGES_EQ(expected.begin(), expected.end(), rho2.begin(), rho2.end());
}

TEST(SourceSink, AppliedToScalarSystemIsIdentity) {
  SourceSink ss(0, 1.0);
  std::vector<Amplitude> rho1(1, 1);
  std::vector<Amplitude> rho2(1, 0);
  ss.apply(1, &rho1[0], &rho2[0]);
  Amplitude expected[] = {0};
  EXPECT_RANGES_EQ(expected, expected + 1, rho2.begin(), rho2.end());
}

TEST(SourceSink, IsTraceLess) {
  SourceSink ss(1, 2.6);
  CheckLindbladTraceProperty(ss, 3);
}

TEST(SourceSink, IsHermitian) {
  SourceSink ss(1, 2.6);
  CheckLindbladHermiticityProperty(ss, 3);
}

TEST(SourceSink, LeadsToNetGainInInputState) {
  int m = 1;
  double gamma = 2.6;
  SourceSink ss(m, gamma);
  int dim = 3;
  std::vector<Amplitude> rhoIn = hermitianMatrix(dim);
  std::vector<Amplitude> rhoOut(dim * dim, 0);
  ss.apply(dim, &rhoIn[0], &rhoOut[0]);
  EXPECT_GT(std::abs(rhoOut[m + m * dim]), 0);
}

TEST(SourceSink, LeadsToDecayingCoherence) {
  int m = 1;
  double gamma = 2.6;
  SourceSink ss(m, gamma);
  int dim = 3;
  std::vector<Amplitude> rhoIn = hermitianMatrix(dim);
  std::vector<Amplitude> rhoOut(dim * dim, 0);
  ss.apply(dim, &rhoIn[0], &rhoOut[0]);
  EXPECT_DOUBLE_EQ(-gamma,
                  std::real(rhoOut[0 + 1 * dim] / rhoIn[0 + 1 * dim]));
}
