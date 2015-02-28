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
#include <MasterEqn.hpp>
using namespace Lindblad;
#include <testUtils.hpp>

TEST(MasterEqn, Constructor) {
  MasterEqn* meqn = new MasterEqn(2);
  meqn->addCoupling(1, 0, 1.0);
  EXPECT_TRUE(meqn != 0);
  delete meqn;
}

TEST(MasterEqn, AddCoupling) {
  MasterEqn meqn(4);
  EXPECT_NO_THROW(meqn.addCoupling(2, 3, Amplitude(1.0, 11.0)));
}

TEST(MasterEqn, AddDecay) {
  MasterEqn meqn(4);
  EXPECT_NO_THROW(meqn.addDecay(2, 3, 2.0));
}

TEST(MasterEqn, AddSourceSink) {
  MasterEqn meqn(4);
  EXPECT_NO_THROW(meqn.addSourceSink(1, 2.0));
}

MasterEqn buildME(int dim) {
  MasterEqn meqn(dim);
  meqn.addCoupling(2, 1, Amplitude(2.0, -1.0));
  meqn.addCoupling(1, 1, Amplitude(1.0, -7.0));
  meqn.addCoupling(2, 0, Amplitude(3.0, -2.0));
  meqn.addCoupling(3, 4, Amplitude(4.0, -3.0));
  meqn.addCoupling(4, 2, Amplitude(7.0, -3.0));
  meqn.addDecay(0, 0, 2.0);
  meqn.addDecay(2, 3, 1.0);
  meqn.addDecay(2, 1, 3.0);
  meqn.addDecay(5, 7, 4.0);
  return meqn;
}

TEST(MasterEqn, TraceProperty) {
  MasterEqn meqn = buildME(8);
  CheckLindbladTraceProperty(meqn, 8);
}

TEST(MasterEqn, HermiticityProperty) {
  MasterEqn meqn = buildME(8);
  CheckLindbladHermiticityProperty(meqn, 8);
}

TEST(MasterEqn, NullOperatorHasZeroLevelShifts) {
  MasterEqn nullME(3);
  Amplitude omegas[3];
  nullME.getEnergyLevels(omegas);
  Amplitude zero(0);
  MY_EXPECT_EQ(zero, omegas[0], 0);
  MY_EXPECT_EQ(zero, omegas[1], 1);
  MY_EXPECT_EQ(zero, omegas[2], 2);
}

TEST(MasterEqn, DiagonalCouplingLeadsToLevelShift) {
  MasterEqn me(3);
  me.addCoupling(0, 0, 3.0);
  Amplitude omegas[3];
  me.getEnergyLevels(omegas);
  MY_EXPECT_EQ(Amplitude(3.0), omegas[0], 0);
  MY_EXPECT_EQ(Amplitude(0), omegas[1], 1);
  MY_EXPECT_EQ(Amplitude(0), omegas[2], 2);
}

TEST(MasterEqn, OffDiagonalCouplingDoesNotLeadToLevelShift) {
  MasterEqn me(3);
  me.addCoupling(0, 1, 3.0);
  Amplitude omegas[3];
  me.getEnergyLevels(omegas);
  MY_EXPECT_EQ(Amplitude(0), omegas[0], 0);
  MY_EXPECT_EQ(Amplitude(0), omegas[1], 1);
  MY_EXPECT_EQ(Amplitude(0), omegas[2], 2);
}

TEST(MasterEqn, BuildMatrix) {
  MasterEqn me(3);
  std::vector<Amplitude> matrix(9 * 9, 3.0);
  me.buildMatrix(&matrix[0]);
  MY_EXPECT_EQ(Amplitude(0), matrix[0], 0);
}

static int delta(int m, int n) { return (m == n) ? 1 : 0; }

TEST(MasterEqn, BuildMatrixForALevelShift) {
  MasterEqn me(2);
  double g = 2.3;
  me.addCoupling(0, 0, g);
  std::vector<Amplitude> matrix(4 * 4, 3.0);
  me.buildMatrix(&matrix[0]);
  int i = 0;
  for (int m = 0; m < 2; ++m) {
    for (int n = 0; n < 2; ++n) {
      for (int mp = 0; mp < 2; ++mp) {
        for (int np = 0; np < 2; ++np) {
          Amplitude expectedResult =
              g * Amplitude(0, -1) * Amplitude((delta(m, 0) - delta(n, 0)) *
                                               delta(m, mp) * delta(n, np));
          MY_EXPECT_EQ(expectedResult, matrix[i], i);
          ++i;
        }
      }
    }
  }
}

TEST(MasterEqn, BuildTransposedMatrix) {
  MasterEqn me(3);
  std::vector<Amplitude> matrix(9 * 9, 3.0);
  me.buildTransposedMatrix(&matrix[0]);
  MY_EXPECT_EQ(Amplitude(0), matrix[0], 0);
}

TEST(MasterEqn, BuildTransposedMatrixForALevelShift) {
  MasterEqn me(2);
  double g = 2.3;
  me.addCoupling(0, 0, g);
  std::vector<Amplitude> matrix(4 * 4, 3.0);
  me.buildTransposedMatrix(&matrix[0]);
  int i = 0;
  for (int m = 0; m < 2; ++m) {
    for (int n = 0; n < 2; ++n) {
      for (int mp = 0; mp < 2; ++mp) {
        for (int np = 0; np < 2; ++np) {
          Amplitude expectedResult =
              g * Amplitude(0, -1) * Amplitude((delta(m, 0) - delta(n, 0)) *
                                               delta(m, mp) * delta(n, np));
          MY_EXPECT_EQ(expectedResult, matrix[i], i);
          ++i;
        }
      }
    }
  }
}
