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
