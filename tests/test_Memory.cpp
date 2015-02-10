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
#include <Memory.hpp>

TEST(AlignedAlloca, ReturnsNonNullPtr) {
  double *a = (double *)LINDBLAD_ALIGNED_ALLOCA(8 * sizeof(*a), 64);
  EXPECT_TRUE(a != 0);
}

TEST(AlignedAlloca, ReturnsSufficientlySizedMemory) {
  int len = 8;
  double *a = (double *)LINDBLAD_ALIGNED_ALLOCA(8 * sizeof(*a), 16);
  for (int i = 0; i < len; ++i) {
    a[i] = 2.8;
  }
  EXPECT_DOUBLE_EQ(a[len - 1], 2.8);
}

TEST(AlignedAlloca, ReturnsCorrectlyAlignedMemory) {
  double *a = (double *)LINDBLAD_ALIGNED_ALLOCA(8 * sizeof(*a), 32);
  EXPECT_EQ((long long)a % 32, 0);
}

TEST(AlignedAlloca, ReturnsCorrectlyAlignedMemoryForNonPowerTwoAlignment) {
  double *a = (double *)LINDBLAD_ALIGNED_ALLOCA(8 * sizeof(*a), 9);
  EXPECT_EQ((long long)a % 9, 0);
}

TEST(AlignedAlloca, ReturnsCorrectlyAlignedMemoryForPrimeAlignment) {
  double *a = (double *)LINDBLAD_ALIGNED_ALLOCA(8 * sizeof(*a), 17);
  EXPECT_EQ((long long)a % 17, 0);
}

TEST(AlignedAlloca, IsAbleToDoMultipleAllocations) {
  double *a = (double *)LINDBLAD_ALIGNED_ALLOCA(4 * sizeof(*a), 64);
  double *b = (double *)LINDBLAD_ALIGNED_ALLOCA(4 * sizeof(*b), 64);
  EXPECT_GT(a, b);
}
