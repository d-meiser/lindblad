/*
Copyright 2014 Dominic Meiser

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
#include <SparseMatrix.hpp>
using namespace Lindblad;

TEST(SparseMatrix, DefaultConstructible) {
  SparseMatrix a;
  EXPECT_TRUE(a.entries.empty());
}

TEST(SparseMatrix, ConstructFromMatrixEntry) {
  SparseMatrix a(SparseMatrixEntry(0, 1, 2.0));
  EXPECT_EQ(1u, a.entries.size());
}

TEST(SparseMatrix, CanAddTo) {
  SparseMatrix a;
  a.add(SparseMatrixEntry(2, 3, 5.0));
  EXPECT_EQ(1u, a.entries.size());
}
