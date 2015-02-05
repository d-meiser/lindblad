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
#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

#include <Amplitude.hpp>
#include <vector>

struct SparseMatrixEntry {
  SparseMatrixEntry(int row, int col, Amplitude element)
      : row(row), col(col), element(element) {}
  int row;
  int col;
  Amplitude element;
};

struct SparseMatrix {
  SparseMatrix() : entries() {}
  SparseMatrix(SparseMatrixEntry entry) : entries(1, entry) {}
  void add(SparseMatrixEntry entry) { entries.push_back(entry); }
  std::vector<SparseMatrixEntry> entries;
};

#endif


