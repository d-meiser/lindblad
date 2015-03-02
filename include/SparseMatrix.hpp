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

namespace Lindblad {

/**
@brief A non-zero entry in a sparse matrix
*/
struct SparseMatrixEntry {
/// @brief Constructor
  SparseMatrixEntry(int row, int col, Amplitude element)
      : row(row), col(col), element(element) {}
  int row;           /**< Row index of the entry */
  int col;           /**< Column index of the entry */
  Amplitude element; /**< Value of the entry */
};

/**
@brief Simplistic representation of sparse matrices.

Internally this sparse matrix type uses a coordinate format to represent
the non-zero entries of the matrix.
*/
struct SparseMatrix {
/**
@brief Build an emtpy matrix
*/
  SparseMatrix() : entries() {}

/**
@brief Build a matrix with one non-zero entry
*/
  SparseMatrix(SparseMatrixEntry entry) : entries(1, entry) {}

/**
@brief Add a non-zero entry to the sparse matrix
*/
  void add(SparseMatrixEntry entry) { entries.push_back(entry); }

/**
@brief Collection of non-zero entries
*/
  std::vector<SparseMatrixEntry> entries;
};

}

#endif
