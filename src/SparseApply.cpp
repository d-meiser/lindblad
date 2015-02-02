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
#include <cstdlib>
#include <iostream>
#include <detail/SparseApply.hpp>
#include <detail/Memory.hpp>
#include <config.h>

#if !defined(LINDBLAD_STACK_ARRAY_ALIGNMENT)
#define LINDBLAD_STACK_ARRAY_ALIGNMENT 64
#endif
#if !defined(LINDBLAD_SMALL_DIM)
#define LINDBLAD_SMALL_DIM 1024
#endif


static void extractStrided(const Amplitude *x, int offset, int stride, int n,
                           Amplitude *y);
static void insertStrided(const Amplitude *x, int offset, int stride, int n,
                          Amplitude *y);
static void rightApplyLargeDim(int row, int col, Amplitude alpha, int dim,
                               const Amplitude *A, Amplitude *B);

void leftApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
               Amplitude *B) {
  B += row * dim;
  A += col * dim;
  for (int c = 0; c < dim; ++c) {
    B[c] += alpha * A[c];
  }
}

void rightApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
                Amplitude *B) {
  if (dim > LINDBLAD_SMALL_DIM) {
    rightApplyLargeDim(row, col, alpha, dim, A, B);
  }
#if 0
  Amplitude *Acolumn = (Amplitude *)alloca(
      dim * sizeof(*Acolumn));
  Amplitude *Bcolumn = (Amplitude *)alloca(
      dim * sizeof(*Bcolumn));
#endif
  Amplitude *Acolumn = (Amplitude *)LINDBLAD_ALIGNED_ALLOCA(
      dim * sizeof(*Acolumn), LINDBLAD_STACK_ARRAY_ALIGNMENT);
  Amplitude *Bcolumn = (Amplitude *)LINDBLAD_ALIGNED_ALLOCA(
      dim * sizeof(*Bcolumn), LINDBLAD_STACK_ARRAY_ALIGNMENT);
  std::cout << Acolumn << std::endl;
  std::cout << Bcolumn << std::endl;
  extractStrided(A, row, dim, dim, &Acolumn[0]);
  extractStrided(B, col, dim, dim, &Bcolumn[0]);
#ifdef HAVE_PRAGMA_OMP_SIMD
#pragma omp simd aligned(Acolumn, Bcolumn: LINDBLAD_STACK_ARRAY_ALIGNMENT)
#endif
  for (int r = 0; r < dim; ++r) {
    Bcolumn[r] += alpha * Acolumn[r];
  }
  insertStrided(&Bcolumn[0], col, dim, dim, B);
}

void rightApplyLargeDim(int row, int col, Amplitude alpha, int dim,
                        const Amplitude *A, Amplitude *B) {
  B += col;
  A += row;
  for (int r = 0; r < dim; ++r) {
    B[r * dim] += alpha * A[r * dim];
  }
}


void extractStrided(const Amplitude *x, int offset, int stride, int n,
                           Amplitude *y) {
#ifdef HAVE_PRAGMA_OMP_SIMD
#pragma omp simd aligned(y : LINDBLAD_STACK_ARRAY_ALIGNMENT)
#endif
  for (int i = 0; i < n; ++i) {
    y[i] = x[offset + i * stride];
  }
}

void insertStrided(const Amplitude *x, int offset, int stride, int n,
                              Amplitude *y) {
#ifdef HAVE_PRAGMA_OMP_SIMD
#pragma omp simd aligned(x : LINDBLAD_STACK_ARRAY_ALIGNMENT)
#endif
  for (int i = 0; i < n; ++i) {
    y[offset + i * stride] = x[i];
  }
}

