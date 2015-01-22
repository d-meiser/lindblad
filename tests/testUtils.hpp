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
#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include <gtest/gtest.h>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <MasterEqn.hpp>

template <typename T>
void MY_EXPECT_EQ(T a, T b, size_t i) {
  EXPECT_EQ(a, b) << " i == " << i;
}

template <>
void MY_EXPECT_EQ(Amplitude a, Amplitude b, size_t i) {
  EXPECT_FLOAT_EQ(a.real(), b.real()) << "real, i == " << i;
  EXPECT_FLOAT_EQ(a.imag(), b.imag()) << "imag, i == " << i;
}

template <typename IterA, typename IterB>
void EXPECT_RANGES_EQ(IterA beginA, IterA endA, IterB beginB, IterB endB) {
  size_t i = 0;
  while (beginA != endA && beginB != endB) {
    MY_EXPECT_EQ(*beginA, *beginB, i);
    ++beginA;
    ++beginB;
    ++i;
  }
}

Amplitude RandomNumber() {
  return Amplitude(std::rand() / (double)RAND_MAX,
                   std::rand() / (double)RAND_MAX);
}

std::vector<Amplitude> randomVector(int dim) {
  std::vector<Amplitude> A(dim);
  std::generate(A.begin(), A.end(), RandomNumber);
  return A;
}

std::vector<Amplitude> hermitianMatrix(int dim) {
  std::vector<Amplitude> A = randomVector(dim * dim);
  std::vector<Amplitude> B(dim * dim, 0);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      B[i * dim + j] = A[i * dim + j] + conj(A[j * dim + i]);
    }
  }
  return B;
}

template <typename T>
void CheckLindbladTraceProperty(const T& op, int dim) {
  std::vector<Amplitude> A = hermitianMatrix(dim);
  std::vector<Amplitude> B(dim * dim, 0);
  op.apply(dim, &A[0], &B[0]);
  Amplitude trace;
  for (int i = 0; i < dim; ++i) {
    trace += B[i + i * dim];
  }
  EXPECT_TRUE(abs(trace) < 1.0e-12);
}

template <>
void CheckLindbladTraceProperty(const MasterEqn& op, int dim) {
  std::vector<Amplitude> A = hermitianMatrix(dim);
  std::vector<Amplitude> B(dim * dim, 0);
  op.apply(&A[0], &B[0]);
  Amplitude trace;
  for (int i = 0; i < dim; ++i) {
    trace += B[i + i * dim];
  }
  EXPECT_TRUE(abs(trace) < 1.0e-12);
}

template <typename T>
void CheckLindbladHermiticityProperty(const T& op, int dim) {
  std::vector<Amplitude> A = hermitianMatrix(dim);
  std::vector<Amplitude> B(dim * dim, 0);
  op.apply(dim, &A[0], &B[0]);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j <= i; ++j) {
      Amplitude bij = B[i * dim + j];
      Amplitude bji = B[j * dim + i];
      EXPECT_FLOAT_EQ(bij.real(), bji.real()) << "(i,j) == (" << i << "," << j
                                              << ") [real]";
      EXPECT_FLOAT_EQ(bij.imag(), -bji.imag()) << "(i,j) == (" << i << "," << j
                                              << ") [imag]";
    }
  }
}

template <>
void CheckLindbladHermiticityProperty(const MasterEqn& op, int dim) {
  std::vector<Amplitude> A = hermitianMatrix(dim);
  std::vector<Amplitude> B(dim * dim, 0);
  op.apply(&A[0], &B[0]);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j <= i; ++j) {
      Amplitude bij = B[i * dim + j];
      Amplitude bji = B[j * dim + i];
      EXPECT_FLOAT_EQ(bij.real(), bji.real()) << "(i,j) == (" << i << "," << j
                                              << ") [real]";
      EXPECT_FLOAT_EQ(bij.imag(), -bji.imag()) << "(i,j) == (" << i << "," << j
                                              << ") [imag]";
    }
  }
}

#endif
