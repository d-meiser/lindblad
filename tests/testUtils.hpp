#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include <gtest/gtest.h>
#include <cstdlib>
#include <algorithm>
#include <vector>

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

template <typename T>
void LindbladTraceProperty(T op, int dim) {
  std::vector<Amplitude> A = randomVector(dim * dim);
  std::vector<Amplitude> B(dim * dim, 0);
  op.apply(dim, &A[0], &B[0]);
  Amplitude trace;
  for (int i = 0; i < dim; ++i) {
    trace += B[i + i * dim];
  }
  EXPECT_TRUE(abs(trace) < 1.0e-12);
}

template <typename T>
void LindbladHermiticityProperty(T op, int dim) {
  std::vector<Amplitude> A = randomVector(dim * dim);
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

#endif
