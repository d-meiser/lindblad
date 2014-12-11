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
#include <GMRES.hpp>
#include <cmath>
#include <algorithm>

void printVector(int dim, Amplitude* v) {
  for (int i = 0; i < dim; ++i) {
    std::cout << v[i] << std::endl;
  }
}

struct GenRandAmplitude {
public:
  Amplitude operator()() {
    Amplitude a(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);
    return a;
  }
};


TEST(GMRES, Constructor) {
  GMRES *gmres = new GMRES(1); 
  ASSERT_TRUE(gmres != 0);
  delete gmres;
}

static void fDouble(int dim, const Amplitude *x, Amplitude *result,
                    void *ctx) {
  for (int i = 0; i < dim; ++i) {
    result[i] = 2.0 * x[i];
  }
}

TEST(GMRES, axpy) {
  GMRES gmres(2); 
  std::vector<Amplitude> y(2);
  y[0] = 2.0;
  y[1] = 3.0;
  std::vector<Amplitude> x(2);
  x[0] = -1.0;
  x[1] = 2.7;
  std::vector<Amplitude> result(2);
  double alpha = -2.4;
  gmres.axpy(alpha, &fDouble, &x[0], &y[0], &result[0], 0);
  EXPECT_FLOAT_EQ((alpha * 2.0 * x[0] + y[0]).real(), result[0].real());
}

TEST(GMRES, norm) {
  GMRES gmres(2); 
  std::vector<Amplitude> x(2, 1.3);
  double nrm = gmres.norm(2, &x[0]);
  EXPECT_FLOAT_EQ(sqrt(2.0 * 1.3 * 1.3), nrm);
}

TEST(GMRES, dot) {
  GMRES gmres(2); 
  std::vector<Amplitude> x(2, 1.3);
  double nrm = gmres.norm(2, &x[0]);
  Amplitude dot = gmres.dot(x.size(), &x[0], &x[0]);
  EXPECT_FLOAT_EQ(nrm * nrm, dot.real());
  EXPECT_FLOAT_EQ(0.0, dot.imag());
}

TEST(GMRES, setAbsTol) {
  GMRES gmres(2);
  double tol = 1.0e-6;
  gmres.setAbsTol(tol);
  EXPECT_FLOAT_EQ(tol, gmres.getAbsTol());
}

TEST(GMRES, smallEnoughTrue) {
  GMRES gmres(2);
  double tol = 1.0e-6;
  gmres.setAbsTol(tol);
  EXPECT_TRUE(gmres.smallEnough(tol / 2.0));
}

TEST(GMRES, smallEnoughFalse) {
  GMRES gmres(2);
  double tol = 1.0e-6;
  gmres.setAbsTol(tol);
  EXPECT_FALSE(gmres.smallEnough(tol * 2.0));
}

TEST(GMRES, setKrylovDim) {
  GMRES gmres(5);
  gmres.setKrylovDim(12);
  EXPECT_EQ(12, gmres.getKrylovDim());
}

TEST(GMRES, solveDouble) {
  int dim = 2;
  GMRES gmres(dim); 
  std::vector<Amplitude> b(dim);
  std::generate_n(b.begin(), dim, GenRandAmplitude());
  std::vector<Amplitude> x(dim);
  gmres.solve(&fDouble, &b[0], &x[0], 0);
  for (int i = 0; i < dim; ++i) {
    EXPECT_FLOAT_EQ(b[i].real() / 2.0, x[i].real()) << "i == " << i;
    EXPECT_FLOAT_EQ(b[i].imag() / 2.0, x[i].imag()) << "i == " << i;
  }
}

struct DiagOpCtx {
  double *diagonal;
  int dim;
};

static void diagonalOperator(int dim, const Amplitude *x, Amplitude *result,
                             void *ctx) {
  struct DiagOpCtx *cont = (struct DiagOpCtx*)ctx;
  assert(dim == cont->dim);
  for (int i = 0; i < dim; ++i) {
    result[i] = cont->diagonal[i] * x[i];
  }
}

TEST(GMRES, solveDiagonal) {
  int dim = 5;
  GMRES gmres(dim); 
  std::vector<Amplitude> b(dim);
  std::generate_n(b.begin(), dim, GenRandAmplitude());
  std::vector<Amplitude> x(dim);
  double diagonal[] = {1.0, 2.0, 3.0, 3.0, 3.0};
  struct DiagOpCtx ctx;
  ctx.diagonal = diagonal;
  ctx.dim = dim;
  gmres.solve(&diagonalOperator, &b[0], &x[0], &ctx);
  for (int i = 0; i < dim; ++i) {
    EXPECT_FLOAT_EQ(b[i].real() / diagonal[i], x[i].real()) << "i == " << i;
    EXPECT_FLOAT_EQ(b[i].imag() / diagonal[i], x[i].imag()) << "i == " << i;
  }
}
